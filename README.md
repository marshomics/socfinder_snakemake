# SOCfinder-Snakemake

A Snakemake re-implementation of [SOCfinder](https://github.com/lauriebelch/SOCfinder)
built to run hundreds of thousands of genomes on an SGE/PBS cluster. It produces the
same social-gene calls as upstream SOCfinder but processes genomes the way the problem
actually wants to be processed: many at once, one core each, with checkpoint/resume so a
failure at genome 150,000 doesn't restart the run.

## Why the original is slow at 300k genomes

`SOC_mine.py` runs three heavy tools per genome (KOFAMscan, three BLASTp searches,
antiSMASH) and hardcodes 32 threads for each. Run one genome at a time and you get two
problems at once. The tools oversubscribe the node (KOFAMscan + three BLASTp jobs each
ask for 32 threads *while* antiSMASH runs and grabs every core by default), and genomes
are processed serially, so total wall time is roughly `300,000 × (one genome's time)`.
Neither BLAST against a 90 MB database nor antiSMASH gets 32× faster on 32 threads, so
those threads are mostly wasted contention.

Three smaller things also cost time and, more importantly, robustness:

- `SOC_mine.py` builds a gffutils SQLite database for every genome (`gffutils.create_db`)
  and then never uses it — it re-parses the GFF by hand right afterwards. That's a wasted
  disk write per genome, and on a shared filesystem with 300k genomes those writes hurt.
- The GFF cleanup uses two `DataFrame.iterrows()` passes plus a per-group loop. Fine once,
  slow 300k times.
- `SOC_parse.py` parses antiSMASH GenBank files with `exec()` to build `vec0, vec1, …`
  variables it reads back out of `locals()`. It's O(n²) and breaks on large genomes.

## What this pipeline changes

The work across genomes is embarrassingly parallel, so that's where the speed comes from.

Each genome runs single-threaded by default, and the cluster runs as many genomes
concurrently as you have cores. On a few hundred 16-core nodes that's thousands of genomes
in flight at once, versus one at a time. This alone is the dominant win.

KOFAMscan — the slowest step in the original — isn't run here at all. You've already
computed it, so its table is just an input column and gets parsed directly.

antiSMASH stays in its default configuration so the calls match upstream exactly. A faster
`minimal` mode is available and described below, with its caveat.

The GFF fix is vectorised and drops the unused SQLite build. The antiSMASH GenBank parser
is rewritten with Biopython instead of `exec()`. Both produce identical gene sets to the
original (there's a test that proves the GFF rewrite keeps exactly the same features), just
faster and without the failure modes.

Snakemake adds the part that matters most at this scale: resume. Batches that finished stay
finished, and within a re-run batch only the genomes still missing a `SOCKS.csv` get redone.

## How it scales: batching

A naive "one Snakemake job per genome" workflow would build a multi-million-node DAG and
flood the SGE queue with 300k jobs. Neither Snakemake nor SGE enjoys that. Instead, genomes
are processed in **batches**. One batch is one cluster job that runs `batch_size` genomes,
`cores_per_batch` at a time, inside a thread pool. With `batch_size: 200` that's 1,500
cluster jobs for 300k genomes instead of 300,000 — small enough for both the scheduler and
the DAG, while every core stays busy.

```
samples.tsv ──► [batch 0] ─┐
                [batch 1] ─┤  each batch job runs run_one_genome.py per genome
                  ...      ├─►  results/genomes/<id>/SOCKS.csv ──► aggregate ──► master_SOCKS.tsv
                [batch N] ─┘
```

## Requirements

You need a working SOCfinder checkout — this pipeline drives the same tools, it doesn't
replace them. Follow the upstream README to:

1. create the `SOCfinder` conda environment (it already includes BLAST, antiSMASH deps,
   Biopython, pandas, numpy);
2. build the BLAST databases (`SOC_MakeBlastDB.py`), giving you `blast_databases/`;
3. install antiSMASH and its databases via the upstream `helper_script` so `antismash` is
   on your `PATH`.

Then add Snakemake to that environment:

```bash
conda activate SOCfinder
pip install snakemake               # Snakemake 8.x recommended
pip install snakemake-executor-plugin-cluster-generic   # for SGE/PBS submission
```

You also need your precomputed KOFAMscan tables (see the next section).

## Input: the sample sheet

The pipeline reads a tab-separated `samples.tsv` with one row per genome:

| column  | meaning                                                            |
|---------|-------------------------------------------------------------------|
| `id`    | unique genome id (used as the output folder name)                 |
| `faa`   | protein FASTA                                                     |
| `fna`   | nucleotide FASTA                                                  |
| `gff`   | GFF3 annotation                                                  |
| `kofam` | precomputed KOFAMscan table (see below)                          |
| `gram`  | `p`, `n`, or `both` (optional; falls back to `gram_default`)      |

`config/samples.tsv.example` shows the format. To build one from directories of inputs:

```bash
python workflow/scripts/make_samplesheet.py \
    --faa-dir proteins --fna-dir genomes --gff-dir gffs --kofam-dir kofam \
    --gram both --out config/samples.tsv
```

### KOFAMscan table format

Tables must be KOFAMscan's `detail-tsv` output, produced the way SOCfinder expects:

```bash
exec_annotation INPUT.faa -o GENOME.txt -f detail-tsv --threshold-scale 0.75
```

The `--threshold-scale 0.75` matters — SOCfinder relies on it, and the parser keeps only
rows KOFAMscan marked significant (a `*` in the first column). The gene-name column has to
use the same protein identifiers as the `.faa` headers and the GFF, otherwise the final
deduplication across modules won't line up. If your inputs all come from one NCBI Datasets
download per genome, they already match.

### Running KOFAMscan inside the pipeline

If you'd rather have the pipeline run KOFAMscan instead of supplying tables, set it in
`config/config.yaml`:

```yaml
run_kofamscan:  true
kofam_profiles: "/path/to/SOCfinder/KOFAM/profiles/prokaryote.hal"
kofam_ko_list:  "/path/to/SOCfinder/KOFAM/ko_list"
kofam_threshold_scale: 0.75   # leave at 0.75 to match SOCfinder
kofam_cpus:     1
```

With this on, the `kofam` column in `samples.tsv` is optional and `make_samplesheet.py`'s
`--kofam-dir` can be dropped. Each genome's table is written to
`results/genomes/<id>/kofam.txt` at the 0.75 scale SOCfinder uses, and it's reused rather
than recomputed if a genome is retried, so a failed batch doesn't repeat the expensive part.

One thing to weigh: KOFAMscan is the slowest step in SOCfinder, and it reloads the full KO
profile set for every genome. Turning this on puts that cost back into the run — which is
exactly why precomputing the tables once (as you've done) and feeding them in is the faster
route for 300k genomes. The toggle is here for when you actually want a fresh run. If you do,
give KOFAMscan more memory (raise `mem_per_core_mb` to ~6000) and expect it, not antiSMASH,
to dominate wall time.

## Running it

Point `config/config.yaml` at your SOCfinder checkout:

```yaml
blast_db_dir:    "/path/to/SOCfinder/blast_databases"
social_ko:       "/path/to/SOCfinder/inputs/SOCIAL_KO.csv"
antismash_types: "/path/to/SOCfinder/inputs/antismash_types.csv"
```

Check the plan without running anything:

```bash
snakemake -n
```

Run on one machine (say 32 cores):

```bash
snakemake --cores 32
```

Run on the cluster:

```bash
mkdir -p logs/sge
snakemake --profile profiles/sge
```

The SGE profile is in `profiles/sge/config.yaml`; it submits one `qsub` per batch with
`-pe smp {threads}` and `-l h_vmem`. PBS/Torque and LSF submit lines are in the comments
there. On **Snakemake 7** there's no executor plugin — submit directly instead:

```bash
snakemake --jobs 300 --latency-wait 60 --rerun-incomplete --restart-times 2 --keep-going \
  --cluster "qsub -cwd -V -pe smp {threads} -l h_vmem={resources.mem_mb}M -o logs/sge/ -e logs/sge/"
```

When it finishes you get `results/master_SOCKS.tsv` (every social gene, one row per
`genome_id, protein_id`), `results/master_summary.tsv` (per-genome counts), and
`results/missing.txt` (genomes with no output — empty on a clean run).

## Tuning for your cluster

The knobs are in `config/config.yaml`:

- `batch_size` — genomes per cluster job. Keep it well above `cores_per_batch` (200 is a
  good start). Bigger batches mean fewer scheduler jobs but coarser resume.
- `cores_per_batch` — cores each batch job requests, and how many genomes it runs at once.
  Match it to a sensible slot count on your nodes.
- `threads_per_genome` — leave at 1 for maximum throughput across 300k genomes. Raise it
  only if you have few genomes and want each finished sooner.
- `mem_per_core_mb` — antiSMASH needs roughly 2–4 GB per concurrent genome. The job asks
  for `cores_per_batch × mem_per_core_mb`.

### Disk

Keep `keep_intermediates: false` (the default). antiSMASH writes tens of MB per genome
across many files; at 300k genomes, keeping it means terabytes and hundreds of millions of
inodes. The pipeline extracts what SOCfinder needs and deletes the antiSMASH working
directory and BLAST outputs per genome, leaving only the small CSVs. Even so, budget for
300k small output folders — if your filesystem penalises inode count, raise `batch_size`
and tar finished batches.

### Scratch / filesystem load

`local_scratch: true` (the default) runs the heavy intermediates — BLAST output, the
antiSMASH working directory, and KOFAMscan's temp dir — on each node's local `$TMPDIR`,
and writes only the small result CSVs back to your shared output directory. With a couple
of thousand genomes running at once, all of them hammering antiSMASH temp files on
`/ebio/abt3_scratch` would make the shared filesystem the bottleneck; keeping that churn on
node-local disk avoids it. The node needs enough local scratch for `cores_per_batch`
concurrent antiSMASH runs (a few hundred MB each, so tens of GB — fine on these nodes). Set
`scratch_base` to override where local scratch lives, or `local_scratch: false` to work in
the output directory instead.

### Resume

Just re-run the same command. Snakemake skips completed batches; partial batches redo only
the genomes still missing a `SOCKS.csv`. `--restart-times 2` (in the profile) also auto-
retries a batch that fails, which absorbs transient node failures.

## antiSMASH fidelity (read before using `minimal`)

`antismash_mode: default` reproduces upstream SOCfinder exactly. With KOFAMscan precomputed,
antiSMASH is now the largest remaining per-genome cost, so it's tempting to switch on
`antismash_mode: minimal`, which skips antiSMASH's non-core analyses.

The catch: `minimal` can stop antiSMASH from labelling the non-core cluster genes
(`biosynthetic-additional`, `transport`, `regulatory`), and SOCfinder counts every gene
that carries a `gene_kind` label inside a social cluster. So `minimal` may report fewer
secondary-metabolite genes. Whether that changes your results depends on your genomes.

If you want the speedup, validate it first: run a few thousand genomes both ways and diff
the `secondary_metabolites` column of `master_summary.tsv`. If it's unchanged for your data,
use `minimal`; if not, stay on `default`.

## Checking against upstream

The two genomes shipped with SOCfinder are the reference: *B. aphidicola* has 9 social
genes, *P. salmonis* has 64. Run both through this pipeline (with their KOFAMscan tables)
and confirm `total` in `master_summary.tsv` matches before turning it loose on 300k.

Two tests run without any of the heavy tools installed:

```bash
# vectorised GFF fix keeps exactly the same features as the original code
python test/check_gff_equivalence.py /path/to/SOCfinder/test/B_aphidicola.gff \
                                     /path/to/SOCfinder/test2/P_salmonis.gff

# KOFAM / BLAST / antiSMASH parser ports give the expected calls on synthetic inputs
python test/test_parsers.py
```

## Layout

```
socfinder-snakemake/
├── config/
│   ├── config.yaml            # paths + scaling knobs
│   └── samples.tsv.example    # sample-sheet template
├── profiles/sge/config.yaml   # SGE/PBS submission profile (Snakemake 8)
├── workflow/
│   ├── Snakefile              # batching DAG + resume
│   └── scripts/
│       ├── socfinder_lib.py   # vectorised GFF fix + parser ports (no exec/SQLite)
│       ├── run_one_genome.py  # one genome: gff-fix, blast×3, antismash, parse, combine
│       ├── run_batch.py       # runs a batch of genomes concurrently
│       ├── aggregate.py       # per-genome outputs → master tables
│       └── make_samplesheet.py
└── test/
    ├── check_gff_equivalence.py
    └── test_parsers.py
```

## What this doesn't change

The science is untouched. Same BLAST databases, same e-value thresholds, same KO list, same
antiSMASH detection, same combination logic. If you publish with it, cite the SOCfinder
paper (Belcher et al. 2023, *Microbial Genomics*, doi:10.1099/mgen.0.001171).
