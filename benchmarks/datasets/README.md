# Benchmark datasets

The exact genome sets behind the tables in [`../README.md`](../README.md). These
are the **inputs** — the precise cohorts that were benchmarked — so the comparison
is reproducible against a known genome set rather than a one-off measurement. The
runtime/memory/core-SNP numbers themselves live inline in the parent README.

## Files

- **`SAMPLES.tsv`** — one row per genome staged into a cohort
  (`family`, `cohort`, `n`, `sample_id`). 303 rows across 17 cohorts. The
  `sample_id` is the FASTA leaf label; for the cross-species cohorts it is the
  GCA accession with `.` rewritten to `_` (kSNP4-legal labels).
- **`provenance/`** — for the cross-species family, the exact NCBI accessions and
  their source SNP clusters:
  - `<cohort>.txt` — the 10 GCA accessions per cohort, one per line. Re-fetch with
    `datasets download genome accession --inputfile <cohort>.txt --include genome`.
  - `CLUSTERS.txt` — the NCBI Pathogen Detection SNP clusters (PDS) the close
    cohorts were drawn from; diverse cohorts draw one isolate from each of 10
    distinct clusters.
  - `shigella_diverse.species` — the species label per isolate, showing the
    Shigella diverse cohort spans flexneri / sonnei / boydii / sp.

## Cohort families

| family | cohorts | source |
|---|---|---|
| `cauris_8genome` | 6 × 8 genomes | *C. auris* clade subsets |
| `cauris_scaling` | n25 ⊂ n50 ⊂ n100 | *C. auris* 150-genome pool |
| `xspecies` | 8 × 10 genomes | NCBI Pathogen Detection (E. coli / Salmonella / Shigella / Pseudomonas) |

## How it was produced

Each cohort was run through `../bench.sh <cohort> 21 4 12 <out>` — identical
4 CPU / 12 GB caps, swap off, fixed k=21, pinned images `staphb/ksnp4:4.1` and
`sktree:latest`. The C. auris scaling cells used a dedicated 500 GB scratch volume
so kSNP4 had room for its ~90 GB of Jellyfish intermediates at n=100.
