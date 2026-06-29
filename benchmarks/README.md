# skTree vs kSNP4 — benchmark

How does skTree compare to **kSNP4 v4.1**, the tool it reimplements? We measured
**runtime**, **peak memory**, and **core-SNP concordance** across three test
families, with both tools under the same cap (**4 CPU / 12 GB, swap off, fixed
k=21**) in pinned containers (`staphb/ksnp4:4.1`, `sktree:latest`).

**Headline:** across all 17 cohorts — fungal and bacterial, near-clonal to
cross-genus — skTree's core SNPs are **identical to kSNP4's, exactly**, while it
runs **2–8× faster** on neighbour-joining workloads and uses **3–5× less memory**
at small scale. The exceptions are honest and few: parsimony on a high-SNP cohort
is slow, and skTree's in-RAM matrix overtakes kSNP4's memory past ~50 genomes.

Every cohort below is pinned in [`datasets/`](datasets/) — exact sample IDs
([`SAMPLES.tsv`](datasets/SAMPLES.tsv), 303 genomes) and, for the cross-species
family, the NCBI accessions and source SNP clusters.

## 1. Head-to-head — 8× *C. auris* (~12 Mb fungal genomes)

| cohort | kSNP4 s | skTree s | speedup | kSNP4 GB | skTree GB | core SNPs |
|---|--:|--:|--:|--:|--:|--:|
| Clade1_c200-400 | 146.9 | 34.6 | 4.2× | 5.53 | 1.45 | 274 |
| Clade1_c400-700 | 150.6 | 36.1 | 4.2× | 5.10 | 1.44 | 247 |
| Clade1_c700-1200 | 166.5 | 36.0 | 4.6× | 5.08 | 1.45 | 271 |
| Clade3_c200-400 | 76.1 | 36.3 | 2.1× | 6.07 | 1.54 | 64 |
| Clade3_c400-700 | 152.3 | 34.6 | 4.4× | 5.58 | 1.45 | 139 |
| Clade3_c700-1200 | 147.1 | 34.1 | 4.3× | 5.13 | 1.44 | 83 |

- **Core SNPs identical on every cohort** (274, 247, 271, 64, 139, 83 — exact),
  including the divergent Clade3 cohorts. The `core SNPs` column is one value
  because both tools agreed to the locus.
- **2.1–4.6× faster, ~3.5× leaner.** skTree holds a flat ~35 s and ~1.4 GB; kSNP4
  takes 1.3–2.8 min and 5–6 GB.

## 2. Scaling — *C. auris* 25 → 50 → 100 genomes

Nested cohorts (`n25 ⊂ n50 ⊂ n100`) from a 150-genome pool. skTree ran NJ-only here
(its Python parsimony is `O(n³·sites)`, impractical past ~30 samples).

|   n | kSNP4 s | skTree s | speedup | kSNP4 GB | skTree GB | core SNPs |
|--:|--:|--:|--:|--:|--:|--:|
|  25 |  481.1 |  72.1 | 6.7× | 8.81 |  3.48 | 27972 |
|  50 |  964.7 | 117.9 | 8.2× | 9.55 |  7.68 | 27941 |
| 100 | 2010.3 | 251.4 | 8.0× | 9.42 | 10.08 | 27874 |

- **Core SNPs match to the locus at every n** (27972, 27941, 27874 — exact), and
  skTree stays **~7–8× faster**.
- **Memory crosses over by n≈50.** skTree holds the SNP matrix in RAM (grows
  ~linearly with n); kSNP4 streams k-mers and stays flat near 9–10 GB. They draw
  even at n=50 (7.7 vs 9.5 GB) and skTree edges *past* kSNP4 at n=100 (10.1 vs
  9.4 GB). Both still fit the 12 GB cap here — but skTree, not kSNP4, is the one to
  watch on RAM past n≈100.
- **kSNP4's real limit is disk, not RAM.** Its Jellyfish + `all_SNPs` scratch
  peaked near **90 GB** at n=100; this run used a dedicated 500 GB volume so it
  could finish. skTree used a few GB.

## 3. Cross-species — *E. coli*, *Salmonella*, *Shigella*, *Pseudomonas*

Eight 10-genome cohorts from NCBI Pathogen Detection: one **close** (single SNP
cluster) and one **diverse** per species. skTree ran parsimony on the close
cohorts, NJ-only on the high-SNP diverse ones. Exact accessions and source
clusters are in [`datasets/provenance/`](datasets/provenance/).

| cohort | kSNP4 s | skTree s | speedup | kSNP4 GB | skTree GB | core SNPs |
|---|--:|--:|--:|--:|--:|--:|
| ecoli_close | 99.1 | 180.6¹ | 0.5× | 6.86 | 1.35 | 151 |
| ecoli_diverse | 174.7 | 33.4 | 5.2× | 9.13 | 2.14 | 48939 |
| salmonella_close | 92.0 | 23.6 | 3.9× | 6.05 | 1.29 | 93 |
| salmonella_diverse | 158.8 | 30.3 | 5.2× | 7.68 | 2.14 | 55620 |
| shigella_close | 89.2 | 53.4¹ | 1.7× | 5.86 | 1.24 | 137 |
| shigella_diverse | 129.0 | 25.0 | 5.2× | 6.67 | 1.59 | 36697 |
| pseudomonas_close | 119.7 | 207.2¹ | 0.6× | 8.70 | 1.94 | 1064 |
| pseudomonas_diverse | 175.5 | 35.1 | 5.0× | 9.02 | 2.47 | 65382 |

¹ skTree's optional parsimony tree search dominates these close cohorts (it grows
`O(n³·sites)`); the slowdown tracks the SNP count — mild at 137 SNPs (shigella),
heavy at 1064 (pseudomonas). The diverse rows are NJ-only.

- **Concordance is exact across four genera.** Every cohort matches to the locus —
  close and diverse alike, from 93 core SNPs up to 65 382, including the Shigella
  cohort spanning four species. No edge-effect drift, no core-definition split.
- **~5× faster and 3–5× leaner on the diverse (NJ) cohorts;** the only losses are
  the parsimony close cohorts.
- **The lesson is tree-method selection, not SNP accuracy.** skTree's parsimony is
  the wrong tool whenever SNP counts climb — use NJ or `--ml` past a few thousand
  SNPs, and prefer NJ even on close cohorts once parsimony's tree search starts to
  dominate (≈1000 SNPs here).

## A note on tree topology

Where core SNPs agreed the trees matched (Robinson-Foulds **0** on the
best-signal cohorts). High RF values appear only on near-clonal cohorts whose
core-SNP trees are full of zero-length branches: with no signal to resolve the
backbone, NJ (skTree) and core-parsimony (kSNP4) break unsupported ties
differently. The disagreement is always in the *unresolved* part of the tree,
never in supported splits.

## How to run

```bash
# one cohort (runs both tools, writes per-tool metrics + outputs)
benchmarks/bench.sh <cohort_dir> <k> <cpus> <mem_gb> <out_root>

# topology agreement between two trees
benchmarks/compare_trees.py tree_a.nwk tree_b.nwk

# roll every cohort under <out_root> into one table
benchmarks/aggregate.py <out_root> -o results.md
```

Needs Docker, the two images above, and `pip install dendropy`.

## Caveats

- **Single replicate** per cell (skTree's flat ~35 s suggests low noise).
- **Fixed k=21** disables both tools' native k-selection (Kchooser), for a
  controlled comparison.
- **Matched invocations:** `kSNP4 -k 21 -core -vcf -CPU 4` vs
  `sktree run -k 21 --threads 4` (+`--parsimony` on the close/low-SNP cohorts).
- Trees compared are skTree **NJ** vs kSNP4 **core parsimony** (its headline tree).
- **kSNP4's scaling cells used a dedicated ~500 GB scratch volume** so it could
  finish; RAM caps were identical for both tools.
- GCA accessions are sanitized (`.` → `_`) so both tools get identical,
  kSNP4-legal leaf labels.
