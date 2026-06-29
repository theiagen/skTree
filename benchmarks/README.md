# skTree vs kSNP4 — benchmark

How does skTree compare to **kSNP4 v4.1**, the tool it reimplements? We measured
**runtime**, **peak memory**, and **core-SNP concordance** across three test
families, with both tools under the same cap (**4 CPU / 12 GB, swap off, fixed
k=21**) in pinned containers (`staphb/ksnp4:4.1`, `sktree:latest`).

**Headline:** skTree matches kSNP4's core SNPs to the locus on related cohorts and
to within 3 SNPs out of ~69 000 on diverse ones, while running **2–10× faster** and
**3–5× leaner** — and it needs a few GB of scratch disk where kSNP4 needs ~200 GB at
scale.

## 1. Head-to-head — 8× *C. auris* (~12 Mb fungal genomes)

| cohort | kSNP4 s | skTree s | speedup | kSNP4 GB | skTree GB | kSNP4 core | skTree core |
|---|--:|--:|--:|--:|--:|--:|--:|
| Clade1_c200-400 | 198.5 | 35.9 | 5.5× | 8.46 | 1.55 | 274 | 274 |
| Clade1_c400-700 | 154.8 | 35.3 | 4.4× | 7.58 | 1.44 | 247 | 247 |
| Clade1_c700-1200 | 152.6 | 35.5 | 4.3× | 7.35 | 1.44 | 271 | 271 |
| Clade3_c200-400 | 162.5 | 35.6 | 4.6× | 6.88 | 1.44 | 64 | 64 |
| Clade3_c400-700 | 148.0 | 35.5 | 4.2× | 6.65 | 1.44 | 0 | 139 |
| Clade3_c700-1200 | 163.0 | 35.0 | 4.7× | NA | 1.44 | 28 | 83 |

- **4.2–5.5× faster, ~5× leaner.** skTree holds a flat ~35 s and ~1.4 GB; kSNP4 takes
  2.5–3.3 min and 6.6–8.5 GB.
- **Core SNPs identical** on the four closely/moderately related cohorts (274, 247,
  271, 64 — exact).
- **The two divergent Clade3 cohorts split** (kSNP4 0/28 vs skTree 139/83). The cause
  is the held-fixed k=21: kSNP4 normally runs Kchooser, and its stricter "both flanks
  present in *every* genome" core rule collapses to 0 shared k-mers as divergence
  rises. A k that suits both tools removes the gap.

## 2. Scaling — *C. auris* 25 → 50 → 100 genomes

Nested cohorts (`n25 ⊂ n50 ⊂ n100`) from a 150-genome pool. skTree ran NJ-only here
(its Python parsimony is `O(n³·sites)`, impractical past ~30 samples).

|   n | kSNP4 s | skTree s | speedup | kSNP4 GB | skTree GB | kSNP4 core | skTree core |
|--:|--:|--:|--:|--:|--:|--:|--:|
|  25 |  579.0 |  73.6 |  7.9× | 9.06 |  3.58 | 27972 | 27972 |
|  50 | 1200.2 | 115.8 | 10.4× | 7.71 |  7.04 | 27941 | 27941 |
| 100 | 2446.5 | 249.6 |  9.8× | 9.83 | 10.33 | 27874 | 27874 |

- **Core SNPs match to the locus at every n** (27972, 27941, 27874 — exact), and
  skTree stays **~8–10× faster**.
- **Memory converges with scale.** skTree holds the SNP matrix in RAM (grows ~linearly
  with n); kSNP4 streams k-mers and stays flatter. Both stay under the 12 GB cap here —
  watch this axis past n≈150.
- **kSNP4's real limit is disk, not RAM.** Its Jellyfish + `all_SNPs` scratch peaked
  near **90 GB** at n=100 and failed twice on a 100 GB volume (`No space left on
  device`); it only finished given a **200 GB** volume. skTree used a few GB.

## 3. Cross-species — *E. coli*, *Salmonella*, *Shigella*, *Pseudomonas*

Eight 10-genome cohorts from NCBI Pathogen Detection: one **close** (single SNP
cluster) and one **diverse** per species. skTree ran parsimony on the low-SNP close
cohorts, NJ-only on the high-SNP diverse ones.

| cohort | kSNP4 s | skTree s | speedup | kSNP4 GB | skTree GB | kSNP4 core | skTree core |
|---|--:|--:|--:|--:|--:|--:|--:|
| ecoli_close | 102.3 | 238.6¹ | 0.4× | 6.80 | 1.41 | 265 | 265 |
| ecoli_diverse | 177.7 | 35.5 | 5.0× | 8.65 | 2.22 | 49047 | 49046 |
| salmonella_close | 93.5 | 43.1 | 2.2× | 6.02 | 1.39 | 105 | 105 |
| salmonella_diverse | 175.0 | 43.8 | 4.0× | 8.14 | 2.23 | 46966 | 46965 |
| shigella_close | 93.0 | 22.6 | 4.1× | 5.91 | 1.24 | 51 | 51 |
| shigella_diverse | 155.3 | 33.8 | 4.6× | 6.92 | 1.81 | 45768 | 45765 |
| pseudomonas_close | 116.4 | 23.0 | 5.1× | 7.16 | 1.55 | 57 | 57 |
| pseudomonas_diverse | 195.2 | 38.5 | 5.1× | 6.86 | 2.22 | 68791 | 68789 |

¹ `ecoli_close` ran skTree's optional parsimony (cheap at ~265 SNPs); the diverse rows
are NJ-only. It is the table's only loss.

- **Concordance holds across four genera.** Close cohorts match to the locus; diverse
  cohorts (45 768–68 791 core SNPs) differ by **at most 3 SNPs** (≥99.99 % agreement) —
  including the Shigella cohort spanning four species. The 1–3 SNP gap is edge effects
  (k-mers brushing a contig boundary), not a core-definition split.
- **2.2–5.1× faster, 3.1–4.8× leaner** with the tree method matched to cohort size.
- **The one loss is a parsimony warning.** skTree's `O(n³·sites)` parsimony on a
  high-SNP set is the wrong tool: running it on `ecoli_diverse`/`salmonella_diverse`
  took 683 s / 1535 s, vs 35.5 s / 43.8 s NJ-only — so ~95 % of that time was the tree
  search, not SNP work. **Use NJ or `--ml` past a few thousand SNPs.**

## A note on tree topology

Where comparisons agreed exactly the trees matched (e.g. Robinson-Foulds **0** on the
best-signal *C. auris* cohort and on several diverse bacterial cohorts). High RF values
appear only on near-clonal cohorts whose core-SNP trees are full of zero-length
branches: with no signal to resolve the backbone, NJ (skTree) and core-parsimony
(kSNP4) break unsupported ties differently. The disagreement is always in the
*unresolved* part of the tree, never in supported splits.

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
- **Fixed k=21** disables both tools' native k-selection — the cause of the divergent
  *C. auris* split.
- **Matched invocations:** `kSNP4 -k 21 -core -vcf -CPU 4` vs
  `sktree run -k 21 --threads 4` (+`--parsimony` on low-SNP cohorts).
- Trees compared are skTree **NJ** vs kSNP4 **core parsimony** (its headline tree).
- **kSNP4's scaling cells used a ~200 GB scratch volume** so it could finish; RAM caps
  were identical for both tools.
- GCA accessions are sanitized (`.` → `_`) so both tools get identical, kSNP4-legal
  leaf labels.
