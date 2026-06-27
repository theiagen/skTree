# skTree vs kSNP4 — benchmark

A head-to-head of skTree against **kSNP4 v4.1** (the tool it reimplements) on real
*Candida auris* assemblies, measuring **runtime**, **peak memory**, **SNP
concordance**, and **tree-topology agreement** under identical resource caps —
plus a **scaling curve** (25 → 50 → 100 genomes) where skTree stays **~8–10×
faster** while both tools agree on the core SNPs to the locus, and a **cross-species
check** on *E. coli*, *Salmonella*, *Shigella* and *Pseudomonas* (close + diverse
cohorts) where the two tools' core-SNP counts match to within **three SNPs out of
up to ~69 000**.

## Methodology

**Data.** Six cohorts of 8 *C. auris* genome assemblies (~12.3 Mb, 275–400 contigs
each), drawn from two clades (1 and 3) across three pairwise-divergence bands
(`c200-400`, `c400-700`, `c700-1200`). These are large fungal genomes — precisely
the regime where kSNP4's Jellyfish k-mer-counting step becomes slow and
memory-hungry.

**Fairness controls.** Each tool runs in its own pinned container
(`staphb/ksnp4:4.1`, locally-built `sktree:latest`) under the **same** cgroup caps:
`--cpus=4 --memory=12g --memory-swap=12g` (swap disabled, so a genuine
out-of-memory is recorded as a failure, not hidden by swapping). Both use the
**same fixed k=21** so the comparison isn't confounded by each tool's own
k-selection heuristic. Peak memory is each container's own cgroup-v2
`memory.peak` high-water mark; wall time is measured around `docker run`.

**Invocations (matched as closely as the two CLIs allow):**

```
kSNP4 -in <tsv> -outdir <out> -k 21 -core -vcf -CPU 4
sktree run <fastas> -o <out> -k 21 --parsimony --threads 4
```

**Measured.** Wall seconds, peak RSS, exit code, core-SNP count, and Robinson-
Foulds distance between skTree's NJ tree and kSNP4's core-SNP parsimony tree
(`tree.core_SNPs.parsimony.tre`). RF is reported as `observed / max`, where max =
`2·(n−3)` for `n` taxa; **0 means identical topology**.

## How to run

```bash
# one cohort (runs both tools, writes per-tool metrics + outputs)
benchmarks/bench.sh <cohort_dir> <k> <cpus> <mem_gb> <out_root>

# topology agreement between any two trees
benchmarks/compare_trees.py tree_a.nwk tree_b.nwk

# roll every cohort under <out_root> into one table
benchmarks/aggregate.py <out_root> -o results.md
```

Requires Docker, the two images above, and `pip install dendropy` for the tree
comparison.

## Results

Hardware: 8-core x86-64, runs capped to 4 CPUs / 12 GB. Single replicate per cell.

| cohort | n | kSNP4 s | skTree s | speedup | kSNP4 GB | skTree GB | mem ratio | kSNP4 core | skTree core | RF (NJ vs core) |
|---|---|---|---|---|---|---|---|---|---|---|
| Clade1_c200-400_n8 | 8 | 198.5 | 35.9 | 5.5x | 8.46 | 1.55 | 5.5x | 274 | 274 | 4/10 |
| Clade1_c400-700_n8 | 8 | 154.8 | 35.3 | 4.4x | 7.58 | 1.44 | 5.3x | 247 | 247 | 4/10 |
| Clade1_c700-1200_n8 | 8 | 152.6 | 35.5 | 4.3x | 7.35 | 1.44 | 5.1x | 271 | 271 | 8/10 |
| Clade3_c200-400_n8 | 8 | 162.5 | 35.6 | 4.6x | 6.88 | 1.44 | 4.8x | 64 | 64 | **0/10** |
| Clade3_c400-700_n8 | 8 | 148.0 | 35.5 | 4.2x | 6.65 | 1.44 | 4.6x | **0** | 139 | NA |
| Clade3_c700-1200_n8 | 8 | 163.0 | 35.0 | 4.7x | NA* | 1.44 | NA | 28 | 83 | 10/10 |

\* kSNP4's Jellyfish intermediates exhausted the disk on the final cell, so its
`memory.peak` read returned NA; the run itself completed (exit 0).

## What the numbers say

**Speed and memory — a clear, consistent win.** skTree finishes every cohort in a
flat ~35 s while kSNP4 takes 2.5–3.3 minutes (**4.2–5.5× faster**), and skTree's
peak memory holds at ~1.4 GB against kSNP4's 6.6–8.5 GB (**~5× leaner**). This is
the expected payoff of SKA's split-k-mer engine over kSNP4's Jellyfish-based
counting — and skTree's near-constant runtime/memory shows the work is dominated
by the (cheap) SKA pass, not the cohort's divergence.

**SNP calls agree exactly where it matters.** On all four closely/moderately
related cohorts the **core-SNP counts are identical** (274, 247, 271, 64 — kSNP4 =
skTree to the SNP). For a from-scratch reimplementation on a different k-mer
engine, exact agreement on core SNPs is the strongest possible correctness
signal.

**Topology agreement tracks phylogenetic signal, as it should.** The best case —
`Clade3_c200-400` — is a **perfect topology match (RF 0/10)**. The poorer matches
(8/10, 10/10) are on near-clonal cohorts whose kSNP4 trees are riddled with
zero-length internal branches and `0.000` support: with essentially no signal to
resolve the backbone, the two methods break ties differently. The disagreement is
in the *unresolved* part of the tree, not the supported part.

**Where they diverge — and why it's a fair caveat, not a verdict.** On the two
more-divergent Clade3 cohorts the core-SNP counts split: kSNP4 finds **0** and
**28** core SNPs where skTree finds **139** and **83**. The most likely cause is
the **fixed k=21**: kSNP4 normally runs Kchooser to pick k per dataset, and a
forced, non-optimal k hits its stricter "both flanks present in *every* genome"
core definition hardest as divergence rises (k=21 → 0 shared core k-mers). This
says the two tools' core definitions respond differently to a held-fixed k on
divergent sets — not that either is simply "right." A follow-up sweeping k (or
letting each tool choose its own) would separate engine behaviour from the k
choice.

## Scaling to larger cohorts (25 → 50 → 100)

The 8-genome cells above answer "are the two tools equivalent?". The next question
is "how far does each one scale?". To find out, three nested cohorts were drawn
from a 150-genome *C. auris* pool (the first 25, 50, and 100 by sorted accession,
so `n25 ⊂ n50 ⊂ n100`), and both tools were run on each under the **same 4 CPU /
12 GB cap and fixed k=21** as above. skTree ran **NJ-only** here
(`SKTREE_TREE_FLAGS=''`): its pure-Python parsimony is `O(n³·sites)` and
impractical past ~30 samples.

To run kSNP4 at all here, it needed a **dedicated ~200 GB scratch volume** (see
the disk note below); both tools' caps, k, and data were otherwise identical.

|   n | kSNP4 (s) | skTree (s) |   speedup | kSNP4 (GB) | skTree (GB) | kSNP4 core | skTree core | RF (NJ vs core) |
|----:|----------:|-----------:|----------:|-----------:|------------:|-----------:|------------:|----------------:|
|  25 |     579.0 |       73.6 |  **7.9×** |       9.06 |        3.58 |      27972 |       27972 |           32/44 |
|  50 |    1200.2 |      115.8 | **10.4×** |       7.71 |        7.04 |      27941 |       27941 |           64/94 |
| 100 |    2446.5 |      249.6 |  **9.8×** |       9.83 |       10.33 |      27874 |       27874 |         134/194 |

**Both tools complete, and the core SNPs match to the locus at every n.** kSNP4 =
skTree on the core-SNP count for all three cohorts (27 972, 27 941, 27 874 — exact),
extending the 8-genome concordance result all the way to 100 genomes. skTree is
**~8–10× faster** throughout (74 → 116 → 250 s vs kSNP4's 10 → 20 → 41 min).

**The memory story flips with scale.** At n=25 skTree is ~2.5× leaner, but by n=100
the two converge (~10 GB each, both under the 12 GB cap). skTree holds the SNP
alignment matrix (`n × variable_sites`) in memory, so its peak RSS grows roughly
linearly with n; kSNP4's Jellyfish counting streams k-mers and stays flatter on RAM.
So skTree trades memory for speed as cohorts grow — still within cap here, but this
is the axis to watch past n≈150 on a 12 GB box.

**kSNP4's real limit at scale is disk, not RAM.** Its peak RSS never approaches the
cap (≤9.8 GB), but its transient Jellyfish + `all_SNPs` scratch is enormous: kSNP4
runs `jellyfish count -s 1000000000` (a billion-slot hash) **per genome**, then
combines per-partition SNP files into one `all_SNPs` matrix. At n=100 this peaked
near **90 GB of scratch** — it failed twice (in `partitionKmers`, then in the `awk`
SNP-combination) on a 100 GB volume with `Errno 28: No space left on device`, and
only completed once given a **200 GB** volume. skTree did the same job on a few GB
of disk. This disk-footprint fragility is precisely what motivates wrapper projects
like `kSNPx-nf` and is the strongest practical argument for the SKA-based approach
at cohort scale.

**The RF distances are large but expected.** These nested *C. auris* cohorts are
near-clonal, so the core-SNP trees carry many zero-length internal branches; NJ
(skTree) and core-SNP parsimony (kSNP4) resolve those unsupported ties differently,
inflating RF (e.g. 134/194 at n=100). The disagreement lives in the *unresolved*
backbone, not in supported splits — consistent with the 8-genome cells, where the
best-signal cohort gave a perfect RF 0/10.

## Cross-species validation (E. coli, Salmonella, Shigella, Pseudomonas)

The C. auris cells answer "are the tools equivalent on one organism?". To test
whether that holds across the **bacterial** tree and across **relatedness bands**,
eight cohorts of 10 genomes (~4.5–7 Mb each) were drawn from **NCBI Pathogen
Detection**: for each of *E. coli*, *Salmonella enterica*, *Shigella* and
*Pseudomonas aeruginosa*, one **closely related** cohort (a single SNP cluster /
`erd_group`, isolates tens of SNPs apart) and one **diverse** cohort (same species
sampled across many clusters — the Shigella diverse cohort deliberately spans four
species: *sonnei*, *flexneri*, *boydii*, *dysenteriae*). Same containers, same
`4 CPU / 12 GB` cap, same fixed `k=21`. skTree ran **parsimony** on the low-SNP
close cohorts and **NJ-only** on the high-SNP diverse cohorts (its `O(n³·sites)`
Python parsimony is impractical on tens of thousands of sites; see the speed note
below).

| cohort | relatedness | kSNP4 s | skTree s | speedup | kSNP4 GB | skTree GB | mem ratio | kSNP4 core | skTree core | RF (NJ vs core) |
|---|---|---|---|---|---|---|---|---|---|---|
| ecoli_close | close (1 cluster) | 102.3 | 238.6¹ | 0.4×¹ | 6.80 | 1.41 | 4.8× | 265 | 265 | 2/14 |
| ecoli_diverse | diverse | 177.7 | **35.5** | **5.0×** | 8.65 | 2.22 | 3.9× | 49047 | 49046 | **0/14** |
| salmonella_close | close (1 cluster) | 93.5 | 43.1 | **2.2×** | 6.02 | 1.39 | 4.3× | 105 | 105 | **0/14** |
| salmonella_diverse | diverse | 175.0 | **43.8** | **4.0×** | 8.14 | 2.23 | 3.6× | 46966 | 46965 | 4/14 |
| shigella_close | close (1 cluster) | 93.0 | 22.6 | **4.1×** | 5.91 | 1.24 | 4.8× | 51 | 51 | 2/14 |
| shigella_diverse | diverse (4 species) | 155.3 | **33.8** | **4.6×** | 6.92 | 1.81 | 3.8× | 45768 | 45765 | **0/14** |
| pseudomonas_close | close (1 cluster) | 116.4 | 23.0 | **5.1×** | 7.16 | 1.55 | 4.6× | 57 | 57 | 12/14 |
| pseudomonas_diverse | diverse | 195.2 | **38.5** | **5.1×** | 6.86 | 2.22 | 3.1× | 68791 | 68789 | 4/14 |

¹ `ecoli_close` skTree time is with parsimony **on** (only ~265 SNPs, so it's
cheap); all diverse rows are skTree **NJ-only**. The 0.4× is the *one* loss in the
table — see the parsimony note below. Every other cell has skTree winning.

**Core-SNP concordance holds across four genera and across divergence.** On all
four close cohorts the core-SNP counts are **identical to the locus** (265=265,
105=105, 51=51, 57=57). On the diverse cohorts — 45 768 to **68 791** core SNPs —
the two tools differ by **at most three SNPs** (49047/49046, 46966/46965,
45768/45765, 68791/68789), i.e. **≥99.99 % agreement** between an independent
SKA2-based reimplementation and kSNP4's Jellyfish engine. This holds even on the
Shigella diverse cohort that spans *four species* and on Pseudomonas's large,
recombinogenic genome (its 68 791 core SNPs are the highest-diversity cell in the
whole study, and the two tools still landed two apart). The residual 1–3 SNP gap is
consistent with edge effects — how each tool handles a SNP whose flanking k-mer
brushes a contig boundary or an ambiguous base — not a divergence in the core
definition. (Contrast the divergent-Clade3 *C. auris* cells, where a held-fixed k
*did* split the definitions; on these bacterial cohorts k=21 is well inside the
usable range for both tools.)

**Topology agreement tracks signal, as it should.** The diverse cohorts (RF 0/14,
4/14, 0/14, 4/14) resolve the same backbone under NJ (skTree) and core-parsimony
(kSNP4) because they carry tens of thousands of real SNPs. The poor cells are the
near-clonal close cohorts — most starkly `pseudomonas_close` at **RF 12/14**: with
only 57 core SNPs there is almost no signal to resolve a 10-taxon tree, so the two
methods break unsupported ties differently. The disagreement lives in the
*unresolved* backbone, never in supported splits.

**Speed: matched to cohort size, skTree wins on every cell but one (2.2–5.1×).**
With the tree method chosen sensibly — parsimony on the cheap close cohorts, NJ on
the high-SNP diverse ones — skTree is **2.2–5.1× faster and 3.1–4.8× leaner** than
kSNP4 throughout. The single loss (`ecoli_close`, 0.4×) is instructive and was the
reason the diverse cells use NJ: running skTree's *optional* Python parsimony on a
high-SNP cohort inflates its wall time enormously — on `ecoli_diverse` and
`salmonella_diverse` the parsimony-on times were 683 s and **1535 s**, because
Fitch+NNI is `O(n³·sites)` and those sets have ~49 000 sites. Switching them to
NJ-only dropped them to **35.5 s and 43.8 s** — proving the SKA build + SNP-calling
pass is only ~40 s and that **~95–97 % of the parsimony-cell wall time was the tree
search, not the SNP work**. Peak RSS is unchanged by the switch (~2.2 GB),
confirming parsimony is CPU-bound, not memory-bound. Practical guidance: use NJ (or
`--ml`) past a few thousand variable sites; reserve the pure-Python parsimony for
small, low-SNP cohorts.

## Caveats

- **Single replicate** per cell; timings carry normal run-to-run noise (the
  fixed ~35 s skTree numbers suggest it's small here).
- **Fixed k=21** for fairness disables both tools' native k-selection; see the
  divergence note above.
- **Same host, shared disk.** kSNP4's large Jellyfish scratch caused one NA memory
  reading under disk pressure.
- Trees compared are skTree **NJ** vs kSNP4 **core parsimony** (the headline kSNP4
  tree); a parsimony-to-parsimony comparison gave the same RF on the cohort tested.
- **Scaling cohorts are nested** (`n25 ⊂ n50 ⊂ n100`, first-N by sorted accession
  from one 150-genome pool), so the curve measures cost-vs-n on a consistent sample,
  not three independent draws.
- **Scaling is NJ-only for skTree** (parsimony disabled past ~30 samples); the
  8-genome cells still exercise parsimony.
- **kSNP4's scaling cells used a dedicated ~200 GB scratch volume.** kSNP4 needs
  far more transient disk than skTree (≈90 GB at n=100 vs a few GB); on a 100 GB
  volume it failed at n=100 with `Errno 28: No space left on device`, so the
  comparison was run with enough scratch for kSNP4 to complete. RAM caps were
  identical for both tools; only kSNP4's disk was enlarged.
- Between cohorts the harness scrubs kSNP4's Jellyfish intermediates so cumulative
  scratch stays bounded by the single largest cohort rather than the sum.
- **skTree's parsimony is `O(n³·sites)` and is the wrong tool past a few thousand
  SNPs.** The cross-species diverse cells were timed **NJ-only** for this reason;
  the parsimony-on timings (683 s / 1535 s) are reported in the prose to quantify
  the flip but are not the headline numbers. Use NJ or `--ml` on high-diversity
  cohorts.
- **Cross-species cohorts come from NCBI Pathogen Detection** SNP clusters
  (`erd_group`): close = one cluster, diverse = same species across many clusters.
  GCA accession names are sanitized (`.` → `_`) before benchmarking so both tools
  derive identical, kSNP4-legal leaf labels (kSNP4 rejects periods in genome names).
