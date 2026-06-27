# skTree — Research Synthesis

A reference-free / alignment-free SNP discovery and phylogenetics tool, built in
Python as a wrapper on top of the **SKA2** split-k-mer engine — it turns SKA2's
raw output into an end-to-end workflow.

> Status: research complete. 25 claims verified 3-0 (adversarial 3-voter), 0
> refuted, across 14 primary sources (kSNP4 MBE 2023; kSNP3 Bioinformatics 2015;
> SKA2 Genome Research 2024; original SKA preprint; ska.rust GitHub + docs.rs;
> skalo preprint; bacpop tree-building guide; OptiK preprint).

---

## 1. How kSNP4 works

- **Reference-free, alignment-free.** kSNP identifies SNPs and builds trees from
  multiple genomes *without* whole-genome alignment or a reference genome
  (true since kSNP3.0).
- **SNP-locus definition.** Uses **odd-length k-mers**. Two k-mers are
  *homologous* if identical at every position **except the central base**. The
  conserved flanks make the locus homologous; the central base is the SNP allele.
  Practically: build the merged cross-genome k-mer list, collapse each k-mer to a
  **canonical form** (min of forward / reverse-complement), then a SNP locus is a
  flank for which **central-base variants exist** across genomes. Intra-genome
  allele conflicts (same flank, two different central bases within one genome) are
  removed.
- **Inputs.** FASTA — complete genomes, assemblies, **or raw reads**.
- **SNP categories.**
  - **Core SNPs** — locus present in **all** genomes.
  - **Majority / pan SNPs** — present in a fraction (e.g. ≥ majority) of genomes.
  - (kSNP also reports "all SNPs".)
- **k selection (Kchooser4).** Increments odd k until **most k-mers are unique**
  in a representative genome — i.e. picks the smallest odd k at which k-mers are
  effectively unique (reported via FCK, the fraction of core k-mers). Too small →
  k-mers not unique (false homology); too large → SNPs lost / sensitive to error.
- **Trees.** Three methods: **parsimony** (its primary/default, via the SNP
  presence matrix), **neighbor-joining (NJ)**, and **maximum-likelihood (ML)**.

## 2. How SKA / SKA2 works

- **Split k-mer.** An odd-length k-mer with the **middle base removed/variable**:
  two equal-length flanks (pattern `XXXX·XXXX`) surrounding one wildcard center.
  The flanks are fixed reference points that locate the variable middle base (the
  SNP). This is *structurally the same idea as kSNP's SNP locus.*
- **Data structure / `.skf`.** A hash map: **split-k-mer flank → middle base**.
  The `.skf` file stores this dictionary; **SKA2's format is snappy-compressed and
  not backward-compatible with SKA1**. Merging two `.skf`s = matching identical
  flanks so middle bases line up (alignment is positional-by-matching, **unordered**
  — no genome coordinates).
- **Variant calling.** Exact match of flanks across strains; if the **same flanks
  carry a different middle base**, that's a variant. Canonical (strand) handling
  applied.
- **Subcommands.**
  - `ska build` — FASTA/FASTQ in → `.skf`.
  - `ska align` — **reference-free SNP alignment** (the core deliverable for trees).
  - `ska map` — map split k-mers onto a reference to get ordered/coordinate output.
  - `ska distance` — pairwise SNP distances + mismatch fraction.
  - `ska lo` (**skalo**) — newer **graph-based** variant caller traversing a
    coloured de Bruijn graph (handles some indels/clustered variants SKA's exact
    flank matching misses).
- **Scope / limitations (the critical constraint).**
  - Designed for **closely related, small haploid genomes** (bacteria, viruses).
  - **Recall drops beyond ~1% sequence divergence** — exact flank matching fails as
    flanks accumulate mutations.
  - **Cannot resolve two SNPs closer together than k** (their flanks overlap the
    other variant).
- **Trees.** SKA itself does not build trees: you take `ska align` output → feed to
  an external tree builder (IQ-TREE, RAxML, FastTree, etc.). No reference bias.

## 3. kSNP4 vs SKA — overlap and gaps

| Capability | kSNP4 | SKA2 | Gap to fill in skTree |
|---|---|---|---|
| Reference-free SNP discovery | ✅ | ✅ (`ska align`) | — (SKA provides core) |
| Odd-k central-base locus model | ✅ | ✅ (split k-mer) | conceptually identical |
| Raw reads input | ✅ | ✅ (`ska build` FASTQ) | — |
| Core SNP set | ✅ explicit | implicit (filter alignment) | **skTree computes** core vs majority partitions from the alignment |
| Majority/pan SNP set | ✅ | ❌ | **skTree adds** (missing-fraction thresholding) |
| Optimal-k selection | ✅ Kchooser4 | ❌ (user picks k) | **skTree adds** optimal-k selection |
| Handles higher divergence | ✅ (better) | ❌ (~1% ceiling) | **document boundary**; warn user |
| Parsimony tree | ✅ built-in | ❌ | **skTree adds** (or wraps) |
| NJ / ML tree | ✅ | external | **skTree wraps** external builders |
| SNP-vs-locus annotation, matrices | ✅ rich outputs | minimal | **skTree adds** SNP-vs-locus reports |

**Net:** SKA2 already delivers the *hardest* piece — fast reference-free SNP
discovery + alignment. skTree's job is the **surrounding workflow**:
optimal-k selection, core/majority SNP partitioning, the SNP matrix, and tree
inference (parsimony + NJ/ML), plus SNP-vs-locus reports — all in Python around the
`ska` binary.

## 4. Python ecosystem

- **SKA2 has no Python API** → wrap the **`ska` CLI as a subprocess** (decided).
  Parse `ska align` (FASTA alignment), `ska distance` (TSV), and optionally `.skf`
  metadata via `ska nk`/`ska weed`.
- **Sequence IO:** Biopython (`SeqIO`) or `pyfastx` (fast, indexed FASTA/FASTQ).
- **k selection:** can implement a uniqueness scan natively; **OptiK**
  (truncated-SVD on k-mer frequency + clustering metrics: Silhouette /
  Calinski-Harabasz / Davies-Bouldin) is a modern, statistics-driven alternative
  to Kchooser — a possible "advanced" mode.
- **Trees:**
  - **DendroPy** — parsimony scoring + NJ, tree IO/manipulation (pure Python).
  - **scikit-bio** — NJ from a distance matrix.
  - **ete3 / Bio.Phylo** — visualization / IO.
  - **External ML:** wrap **IQ-TREE** / **RAxML-NG** / **FastTree** as subprocess
    for ML and bootstrap.
- **Optional annotation:** **pyrodigal** (gene calling) to annotate whether SNPs
  fall in coding regions.
- **k-mer counting in Python is slow** — but we offload counting to SKA2, so Python
  only handles the (small) SNP alignment/matrix. No performance bottleneck.

## 5. Proposed architecture (data flow)

```
FASTA/FASTQ inputs
      │
      ▼
[ ska build ]  ──►  per-sample .skf            (engine: SKA2 subprocess)
      │
      ▼
[ ska merge ]  ──►  combined .skf
      │
      ├──► [ ska align ]   ──► reference-free SNP alignment (FASTA)
      └──► [ ska distance ] ──► pairwise SNP-distance TSV
                 │
                 ▼
        skTree core (Python)
          • k-selection (Kchooser-equiv / OptiK)   → choose odd k
          • parse alignment → SNP matrix (samples × loci)
          • partition: core SNPs / majority SNPs / all SNPs
          • optional pyrodigal annotation
                 │
                 ▼
        skTree phylo (Python)
          • parsimony tree (DendroPy)
          • NJ tree (scikit-bio/DendroPy from distances)
          • ML tree (IQ-TREE/RAxML wrapper, optional + bootstrap)
                 │
                 ▼
        Outputs: SNP alignment (.fasta), SNP matrix (.tsv),
                 core/majority SNP sets, trees (.nwk), summary report
```

### Module layout (proposed)
```
sktree/
  cli.py            # argparse/typer entrypoint, subcommands
  engine/ska.py     # subprocess wrapper for ska build/merge/align/distance
  kselect.py        # optimal odd-k selection (uniqueness scan; OptiK optional)
  snps.py           # parse alignment → matrix; core/majority/all partitions
  annotate.py       # optional pyrodigal SNP-in-gene annotation
  phylo.py          # parsimony / NJ / ML tree inference
  report.py         # SNP summary + matrices
  io.py             # FASTA/FASTQ handling, output writers
```

### Validation strategy
1. **Unit/TDD** on synthetic genomes with known, hand-placed SNPs → assert exact
   SNP loci, core/majority partition, matrix.
2. **Cross-check vs SKA** — our parsed alignment must equal `ska align` output.
3. **Cross-check vs kSNP4** — on a public closely-related bacterial set (e.g. a
   small *Salmonella*/*E. coli* outbreak set), compare core-SNP counts and tree
   topology (Robinson-Foulds distance) against kSNP4 and confirm within expected
   tolerance.
4. **Divergence boundary test** — confirm graceful warning as divergence → >1%.

## Sources
- kSNP4: Hall & Nisbet 2023, *Mol Biol Evol* 40(11):msad235.
- kSNP3.0: Gardner, Slezak & Hall 2015, *Bioinformatics* 31(17):2877.
- kSNP (v2): PMC3857212 (2013).
- kSNP4.1 User Guide (SourceForge).
- SKA2: Tonkin-Hill et al. 2024, *Genome Research* 34(10):1661.
- SKA (original): Harris 2018, bioRxiv 453142.
- ska.rust: github.com/bacpop/ska.rust ; docs.rs/ska.
- skalo: bioRxiv 2024.10.02.616334.
- OptiK: bioRxiv 2025.05.21.655412.
- Building trees with SKA: bacpop.org/guides/building_trees_with_ska.
