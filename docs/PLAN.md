# skTree — Implementation Plan

Goal: a Python CLI that wraps **SKA2** into an end-to-end reference-free,
alignment-free SNP discovery + phylogenetics workflow (SKA2 subprocess-wrapped).
See `RESEARCH.md` for the why.

## Decisions
- **Engine integration:** wrap the `ska` CLI as a subprocess.
- **Trees:** NJ + parsimony (pure Python, DendroPy), ML (wrap IQ-TREE/RAxML).
- **MVP order:** core pipeline first, then k-selection, then parsimony + ML.
- **Testing:** TDD. Unit tests use *recorded/synthetic* SKA outputs (no binary
  required to run the suite); a small set of integration tests run the real `ska`
  binary when present (skipped otherwise).
- **Packaging:** `uv` + `pyproject.toml`, `ruff` for lint, `pytest` for tests.

## Milestones

### M1 — Project scaffold
- `pyproject.toml` (uv), package `sktree/`, `tests/`, ruff config, pytest config.
- `sktree/cli.py` skeleton with `argparse` subcommands: `run`, `build`, `align`.

### M2 — SKA engine wrapper (`engine/ska.py`)  ← TDD starts here
- `SkaRunner` wrapping `ska build`, `ska merge`, `ska align`, `ska distance`.
- Resolve binary path; raise a clear error if missing; capture stdout/stderr.
- Tests mock `subprocess.run` → assert correct argv + output parsing.

### M3 — SNP parsing & partition (`snps.py`)
- Parse `ska align` FASTA alignment → `SnpMatrix` (samples × variable sites).
- Partition loci into **core** (no missing), **majority** (present ≥ threshold),
  **all**. Configurable missing-fraction threshold.
- Variable-site filtering (drop invariant/all-missing columns).
- Tests on hand-built alignments with known partitions.

### M4 — NJ tree + report (`phylo.py`, `report.py`)
- Distance matrix from alignment (or parse `ska distance`).
- NJ tree via DendroPy → Newick.
- Summary report: sample count, SNP counts per category, chosen k, run params.
- End-to-end `sktree run <inputs>` produces alignment, matrix, NJ tree, report.

### M5 — k selection (`kselect.py`)
- uniqueness scan: scan odd k, pick smallest k where k-mers are ~unique (FCK target).
- `--k` to override; `--auto-k` to select.

### M6 — Parsimony + ML trees  ✅ done
- Pure-Python maximum parsimony (Fitch scoring + stepwise addition + NNI).
- ML: wrap IQ-TREE / RAxML-NG (subprocess), bootstrap-capable; skips if absent.
- Both exposed via `--parsimony` / `--ml`; outputs `tree_parsimony.nwk`,
  `tree_ml.treefile`.

### M7 — Polish  ✅ done
- FOR-DEVELOPERS.md (architecture deep-dive) + README outputs/options tables.
- Validation: smoke-tested on deterministic non-repetitive genomes; fixed the
  divergence estimate to use SNP rate, not k-mer mismatch proportion.

### M8 — Reference-based SNP annotation (`annotate.py`)  ✅ done
- `SkaRunner.map()` wraps `ska map -f vcf` for reference-coordinate SNPs.
- pyrodigal predicts genes on the reference (metagenomic mode, no training pass).
- Each SNP labelled intergenic vs coding; coding SNPs get a synonymous /
  non-synonymous call by translating the affected codon (minus-strand aware).
- Exposed via `--reference REF`; outputs `reference_snps.vcf` and
  `snp_annotation.tsv`. Optional `annotate` extra (`pip install sktree[annotate]`).

## Output contract (`sktree run`)
```
outdir/
  combined.skf            # merged split-k-mer file
  alignment.fasta         # reference-free SNP alignment
  snp_matrix.tsv          # samples × loci, with core/majority labels
  core_snps.fasta         # core-only alignment
  tree_nj.nwk             # neighbor-joining tree
  tree_parsimony.nwk      # (M6, --parsimony)
  tree_ml.treefile        # (M6, --ml, if ML binary present)
  reference_snps.vcf      # (M8, --reference) SNPs in reference coordinates
  snp_annotation.tsv      # (M8, --reference) per-SNP gene + codon effect
  summary.txt             # run summary
```
