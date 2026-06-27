# Reference-anchored mapping (`ska map`) feature — Design

**Date:** 2026-06-27
**Status:** Approved for planning

## Summary

skTree's SNP discovery is currently reference-free: `ska build` then `ska align`
produce an unordered, SNP-only alignment that feeds the NJ/parsimony/ML trees.
SKA2 also offers `ska map`, which walks each sample's split k-mers against a
reference genome and emits a **reference-coordinate** result — either a VCF or a
reference-length pseudo-alignment.

The `SkaRunner.map()` wrapper already exists (`sktree/engine/ska.py`) and is used
internally for one purpose only: producing `reference_snps.vcf` to drive gene
annotation when the user passes `--reference`. The reference-anchored
**pseudo-alignment** (`-f aln`) is never exposed.

This feature exposes that capability: when `--reference` is given, skTree also
writes the reference-coordinate pseudo-genome alignment, and a new opt-in
`--map-tree` flag builds reference-anchored tree(s) from it using the existing
tree engines. A README section explains when to prefer `ska map` over read
mapping (bwa/minimap2) for tree building and AMR / gene-of-interest work.

## Goals

- Expose `ska map -f aln` as a first-class output (`ref_aligned.fasta`).
- Build reference-anchored trees from the pseudo-alignment, reusing the existing
  SNP-extraction and tree-building layers (NJ always; parsimony/ML when the
  matching flags are set).
- Document the applications trade-off (`ska map` vs bwa/minimap2) for tree
  building, AMR, and genes of interest.

## Non-goals

- No standalone `sktree map` subcommand. The capability rides on `run`.
- No gene-presence/absence or accessory-genome AMR detection. `ska map` cannot
  do it; that is precisely the limitation the docs call out.
- No new tree algorithms. The map path reuses `neighbor_joining`,
  `parsimony_tree`, and `MlTreeBuilder` unchanged.

## CLI surface

- `--reference PATH` (existing): in addition to its current SNP annotation, now
  also writes `ref_aligned.fasta` (one `ska map -f aln` call).
- `--map-tree` (new, store_true): build reference-anchored tree(s) from the
  pseudo-alignment. Requires `--reference`; if given without it, the CLI exits
  with an argument error before any computation.
- `--map-tree` honors the existing `--ml` and `--parsimony` flags to decide
  which engines run on the reference-anchored alignment.

## Data flow

```
--reference REF ─┬─► ska map -f vcf ─► reference_snps.vcf ─► snp_annotation.tsv   (existing)
                 └─► ska map -f aln ─► ref_aligned.fasta                          (new, always when --reference)
                                            │
                              --map-tree ───┴─► parse_alignment → classify (core/variable)
                                                     ├─► tree_ref_nj.nwk          (always)
                                                     ├─► tree_ref_parsimony.nwk   (if --parsimony)
                                                     └─► tree_ref_ml.*            (if --ml)
```

## Components and files

- **`sktree/cli.py`** — add `--map-tree` (store_true); validate it requires
  `--reference` (argparse error otherwise); pass `map_tree` into `run_pipeline`.
- **`sktree/pipeline.py`**
  - When `--reference` is set, call `runner.map(reference, skf,
    output=outdir / "ref_aligned.fasta", out_format="aln", threads=threads)`.
    Wrap in the same graceful-skip pattern as ML/clustering: on `SkaError`, log a
    warning and continue.
  - New `MapOutputs` dataclass: `alignment: Path`, `nj_tree: Path | None`,
    `parsimony_tree: Path | None`, `ml_tree: Path | None`,
    `ml_report: Path | None`.
  - New `_build_map_trees(ref_alignment, outdir, *, parsimony, ml, threads)`:
    parse the pseudo-alignment with `parse_alignment`, `classify`, extract the
    core subset, write `tree_ref_nj.nwk` (always), and `tree_ref_parsimony.nwk`
    / `tree_ref_ml.*` when the flags are set. Reuses `_build_ml_tree`'s
    collection logic for the ML family (prefix `tree_ref_ml`).
  - Extend `RunResult` with `ref_alignment`, `map_nj_tree`,
    `map_parsimony_tree`, `map_ml_tree`, `map_ml_report` (all `Path | None`).
  - Add the new files to `_output_listing`.
  - Add reference-anchored trees to the HTML `trees` dict under keys
    `ref_nj`, `ref_parsimony`, `ref_ml`.
- **`sktree/engine/ska.py`** — no change; `map()` already supports `aln`.
- **HTML report** — the reference-anchored trees surface in the existing tree
  switcher automatically once present in `trees`; only human-readable label
  strings are added (no JS/CSS structural change).

## Missing-data and reference-row semantics

`ska map -f aln` emits one row per sample at reference coordinates; positions
where a sample's split k-mer does not match (divergence, indels, accessory
regions) are output as missing. `classify` already treats non-ACGT bases as
missing, so the **core** subset is exactly the set of positions resolved in
every sample — the correct filter for a tree.

Two details are confirmed against real `ska map` output on the demo set during
implementation, and the parser normalized if needed:
1. The exact missing/gap character (`-` vs `N`).
2. Whether the reference itself appears as its own row (and therefore as a tree
   tip). If it does and that is undesirable, drop it before tree building.

## Performance consideration

The reference-free alignment fed to the engines is SNP-only and tiny. The
pseudo-alignment is **full reference length** — for a multi-megabase genome
across dozens of samples, pure-Python `classify` scans a large character matrix.
This is the one performance watch-point. Implementation measures it on the demo
set. If unacceptable, the fallback (not taken now, to avoid extra code) is to
reconstruct variable sites from the VCF rather than scanning the full
pseudo-alignment. This caveat is recorded in FOR-DEVELOPERS.md.

## Error handling

- `ska map` failure (either format): log a warning, skip the dependent output,
  continue the run — matching the ML/clustering graceful-skip pattern.
- `--map-tree` without `--reference`: hard CLI argument error, fail fast before
  any compute.

## Testing

- `tests/test_ska_runner.py`: assert `map(..., out_format="aln")` builds the
  expected argv (`map -f aln ... REF inputs`) and returns the output path.
- `tests/test_pipeline.py`:
  - `--reference` writes `ref_aligned.fasta`.
  - `--map-tree` produces `tree_ref_nj.nwk`.
  - `--map-tree` with `--ml` / `--parsimony` produces the matching
    `tree_ref_ml.*` / `tree_ref_parsimony.nwk`.
  - `--map-tree` without `--reference` raises an argument error.
- End-to-end demo run to validate outputs and capture the performance note.

## Documentation

- **README.md**
  - `--reference` note updated: also emits `ref_aligned.fasta`.
  - New `--map-tree` flag documented.
  - Output table rows for `ref_aligned.fasta` and the `tree_ref_*` files.
  - New section **"Applications: `ska map` vs read mapping (bwa/minimap2)"**:
    - *Tree building* — reference-anchored gives genomic coordinates (mask
      repeats/recombination, cross-study comparability) at the cost of reference
      bias on recovered sites; reference-free `align` stays bias-free. Guidance
      on which to pick.
    - *AMR* — `ska map` resolves point-mutation resistance in conserved genes
      (gyrA, rpoB, etc.) fast and alignment-free, but is blind to acquired /
      accessory AMR genes (plasmid-borne, absent from reference), indels, and
      copy number, because split k-mers only map where flanks match. bwa/minimap2
      plus coverage-based callers (AMRFinderPlus, ARIBA) handle presence/absence,
      indels, and novel genes.
    - *Genes of interest* — SNP context within reference-present genes via
      `ska map`; presence/absence and structural variation via read mapping.
    - A short recommendation table.
- **FOR-DEVELOPERS.md** — brief note on why the map-tree path reuses the
  SNP-extraction layer, plus the full-length-alignment performance caveat.

## Implementation order

1. `ska map -f aln` wrapper test (lock the invocation).
2. Pipeline: emit `ref_aligned.fasta` when `--reference` set; tests.
3. Pipeline: `_build_map_trees` + `MapOutputs` + `RunResult` fields; tests.
4. CLI: `--map-tree` flag + `--reference` requirement; tests.
5. HTML report: surface `ref_*` trees + output listing entries.
6. End-to-end demo run; capture performance note.
7. README applications section + FOR-DEVELOPERS note.
