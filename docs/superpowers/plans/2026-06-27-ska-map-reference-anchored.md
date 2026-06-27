# Reference-anchored mapping (`ska map`) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Expose `ska map`'s reference-coordinate pseudo-alignment as a first-class skTree output and build optional reference-anchored trees from it, reusing the existing SNP-extraction and tree engines.

**Architecture:** When `--reference` is set, the pipeline calls the existing `SkaRunner.map(..., out_format="aln")` to write `ref_aligned.fasta`. A new opt-in `--map-tree` flag runs that pseudo-alignment through `parse_alignment`/`classify` and the existing `neighbor_joining` / `parsimony_tree` / `MlTreeBuilder`, producing a parallel `tree_ref_*` set that honors `--ml`/`--parsimony`.

**Tech Stack:** Python 3.12, numpy, dendropy, pytest, argparse. SKA2 `ska` CLI (already wrapped).

## Global Constraints

- English only; self-documenting code; sparing comments.
- Conventional commits, atomic, split by concern; no mention of CLAUDE; no sign-off.
- Graceful-skip pattern for external-tool failures (log warning, continue) — matches ML/clustering.
- `--map-tree` requires `--reference`; argparse error otherwise (fail fast).
- Reuse existing engines; add no new tree algorithm.
- NJ/parsimony run on the *variable* subset of the pseudo-alignment; ML runs on the *core* subset file (mirrors the reference-free path exactly).

---

### Task 1: Lock the `ska map -f aln` invocation with a test

**Files:**
- Test: `tests/test_ska_runner.py` (add one test next to `test_map_argv`)

**Interfaces:**
- Consumes: `SkaRunner.map(reference, inputs, output=, out_format="aln", threads=)` (exists).

- [ ] **Step 1: Add the test**

```python
def test_map_aln_argv(fake_run, tmp_path):
    runner = SkaRunner()
    ref = tmp_path / "ref.fasta"
    skf = tmp_path / "all.skf"
    out = runner.map(ref, skf, output=tmp_path / "ref.aln", out_format="aln", threads=4)
    argv = fake_run[0][0]
    assert argv[1] == "map"
    assert "-f" in argv and "aln" in argv
    assert "--threads" in argv and "4" in argv
    assert argv.index(str(ref)) < argv.index(str(skf))
    assert out == tmp_path / "ref.aln"
```

- [ ] **Step 2: Run it**

Run: `.venv/bin/python -m pytest tests/test_ska_runner.py::test_map_aln_argv -q`
Expected: PASS (the wrapper already supports `aln`).

- [ ] **Step 3: Commit**

```bash
git add tests/test_ska_runner.py
git commit -m "test(ska): cover ska map -f aln invocation"
```

---

### Task 2: Emit `ref_aligned.fasta` when `--reference` is set

**Files:**
- Modify: `sktree/pipeline.py`
- Test: `tests/test_pipeline.py`

**Interfaces:**
- Produces: `_reference_align(runner, skf, reference, outdir, threads) -> Path | None`
  (returns the written `ref_aligned.fasta`, or `None` on `SkaError`).

- [ ] **Step 1: Add the helper** (place after `_annotate_snps`)

```python
def _reference_align(
    runner: SkaRunner,
    skf: Path,
    reference: Path,
    outdir: Path,
    threads: int | None,
) -> Path | None:
    """Write the reference-coordinate pseudo-alignment via ``ska map -f aln``."""
    out = outdir / "ref_aligned.fasta"
    logger.info("Mapping split k-mers onto reference %s (pseudo-alignment)", reference.name)
    try:
        runner.map(reference, skf, output=out, out_format="aln", threads=threads)
    except SkaError as exc:
        logger.warning("Reference mapping failed: %s", exc)
        return None
    return out
```

Add `from .engine.ska import SkaError, SkaRunner` (extend the existing import).

- [ ] **Step 2: Call it in `run_pipeline`** — in the `if reference is not None:` block, before/after annotation:

```python
    ref_alignment_path: Path | None = None
    annotation_path: Path | None = None
    if reference is not None:
        ref_alignment_path = _reference_align(runner, skf, reference, outdir, threads)
        annotation_path = _annotate_snps(runner, skf, reference, outdir, threads)
```

- [ ] **Step 3: Add `ref_alignment` to `RunResult`** (field `ref_alignment: Path | None = None`) and set it in the returned `RunResult`.

- [ ] **Step 4: Add a test**

```python
@needs_ska
def test_reference_pseudo_alignment_written(synthetic_genomes, tmp_path):
    out = tmp_path / "out"
    reference = synthetic_genomes[0]
    result = run_pipeline(synthetic_genomes, out, k=31, reference=reference)
    assert result.ref_alignment is not None
    assert result.ref_alignment.exists()
    assert (out / "ref_aligned.fasta").exists()
    assert result.ref_alignment.read_text().count(">") >= 1
```

- [ ] **Step 5: Run + commit**

Run: `.venv/bin/python -m pytest tests/test_pipeline.py -q`
```bash
git add sktree/pipeline.py tests/test_pipeline.py
git commit -m "feat(pipeline): write reference-anchored pseudo-alignment with --reference"
```

---

### Task 3: Build reference-anchored trees (`_build_map_trees` + `--map-tree`)

**Files:**
- Modify: `sktree/pipeline.py`, `sktree/cli.py`
- Test: `tests/test_pipeline.py`, `tests/test_cli.py`

**Interfaces:**
- Consumes: `parse_alignment`, `SnpMatrix.classify/subset/to_fasta`, `neighbor_joining`, `to_newick`, `parsimony_tree`, `_build_ml_tree`.
- Produces: `MapOutputs` dataclass; `_build_map_trees(ref_alignment, outdir, *, parsimony, ml, threads) -> MapOutputs`; `run_pipeline(..., map_tree=False)`.

- [ ] **Step 1: Parameterize `_build_ml_tree` prefix** — change signature to
  `def _build_ml_tree(alignment, outdir, threads, *, prefix_name="tree_ml")` and use
  `prefix = outdir / prefix_name`. Existing call stays default.

- [ ] **Step 2: Add `MapOutputs` + `_build_map_trees`**

```python
@dataclass
class MapOutputs:
    nj_tree: Path
    parsimony_tree: Path | None
    ml_tree: Path | None
    ml_report: Path | None


def _build_map_trees(
    ref_alignment: Path,
    outdir: Path,
    *,
    parsimony: bool,
    ml: bool,
    threads: int | None,
) -> MapOutputs:
    """Build reference-anchored trees from the pseudo-alignment.

    NJ and parsimony run on the variable subset; ML runs on the core-SNP
    subset file, mirroring the reference-free path. The full-length scan in
    classify() is the performance watch-point on large references.
    """
    full = parse_alignment(ref_alignment)
    cls = full.classify()
    variable = full.subset(cls.variable)

    nj = outdir / "tree_ref_nj.nwk"
    nj.write_text(to_newick(neighbor_joining(variable)))

    parsimony_path: Path | None = None
    if parsimony:
        tree, _ = parsimony_tree(variable)
        parsimony_path = outdir / "tree_ref_parsimony.nwk"
        parsimony_path.write_text(to_newick(tree))

    ml_tree = ml_report = None
    if ml:
        core = full.subset(cls.core)
        core_fasta = outdir / "ref_core_snps.fasta"
        core_fasta.write_text(core.to_fasta())
        ml_out = _build_ml_tree(core_fasta, outdir, threads, prefix_name="tree_ref_ml")
        if ml_out is not None:
            ml_tree = ml_out.tree
            ml_report = ml_out.report
    return MapOutputs(nj_tree=nj, parsimony_tree=parsimony_path, ml_tree=ml_tree, ml_report=ml_report)
```

- [ ] **Step 3: Wire into `run_pipeline`** — add `map_tree: bool = False` param; after `ref_alignment_path` is set:

```python
    map_outputs: MapOutputs | None = None
    if map_tree and ref_alignment_path is not None:
        logger.info("Building reference-anchored tree(s) from the pseudo-alignment")
        map_outputs = _build_map_trees(
            ref_alignment_path, outdir, parsimony=parsimony, ml=ml, threads=threads
        )
```

Add `RunResult` fields: `map_nj_tree`, `map_parsimony_tree`, `map_ml_tree`, `map_ml_report` (all `Path | None = None`); populate from `map_outputs`.

- [ ] **Step 4: CLI flag** — in `sktree/cli.py` add after `--reference`:

```python
    run.add_argument(
        "--map-tree", action="store_true",
        help="build reference-anchored tree(s) from the ska map pseudo-alignment "
             "(requires --reference; honors --ml/--parsimony)",
    )
```

After `args = parser.parse_args(argv)`:

```python
    if getattr(args, "map_tree", False) and args.reference is None:
        parser.error("--map-tree requires --reference")
```

Pass `map_tree=args.map_tree` into `run_pipeline`. Print the map NJ tree path in the summary when present.

- [ ] **Step 5: Tests**

```python
@needs_ska
def test_map_tree_requires_reference_via_cli(synthetic_genomes, tmp_path):
    from sktree.cli import main
    with pytest.raises(SystemExit):
        main(["run", *[str(p) for p in synthetic_genomes], "-o", str(tmp_path / "o"), "--map-tree"])


@needs_ska
def test_map_tree_builds_nj(synthetic_genomes, tmp_path):
    out = tmp_path / "out"
    result = run_pipeline(synthetic_genomes, out, k=31, reference=synthetic_genomes[0], map_tree=True)
    assert result.map_nj_tree is not None and result.map_nj_tree.exists()
    assert (out / "tree_ref_nj.nwk").exists()
```

- [ ] **Step 6: Run + commit**

Run: `.venv/bin/python -m pytest tests/test_pipeline.py tests/test_cli.py -q`
```bash
git add sktree/pipeline.py sktree/cli.py tests/test_pipeline.py tests/test_cli.py
git commit -m "feat(pipeline): build reference-anchored trees behind --map-tree"
```

---

### Task 4: Surface the new outputs in the HTML report + output listing

**Files:**
- Modify: `sktree/pipeline.py` (`_output_listing` block and `trees` dict)

- [ ] **Step 1: Add output-listing rows** (inside the `if html:` block list) for
  `ref_alignment_path` ("Reference-anchored pseudo-alignment (ska map, one row per sample at reference coordinates)"),
  `map_outputs.nj_tree` ("Reference-anchored neighbor-joining tree (Newick)"),
  `map_outputs.parsimony_tree`, `map_outputs.ml_tree`, `map_outputs.ml_report` — each guarded by `map_outputs is not None`.

- [ ] **Step 2: Add reference-anchored trees to the `trees` dict**

```python
        if map_outputs is not None:
            trees["ref_nj"] = map_outputs.nj_tree.read_text().strip()
            if map_outputs.parsimony_tree is not None:
                trees["ref_parsimony"] = map_outputs.parsimony_tree.read_text().strip()
            if map_outputs.ml_tree is not None:
                trees["ref_ml"] = map_outputs.ml_tree.read_text().strip()
```

- [ ] **Step 3: Run the full suite**

Run: `.venv/bin/python -m pytest -q`
Expected: all pass / prior skips unchanged.

- [ ] **Step 4: Commit**

```bash
git add sktree/pipeline.py
git commit -m "feat(report): list reference-anchored outputs and trees in HTML report"
```

---

### Task 5: End-to-end demo run + performance note

**Files:** none (validation), captures timing for docs.

- [ ] **Step 1: Run the demo with the new flags**

Run (reference = first demo input):
`.venv/bin/python -m sktree.cli run <demo inputs> -o /tmp/.../map_demo --reference <first input> --map-tree --html -v`
Expected: `ref_aligned.fasta`, `tree_ref_nj.nwk`, `report.html` written; note wall-clock.

- [ ] **Step 2: Confirm gap char + reference-row behavior** from `ref_aligned.fasta`; if the reference appears as its own row and is unwanted as a tip, note it (no code change unless it breaks the tree taxa).

---

### Task 6: Documentation — README applications section + FOR-DEVELOPERS note

**Files:**
- Modify: `README.md`, `FOR-DEVELOPERS.md`

- [ ] **Step 1: README** — document `--map-tree`; add output-table rows for `ref_aligned.fasta` and `tree_ref_*`; add the **"Applications: `ska map` vs read mapping (bwa/minimap2)"** section (tree building, AMR, genes of interest, recommendation table) per the spec.

- [ ] **Step 2: FOR-DEVELOPERS** — brief note: map-tree reuses the SNP-extraction layer; full-length-alignment performance caveat (`classify` scan).

- [ ] **Step 3: Commit**

```bash
git add README.md FOR-DEVELOPERS.md
git commit -m "docs: document --map-tree and ska map vs read-mapping applications"
```

## Self-Review

- Spec coverage: CLI (T3/T4), pseudo-alignment output (T2), map trees reusing engines (T3), HTML (T4), perf note (T5/T6), applications doc (T6) — all covered.
- No placeholders; all code shown.
- Type consistency: `MapOutputs` fields referenced identically in T3/T4; `_build_ml_tree` `prefix_name` kw used in T3 matches T3 Step 1 signature change; `RunResult` fields named consistently.
