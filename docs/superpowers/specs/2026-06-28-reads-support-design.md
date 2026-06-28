# Design: Build trees from sequencing reads (FASTQ) alongside assemblies

**Date:** 2026-06-28
**Status:** Approved (brainstorming)

## Goal

Let `sktree run` accept paired-end sequencing reads (FASTQ) as samples, freely
mixed with assemblies (FASTA), so a single tree can be built from both. SKA2
already supports this natively; the work is to expose it at skTree's CLI and
input-resolution layer and to error-filter reads sensibly.

## Background — what SKA gives us (verified on ska 0.5.1)

- `ska build -f <list>` reads a TSV where each line is either
  `name<TAB>file` (assembly or single file) or `name<TAB>fwd<TAB>rev`
  (paired-end reads). Assemblies and read pairs may be mixed in one list.
- skTree's `build()` *already* writes this `name<TAB>path` file list — today it
  only ever emits the one-file-per-sample form.
- File type is auto-detected by extension (`.fastq/.fq` ± compression vs FASTA).
  **FASTQ must be paired-end (two deinterleaved files).** A lone FASTQ is
  silently treated as FASTA, producing wrong results.
- Read builds need error filtering: `--min-count` (drop low-count k-mers from
  sequencing error) plus `--min-qual`/`--qual-filter` (already plumbed in
  `build()`). Empirically, `--min-count` is **ignored for assembly inputs**, so a
  global value is safe in mixed runs — but we still only pass it when reads are
  present, to keep pure-assembly runs unchanged.
- A mixed assembly+paired-read build with `--min-count 3` correctly recovered a
  known SNP (assembly base vs read-derived base) in testing.

## The gap

`sktree run` takes a flat positional list where one file = one sample. There is
no way to (a) group two FASTQ files into one paired sample or (b) mix reads with
assemblies. Two secondary touch-points:

- `kselect.select_k` parses inputs as plaintext FASTA; it breaks on FASTQ/gzip.
- `--min-count` is not surfaced on the CLI.

## Design

### 1. Input model — new `sktree/inputs.py`

Decouple *how the user specified files* from *what SKA builds*:

```python
@dataclass(frozen=True)
class Sample:
    name: str
    files: tuple[Path, ...]   # 1 file = assembly/single; 2 = paired reads
    is_reads: bool            # any file is FASTQ (by extension)

def resolve_inputs(
    positional: Sequence[Path],
    manifest: Path | None,
) -> list[Sample]: ...
```

Behaviour:

- **Manifest** (`--manifest samples.tsv`): each non-blank, non-`#` line is
  `name<TAB>file[<TAB>file2]`. Two columns → assembly/single; three → paired
  reads. This is SKA's own `-f` format, so it maps 1:1.
- **Auto-pair positional**: FASTQ files are paired on mate tokens `_R1/_R2`,
  `_1/_2` (including the Illumina `_R1_001` lane suffix); the derived sample name
  is the filename with the mate token and sequence/compression extensions
  stripped. FASTA files are always one-per-sample.
- **Mixing**: manifest samples and positional samples combine into one list.
- **Validation (fail loud):**
  - A lone positional FASTQ with no detectable mate → error: provide both mates
    or use `--manifest` (never silently treat reads as FASTA).
  - Duplicate sample names are disambiguated with a numeric suffix, as today.
- `is_reads` is true when any file in the sample is FASTQ by extension.

`clean_sample_name` and the name-dedup logic move from `engine/ska.py` into this
resolver, which becomes the single owner of sample naming.

### 2. Engine change — `engine/ska.py`

- `build()` takes `samples: Sequence[Sample]` instead of a flat `seq_files`
  list. It writes one file-list line per sample: `name<TAB>f1[<TAB>f2]`.
- `min_count` is passed to `ska build` only when at least one sample
  `is_reads`; pure-assembly runs invoke `build` exactly as before.
- `min_qual` / `qual_filter` are forwarded when provided. (`qual_filter` is a new
  pass-through arg mirroring SKA's `--qual-filter`.)

### 3. Pipeline & CLI

- `run_pipeline` calls `resolve_inputs(...)` first, then passes `Sample`s to
  `build()`. The two-input minimum is enforced on resolved **samples**, not raw
  files (a paired sample is two files but one sample).
- New CLI flags on `sktree run`:
  - `--manifest PATH` — TSV sample sheet.
  - `--min-count N` (default 3) — applied to read samples only.
  - `--min-qual N` — forwarded to `ska build`.
  - `--qual-filter {no-filter,middle,strict}` — forwarded to `ska build`.
- **auto-k with reads**: `select_k` receives only the *assembly* samples'
  files. If a run is reads-only and `--auto-k` is set, log a warning and fall
  back to the `-k` default (FASTQ uniqueness is both unparsed and noisy from
  read error). Mixed runs select k from the assemblies present.
- Reference / annotate / `ska map` paths are unaffected: they operate on the
  produced `.skf`.

### 4. Tests

- **Unit (no ska):**
  - resolver: auto-pair token variants (`_R1/_R2`, `_1/_2`, `_R1_001`), manifest
    parse (2- and 3-column, comments), lone-FASTQ error, mixed positional +
    manifest, duplicate-name dedup.
  - `build()` writes correct 2- and 3-column file-list lines and only adds
    `--min-count` when a read sample is present.
- **Integration (`requires_ska`):** a synthetic paired-FASTQ fixture (tile
  ~150 bp reads off a mutated copy of a genome) plus an assembly → mixed run
  recovers the known SNP between them.

## Out of scope (YAGNI)

- Automatic `ska cov` coverage-cutoff estimation (use fixed `--min-count` + flag).
- Single-end or interleaved FASTQ.
- Per-sample `--min-count` overrides.

## Files touched

- `sktree/inputs.py` (new) — `Sample`, `resolve_inputs`, sample naming.
- `sktree/engine/ska.py` — `build()` takes `Sample`s; `clean_sample_name` /
  dedup move out to `inputs.py`.
- `sktree/pipeline.py` — resolve inputs, plumb `min_count`/`min_qual`/
  `qual_filter`, assembly-only auto-k.
- `sktree/cli.py` — `--manifest`, `--min-count`, `--min-qual`, `--qual-filter`.
- `sktree/kselect.py` — accept the assembly subset (caller-filtered).
- `tests/` — resolver unit tests, mixed-input integration test, paired-FASTQ
  fixture.
- Docs: README + FOR-DEVELOPERS.md note on reads input.
