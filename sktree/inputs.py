"""Resolve user inputs into the samples ``ska build`` consumes.

A *sample* is one tree tip. It is backed by either a single sequence file (an
assembly, or a single-file input) or a pair of FASTQ files (paired-end reads).
This module is the single owner of two concerns the rest of skTree should not
have to think about:

* **Grouping** -- turning a flat list of CLI paths into samples, auto-pairing
  paired-end FASTQ files (``_R1``/``_R2``, ``_1``/``_2`` ...), and parsing a
  manifest TSV in SKA's own ``-f`` format (``name<TAB>file[<TAB>file2]``).
* **Naming** -- deriving tidy, unique tip labels so trees and alignments are
  readable instead of named after verbatim paths.

SKA requires FASTQ to be paired-end (two deinterleaved files); a lone FASTQ is
silently treated as FASTA, so we fail loudly when reads cannot be paired.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

PathLike = str | Path

_COMPRESSION_SUFFIXES = {".gz", ".bz2", ".xz", ".zst", ".zip"}
_SEQUENCE_SUFFIXES = {".fasta", ".fa", ".fna", ".ffn", ".faa", ".fastq", ".fq"}
_FASTQ_SUFFIXES = {".fastq", ".fq"}

# Mate markers, most specific first. Each captures the mate number (1 or 2); the
# part of the stem before the marker is the shared pair key / sample name.
_MATE_PATTERNS = (
    re.compile(r"^(?P<base>.+)_R(?P<mate>[12])(?:_\d+)?$"),  # _R1 / _R2 / _R1_001
    re.compile(r"^(?P<base>.+)\.R(?P<mate>[12])$"),          # .R1 / .R2
    re.compile(r"^(?P<base>.+)_(?P<mate>[12])$"),            # _1 / _2
)


@dataclass(frozen=True)
class Sample:
    """One tree tip: a name plus its backing sequence file(s)."""

    name: str
    files: tuple[Path, ...]
    is_reads: bool


def _strip_known_suffixes(name: str) -> str:
    """Drop trailing compression and sequence extensions (``a.fa.gz`` -> ``a``)."""
    base = name
    for suffix in reversed(Path(name).suffixes):
        if suffix.lower() in _COMPRESSION_SUFFIXES or suffix.lower() in _SEQUENCE_SUFFIXES:
            base = base[: -len(suffix)]
        else:
            break
    return base or name


def clean_sample_name(path: PathLike) -> str:
    """Derive a tidy sample label from an input file path.

    Strips the directory and any trailing compression and FASTA/FASTQ
    extensions, so ``/data/SAMD00052601.fa.gz`` becomes ``SAMD00052601``.
    """
    return _strip_known_suffixes(Path(path).name)


def _is_fastq(path: PathLike) -> bool:
    """True if ``path`` looks like a FASTQ file by extension (ignoring .gz)."""
    suffixes = [s.lower() for s in Path(path).suffixes]
    relevant = [s for s in suffixes if s not in _COMPRESSION_SUFFIXES]
    return bool(relevant) and relevant[-1] in _FASTQ_SUFFIXES


def _mate_of(path: Path) -> tuple[str, int] | None:
    """Return ``(pair_key, mate_number)`` for a FASTQ path, or ``None``.

    The pair key is the sample stem with the mate marker removed, so the two
    files of a pair share it.
    """
    stem = clean_sample_name(path)
    for pattern in _MATE_PATTERNS:
        match = pattern.match(stem)
        if match:
            return match.group("base"), int(match.group("mate"))
    return None


def _pair_positional_reads(fastqs: list[Path]) -> list[Sample]:
    """Group positional FASTQ paths into paired-end samples, fail-loud on gaps."""
    by_key: dict[str, dict[int, Path]] = {}
    order: list[str] = []
    for path in fastqs:
        mate = _mate_of(path)
        if mate is None:
            raise ValueError(
                f"Cannot find a paired mate for FASTQ '{path}': expected a "
                "_R1/_R2, _1/_2 or .R1/.R2 naming. Provide both mates, or use "
                "--manifest to declare the pairing explicitly."
            )
        key, num = mate
        if key not in by_key:
            by_key[key] = {}
            order.append(key)
        by_key[key][num] = path

    samples: list[Sample] = []
    for key in order:
        mates = by_key[key]
        if set(mates) != {1, 2}:
            present = next(iter(mates.values()))
            raise ValueError(
                f"FASTQ '{present}' has no paired mate (need both _1 and _2 for "
                f"sample '{key}'). Provide both mates, or use --manifest."
            )
        samples.append(Sample(key, (mates[1], mates[2]), is_reads=True))
    return samples


def _parse_manifest(manifest: Path) -> list[Sample]:
    """Parse a ``name<TAB>file[<TAB>file2]`` TSV (SKA's own ``-f`` format)."""
    samples: list[Sample] = []
    for lineno, raw in enumerate(manifest.read_text().splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        fields = [f for f in raw.split("\t") if f != ""]
        if len(fields) == 2:
            name, f1 = fields
            samples.append(Sample(name, (Path(f1),), is_reads=_is_fastq(f1)))
        elif len(fields) == 3:
            name, f1, f2 = fields
            samples.append(Sample(name, (Path(f1), Path(f2)), is_reads=True))
        else:
            raise ValueError(
                f"{manifest}:{lineno}: expected 2 columns (name, file) or 3 "
                f"columns (name, fwd, rev), got {len(fields)}: {raw!r}"
            )
    return samples


def _dedupe_names(samples: list[Sample]) -> list[Sample]:
    """Append numeric suffixes so every sample name is unique."""
    seen: dict[str, int] = {}
    result: list[Sample] = []
    for s in samples:
        if s.name in seen:
            seen[s.name] += 1
            new_name = f"{s.name}_{seen[s.name]}"
            result.append(Sample(new_name, s.files, s.is_reads))
        else:
            seen[s.name] = 1
            result.append(s)
    return result


def resolve_inputs(
    positional: Sequence[PathLike],
    manifest: PathLike | None,
) -> list[Sample]:
    """Turn CLI inputs into a de-duplicated list of :class:`Sample`.

    Manifest samples (if any) come first, then positional samples. Positional
    FASTA files become one sample each; positional FASTQ files are auto-paired.
    """
    samples: list[Sample] = []
    if manifest is not None:
        samples.extend(_parse_manifest(Path(manifest)))

    assemblies: list[Path] = []
    fastqs: list[Path] = []
    for raw in positional:
        path = Path(raw)
        (fastqs if _is_fastq(path) else assemblies).append(path)

    samples.extend(_pair_positional_reads(fastqs))
    samples.extend(
        Sample(clean_sample_name(p), (p,), is_reads=False) for p in assemblies
    )
    return _dedupe_names(samples)
