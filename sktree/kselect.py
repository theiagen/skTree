"""Optimal odd-k selection, a uniqueness-based heuristic.

The heuristic picks the smallest odd k at which k-mers in a representative
genome are *mostly unique*: too small and unrelated loci collide (false
homology); too large and real SNPs are lost or swamped by sequencing error.

We reproduce the idea directly: for the median-length input genome, scan odd k
and pick the smallest k whose canonical-k-mer uniqueness reaches a target
(default 0.99). Counting k-mers for a single genome in Python is cheap, so no
external engine is needed for this step.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _reverse_complement(kmer: str) -> str:
    return kmer.translate(_COMPLEMENT)[::-1]


def canonical_kmer(kmer: str) -> str:
    """Strand-canonical form: the lexicographically smaller of k-mer and its RC."""
    rc = _reverse_complement(kmer)
    return kmer if kmer <= rc else rc


def _concatenated_sequence(fasta_path: Path) -> str:
    parts: list[str] = []
    for line in fasta_path.read_text().splitlines():
        if line and not line.startswith(">"):
            parts.append(line.strip().upper())
    return "".join(parts)


def kmer_uniqueness(seq: str, k: int) -> float:
    """Fraction of k-mer occurrences whose canonical k-mer is unique in ``seq``.

    K-mers containing any non-A/C/G/T character are skipped. Returns 0.0 when no
    valid k-mer exists.
    """
    counts: dict[str, int] = {}
    order: list[str] = []
    valid = set("ACGT")
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        if not set(kmer) <= valid:
            continue
        canon = canonical_kmer(kmer)
        counts[canon] = counts.get(canon, 0) + 1
        order.append(canon)
    if not order:
        return 0.0
    unique_positions = sum(1 for c in order if counts[c] == 1)
    return unique_positions / len(order)


@dataclass(frozen=True)
class KSelection:
    k: int
    uniqueness: float
    representative: Path
    scanned: dict[int, float]


def _representative(fasta_paths: list[Path]) -> Path:
    """Pick the median-length genome as the representative."""
    sized = sorted(fasta_paths, key=lambda p: len(_concatenated_sequence(p)))
    return sized[len(sized) // 2]


def select_k(
    fasta_paths: list[Path],
    *,
    k_min: int = 11,
    k_max: int = 31,
    target: float = 0.99,
) -> KSelection:
    """Choose the smallest odd k whose uniqueness meets ``target``.

    Falls back to ``k_max`` (the most specific candidate) if none reaches the
    target.
    """
    rep = _representative(fasta_paths)
    seq = _concatenated_sequence(rep)
    start = k_min if k_min % 2 == 1 else k_min + 1
    scanned: dict[int, float] = {}
    chosen: int | None = None
    for k in range(start, k_max + 1, 2):
        u = kmer_uniqueness(seq, k)
        scanned[k] = u
        if chosen is None and u >= target:
            chosen = k
    if chosen is None:
        chosen = max(scanned) if scanned else start
    return KSelection(
        k=chosen, uniqueness=scanned.get(chosen, 0.0), representative=rep, scanned=scanned
    )
