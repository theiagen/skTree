"""Parse SKA's SNP alignment into a matrix and partition loci by presence.

SKA's ``align`` output is an *unordered* FASTA where every column is one SNP
locus (a split-k-mer with a variable middle base). We partition loci into
several categories:

- **core** SNPs   — variable loci genotyped in *every* sample.
- **majority** SNPs — variable loci present in at least a configurable fraction
  of samples (the soft-core notion).
- **all** SNPs    — every variable locus.

A base is "present" only if it is an unambiguous A/C/G/T; gaps (``-``), ``N`` and
any IUPAC ambiguity code count as missing.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

VALID_BASES = ("A", "C", "G", "T")


def _read_fasta(text: str) -> tuple[list[str], list[str]]:
    names: list[str] = []
    seqs: list[str] = []
    current: list[str] = []
    for line in text.splitlines():
        if not line:
            continue
        if line.startswith(">"):
            if names:
                seqs.append("".join(current))
            names.append(line[1:].split()[0])
            current = []
        else:
            current.append(line.strip())
    if names:
        seqs.append("".join(current))
    return names, seqs


@dataclass(frozen=True)
class LociClassification:
    """Boolean masks (length = n_loci) for each SNP category."""

    variable: np.ndarray
    core: np.ndarray
    majority: np.ndarray

    @property
    def snp(self) -> np.ndarray:
        """All SNPs == all variable loci."""
        return self.variable


@dataclass
class SnpMatrix:
    """A SNP alignment as a ``(n_samples, n_loci)`` array of single chars."""

    sample_names: list[str]
    matrix: np.ndarray  # dtype '<U1'

    @property
    def n_samples(self) -> int:
        return self.matrix.shape[0]

    @property
    def n_loci(self) -> int:
        return self.matrix.shape[1]

    # -- per-locus statistics -------------------------------------------------

    def _present_mask(self) -> np.ndarray:
        """Boolean ``(n_samples, n_loci)``: True where the base is unambiguous."""
        return np.isin(self.matrix, VALID_BASES)

    def present_fraction(self) -> np.ndarray:
        """Fraction of samples with an unambiguous base, per locus."""
        return self._present_mask().mean(axis=0)

    def is_variable(self) -> np.ndarray:
        """True per locus where >= 2 distinct unambiguous alleles occur."""
        present = self._present_mask()
        result = np.zeros(self.n_loci, dtype=bool)
        for j in range(self.n_loci):
            col = self.matrix[present[:, j], j]
            result[j] = np.unique(col).size >= 2
        return result

    # -- partitioning ---------------------------------------------------------

    def classify(self, majority_threshold: float = 0.5) -> LociClassification:
        variable = self.is_variable()
        frac = self.present_fraction()
        core = variable & (frac >= 1.0)
        majority = variable & (frac >= majority_threshold)
        return LociClassification(variable=variable, core=core, majority=majority)

    def subset(self, mask: np.ndarray) -> SnpMatrix:
        return SnpMatrix(sample_names=list(self.sample_names), matrix=self.matrix[:, mask])

    def core_alignment(self, majority_threshold: float = 0.5) -> SnpMatrix:
        return self.subset(self.classify(majority_threshold).core)

    def snp_alignment(self) -> SnpMatrix:
        return self.subset(self.is_variable())

    # -- output ---------------------------------------------------------------

    def to_fasta(self) -> str:
        rows = []
        for name, row in zip(self.sample_names, self.matrix, strict=True):
            rows.append(f">{name}\n{''.join(row)}")
        return "\n".join(rows) + "\n"


def parse_alignment(source: str | Path) -> SnpMatrix:
    """Parse a SKA SNP alignment (FASTA text or path) into a :class:`SnpMatrix`."""
    if isinstance(source, Path):
        text = source.read_text()
    elif "\n" not in source and source and Path(source).is_file():
        text = Path(source).read_text()
    else:
        text = source

    names, seqs = _read_fasta(text)
    if not names:
        raise ValueError("alignment contains no sequences")
    lengths = {len(s) for s in seqs}
    if len(lengths) != 1:
        raise ValueError(f"sequences have unequal length: {sorted(lengths)}")

    matrix = np.array([list(s.upper()) for s in seqs], dtype="<U1")
    return SnpMatrix(sample_names=names, matrix=matrix)
