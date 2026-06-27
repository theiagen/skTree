"""Textual report and SNP-matrix writer."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .snps import LociClassification, SnpMatrix


@dataclass
class RunStats:
    n_samples: int
    k: int
    min_freq: float
    majority_threshold: float
    n_all_snps: int
    n_core_snps: int
    n_majority_snps: int
    sample_names: list[str]


def write_snp_matrix_tsv(path: Path, m: SnpMatrix, cls: LociClassification) -> None:
    """Write loci-by-sample SNP matrix with a category label per locus."""
    header = ["locus", "category", *m.sample_names]
    lines = ["\t".join(header)]
    for j in range(m.n_loci):
        if cls.core[j]:
            cat = "core"
        elif cls.majority[j]:
            cat = "majority"
        elif cls.variable[j]:
            cat = "snp"
        else:
            cat = "constant"
        bases = list(m.matrix[:, j])
        lines.append("\t".join([f"locus_{j}", cat, *bases]))
    path.write_text("\n".join(lines) + "\n")


def build_stats(
    m: SnpMatrix, cls: LociClassification, *, k: int, min_freq: float, majority_threshold: float
) -> RunStats:
    return RunStats(
        n_samples=m.n_samples,
        k=k,
        min_freq=min_freq,
        majority_threshold=majority_threshold,
        n_all_snps=int(cls.snp.sum()),
        n_core_snps=int(cls.core.sum()),
        n_majority_snps=int(cls.majority.sum()),
        sample_names=list(m.sample_names),
    )


def format_summary(stats: RunStats) -> str:
    lines = [
        "skTree run summary",
        "==================",
        f"samples              : {stats.n_samples}",
        f"k-mer size           : {stats.k}",
        f"align min-freq       : {stats.min_freq}",
        f"majority threshold   : {stats.majority_threshold}",
        "",
        "SNP categories",
        "--------------",
        f"all SNPs             : {stats.n_all_snps}",
        f"core SNPs            : {stats.n_core_snps}",
        f"majority SNPs        : {stats.n_majority_snps}",
        "",
        "samples: " + ", ".join(stats.sample_names),
    ]
    return "\n".join(lines) + "\n"
