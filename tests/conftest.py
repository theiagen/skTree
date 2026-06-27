"""Shared fixtures and helpers for the skTree test suite."""

from __future__ import annotations

import random
import shutil

import pytest

SKA_AVAILABLE = shutil.which("ska") is not None
requires_ska = pytest.mark.skipif(not SKA_AVAILABLE, reason="ska binary not installed")


def write_fasta(path, name: str, seq: str) -> None:
    path.write_text(f">{name}_contig1\n{seq}\n")


@pytest.fixture
def synthetic_genomes(tmp_path):
    """Three 300 bp genomes sharing flanks with known SNPs at two positions.

    sample1 = reference; sample2 differs at pos 100 and 200; sample3 at pos 100
    only. Returns the list of FASTA paths.
    """
    rng = random.Random(1)
    base = "".join(rng.choice("ACGT") for _ in range(300))

    def mut(seq, changes):
        s = list(seq)
        for pos, alt in changes.items():
            s[pos] = alt
        return "".join(s)

    alt100 = "A" if base[100] != "A" else "T"
    alt200 = "G" if base[200] != "G" else "T"
    genomes = {
        "sample1": base,
        "sample2": mut(base, {100: alt100, 200: alt200}),
        "sample3": mut(base, {100: ("T" if base[100] != "T" else "G")}),
    }
    paths = []
    for name, seq in genomes.items():
        p = tmp_path / f"{name}.fasta"
        write_fasta(p, name, seq)
        paths.append(p)
    return paths
