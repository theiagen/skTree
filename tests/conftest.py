"""Shared fixtures and helpers for the skTree test suite."""

from __future__ import annotations

import gzip
import random
import shutil

import pytest

SKA_AVAILABLE = shutil.which("ska") is not None
requires_ska = pytest.mark.skipif(not SKA_AVAILABLE, reason="ska binary not installed")

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def write_fasta(path, name: str, seq: str) -> None:
    path.write_text(f">{name}_contig1\n{seq}\n")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def write_paired_reads(
    fwd_path, rev_path, seq: str, *, read_len: int = 150, coverage: int = 40, seed: int = 7
) -> None:
    """Tile deterministic paired-end reads across ``seq`` into gzipped FASTQ.

    Forward reads are taken straight from ``seq``; their mates are the reverse
    complement of the same fragment. Every base gets a high, constant quality so
    the quality filter never drops them; recovery then depends only on coverage
    clearing ``--min-count``.
    """
    rng = random.Random(seed)
    n_reads = max(1, len(seq) * coverage // read_len)
    qual = "I" * read_len
    fwd_records, rev_records = [], []
    for i in range(n_reads):
        start = rng.randint(0, max(0, len(seq) - read_len))
        frag = seq[start : start + read_len]
        fwd_records.append(f"@read{i}/1\n{frag}\n+\n{qual[: len(frag)]}\n")
        rc = _revcomp(frag)
        rev_records.append(f"@read{i}/2\n{rc}\n+\n{qual[: len(rc)]}\n")
    with gzip.open(fwd_path, "wt") as fh:
        fh.write("".join(fwd_records))
    with gzip.open(rev_path, "wt") as fh:
        fh.write("".join(rev_records))


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


@pytest.fixture
def assembly_plus_reads(tmp_path):
    """One assembly and one paired-read sample differing by a single known SNP.

    ``asm`` is a 2 kb genome written as FASTA; ``reads`` is the same genome with
    one base changed at ``snp_pos``, emitted as paired-end FASTQ at high
    coverage. A mixed build must recover exactly this difference. Returns a dict
    with the positional input paths plus the planted SNP for assertions.
    """
    rng = random.Random(3)
    base = "".join(rng.choice("ACGT") for _ in range(2000))
    snp_pos = 1000
    alt = "A" if base[snp_pos] != "A" else "T"
    variant = base[:snp_pos] + alt + base[snp_pos + 1 :]

    asm = tmp_path / "asm.fasta"
    write_fasta(asm, "asm", base)
    fwd = tmp_path / "strainB_R1.fastq.gz"
    rev = tmp_path / "strainB_R2.fastq.gz"
    write_paired_reads(fwd, rev, variant)

    return {
        "inputs": [asm, fwd, rev],
        "assembly": asm,
        "reads": (fwd, rev),
        "snp_pos": snp_pos,
        "ref_base": base[snp_pos],
        "alt_base": alt,
    }
