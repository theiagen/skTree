"""Integration tests exercising the real ``ska`` binary.

Skipped automatically when ska is not on PATH.
"""

from __future__ import annotations

import shutil

import pytest

from sktree.engine.ska import SkaRunner
from sktree.inputs import resolve_inputs

pytestmark = pytest.mark.skipif(shutil.which("ska") is None, reason="ska binary not installed")


def test_build_align_distance_roundtrip(synthetic_genomes, tmp_path):
    runner = SkaRunner()
    samples = resolve_inputs(synthetic_genomes, None)
    skf = runner.build(samples, out_prefix=tmp_path / "all", k=31)
    assert skf.exists()

    aln_path = runner.align(skf, output=tmp_path / "aln.fasta")
    text = aln_path.read_text()
    # three samples present in the alignment
    assert text.count(">") == 3
    for name in ("sample1", "sample2", "sample3"):
        assert f">{name}" in text

    dist_path = runner.distance(skf, output=tmp_path / "d.tsv")
    rows = [r for r in dist_path.read_text().splitlines() if r and not r.startswith("Sample1")]
    # 3 samples -> 3 unordered pairs
    assert len(rows) == 3


def test_align_to_stdout_is_clean_fasta(synthetic_genomes, tmp_path):
    runner = SkaRunner()
    samples = resolve_inputs(synthetic_genomes, None)
    skf = runner.build(samples, out_prefix=tmp_path / "all", k=31)
    text = runner.align(skf, output=None)
    # stdout must be pure FASTA: starts with '>' and carries no SKA banner/log noise
    assert text.lstrip().startswith(">")
    assert "SKA" not in text


def test_mixed_assembly_and_reads_recovers_snp(assembly_plus_reads, tmp_path):
    # A build mixing an assembly with paired-end reads must place both as samples
    # and recover the single planted SNP between them.
    runner = SkaRunner()
    samples = resolve_inputs(assembly_plus_reads["inputs"], None)
    assert {s.name for s in samples} == {"asm", "strainB"}
    assert any(s.is_reads for s in samples)  # reads were auto-paired

    skf = runner.build(samples, out_prefix=tmp_path / "mixed", k=31, min_count=3)
    text = runner.align(skf, output=tmp_path / "aln.fasta").read_text()

    records = {}
    name = None
    for line in text.splitlines():
        if line.startswith(">"):
            name = line[1:].strip()
            records[name] = ""
        elif name is not None:
            records[name] += line.strip()
    assert set(records) == {"asm", "strainB"}
    # both samples genotyped at every variable site, and the two differ (the SNP).
    cols = list(zip(records["asm"], records["strainB"], strict=True))
    assert cols, "alignment had no variable sites"
    assert any(a != b for a, b in cols)
    # the assembly carries the reference allele and reads carry the alt at the SNP.
    alleles = {a for a, _ in cols} | {b for _, b in cols}
    assert assembly_plus_reads["alt_base"] in alleles
