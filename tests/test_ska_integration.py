"""Integration tests exercising the real ``ska`` binary.

Skipped automatically when ska is not on PATH.
"""

from __future__ import annotations

import shutil

import pytest

from sktree.engine.ska import SkaRunner

pytestmark = pytest.mark.skipif(shutil.which("ska") is None, reason="ska binary not installed")


def test_build_align_distance_roundtrip(synthetic_genomes, tmp_path):
    runner = SkaRunner()
    skf = runner.build(synthetic_genomes, out_prefix=tmp_path / "all", k=31)
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
    skf = runner.build(synthetic_genomes, out_prefix=tmp_path / "all", k=31)
    text = runner.align(skf, output=None)
    # stdout must be pure FASTA: starts with '>' and carries no SKA banner/log noise
    assert text.lstrip().startswith(">")
    assert "SKA" not in text
