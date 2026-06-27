"""Tests for the SNP matrix writer and summary formatting."""

from __future__ import annotations

from sktree.report import build_stats, format_summary, write_snp_matrix_tsv
from sktree.snps import parse_alignment

ALN = ">s1\nAAAGA\n>s2\nACCGC\n>s3\nAA--T\n>s4\nACC--\n"


def test_write_snp_matrix_tsv(tmp_path):
    m = parse_alignment(ALN)
    cls = m.classify()
    out = tmp_path / "snp_matrix.tsv"
    write_snp_matrix_tsv(out, m, cls)
    rows = out.read_text().splitlines()
    assert rows[0].split("\t") == ["locus", "category", "s1", "s2", "s3", "s4"]
    # one header + 5 loci
    assert len(rows) == 6
    categories = [r.split("\t")[1] for r in rows[1:]]
    assert categories[0] == "constant"  # locus 0 AAAA
    assert categories[1] == "core"      # locus 1 ACAC, present in all
    assert "majority" in categories


def test_format_summary_mentions_counts():
    m = parse_alignment(ALN)
    cls = m.classify()
    stats = build_stats(m, cls, k=31, min_freq=0.9, majority_threshold=0.5)
    text = format_summary(stats)
    assert "core SNPs" in text
    assert "samples              : 4" in text
    assert stats.n_core_snps == 1
