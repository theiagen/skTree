"""Tests for report-data assembly and the embedded JSON payload."""

from __future__ import annotations

import numpy as np

from sktree.html_report import (
    MAX_ALL_SNP_COLUMNS,
    build_report_data,
    render_html,
    select_informative_columns,
    to_payload,
    write_html_report,
)
from sktree.report import build_stats
from sktree.snps import parse_alignment

ALN = ">s1\nAAAGA\n>s2\nACCGC\n>s3\nAA--T\n>s4\nACC--\n"


def _stats_and_matrix():
    m = parse_alignment(ALN)
    cls = m.classify()
    stats = build_stats(m, cls, k=31, min_freq=0.9, majority_threshold=0.5)
    return stats, m, cls


def test_select_informative_columns_caps_and_sorts():
    m = parse_alignment(">a\nACGTACGT\n>b\nAGGTTCGA\n>c\nACGAACGT\n")
    idx = select_informative_columns(m, cap=3)
    assert len(idx) == 3
    assert list(idx) == sorted(idx)


def test_build_report_data_caps_all_matrix():
    stats, m, cls = _stats_and_matrix()
    big = m.subset(np.ones(m.n_loci, dtype=bool))
    data = build_report_data(
        stats=stats, trees={"nj": "(s1,s2,(s3,s4));"},
        core_matrix=m.subset(cls.core), all_matrix=big,
    )
    # tiny here, so no truncation
    assert data.all_truncated_to is None
    assert data.all_matrix is not None


def test_to_payload_shape():
    stats, m, cls = _stats_and_matrix()
    data = build_report_data(
        stats=stats, trees={"nj": "(s1,s2,(s3,s4));"},
        core_matrix=m.subset(cls.core),
        clusters={"s1": 1, "s2": 1, "s3": 2, "s4": 2},
        annotation_summary={"coding": 3, "nonsynonymous": 1, "intergenic": 2},
        max_divergence=0.02, divergence_warn_threshold=0.01,
    )
    payload = to_payload(data)
    assert payload["samples"] == ["s1", "s2", "s3", "s4"]
    assert "nj" in payload["trees"]
    assert payload["alignment"]["core"]["columns"] == data.core_matrix.n_loci
    assert payload["alignment"]["all"] is None
    assert payload["clusters"] == {"s1": 1, "s2": 1, "s3": 2, "s4": 2}
    assert payload["annotation"]["coding"] == 3
    assert payload["divergence"]["max"] == 0.02
    assert set(payload["alignment"]["core"]["rows"]) == {"s1", "s2", "s3", "s4"}


def test_to_payload_carries_outputs():
    stats, m, cls = _stats_and_matrix()
    outputs = [
        {"name": "alignment.fasta", "description": "SNP alignment"},
        {"name": "summary.txt", "description": "Run summary"},
    ]
    data = build_report_data(
        stats=stats, trees={"nj": "(s1,s2,(s3,s4));"},
        core_matrix=m.subset(cls.core), outputs=outputs,
    )
    assert to_payload(data)["outputs"] == outputs


def test_to_payload_outputs_default_none():
    stats, m, cls = _stats_and_matrix()
    data = build_report_data(
        stats=stats, trees={"nj": "(s1,s2,(s3,s4));"},
        core_matrix=m.subset(cls.core),
    )
    assert to_payload(data)["outputs"] is None


def test_max_all_snp_columns_default():
    assert MAX_ALL_SNP_COLUMNS == 20_000


def _data():
    stats, m, cls = _stats_and_matrix()
    return build_report_data(
        stats=stats, trees={"nj": "(s1,s2,(s3,s4));"},
        core_matrix=m.subset(cls.core),
        clusters={"s1": 1, "s2": 1, "s3": 2, "s4": 2},
    )


def test_render_html_is_self_contained():
    html = render_html(_data())
    assert "<canvas id=\"aln\">" in html
    assert "SKTREE_RENDER" in html
    # No external resources: proves single-file portability.
    assert "src=\"http" not in html
    assert "href=\"http" not in html
    # Embedded data + tree present.
    assert "(s1,s2,(s3,s4));" in html
    assert "\"s3\": 2" in html or "\"s3\":2" in html


def test_render_html_graceful_without_clusters_or_annotation():
    stats, m, cls = _stats_and_matrix()
    data = build_report_data(stats=stats, trees={"nj": "(s1,s2,(s3,s4));"},
                             core_matrix=m.subset(cls.core))
    html = render_html(data)
    assert "\"clusters\": null" in html or "\"clusters\":null" in html
    assert "\"annotation\": null" in html or "\"annotation\":null" in html


def test_render_html_escapes_script_breakout():
    stats, m, cls = _stats_and_matrix()
    data = build_report_data(
        stats=stats,
        trees={"nj": "(evil</script><img src=x>:1,s2:1);"},
        core_matrix=m.subset(cls.core),
    )
    html = render_html(data)
    assert "</script><img src=x>" not in html
    assert "\\u003c/script\\u003e\\u003cimg src=x\\u003e" in html


def test_build_report_data_surfaces_truncation(monkeypatch):
    import sktree.html_report as hr

    monkeypatch.setattr(hr, "MAX_ALL_SNP_COLUMNS", 2)
    stats, m, cls = _stats_and_matrix()
    big = m.subset(np.ones(m.n_loci, dtype=bool))
    assert big.n_loci > 2
    data = hr.build_report_data(
        stats=stats, trees={"nj": "(s1,s2,(s3,s4));"},
        core_matrix=m.subset(cls.core), all_matrix=big,
    )
    assert data.all_truncated_to == 2
    assert data.all_matrix.n_loci == 2
    payload = hr.to_payload(data)
    assert payload["alignment"]["all"]["truncatedTo"] == 2
    assert payload["alignment"]["all"]["columns"] == 2


def test_write_html_report(tmp_path):
    out = tmp_path / "report.html"
    write_html_report(out, _data())
    assert out.exists() and out.stat().st_size > 1000
