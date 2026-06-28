"""CLI argument wiring tests."""

from __future__ import annotations

from sktree.cli import _build_parser


def test_run_parser_has_html_and_cluster_flags():
    parser = _build_parser()
    args = parser.parse_args(["run", "a.fa", "b.fa", "-o", "out", "--html", "--cluster"])
    assert args.html is True
    assert args.cluster is True


def test_html_and_cluster_default_false():
    parser = _build_parser()
    args = parser.parse_args(["run", "a.fa", "b.fa", "-o", "out"])
    assert args.html is False
    assert args.cluster is False


def test_read_filter_flags_parse():
    parser = _build_parser()
    args = parser.parse_args([
        "run", "a.fa", "S_R1.fq.gz", "S_R2.fq.gz", "-o", "out",
        "--min-count", "5", "--min-qual", "25", "--qual-filter", "middle",
    ])
    assert args.min_count == 5
    assert args.min_qual == 25
    assert args.qual_filter == "middle"


def test_min_count_defaults_to_three():
    parser = _build_parser()
    args = parser.parse_args(["run", "a.fa", "b.fa", "-o", "out"])
    assert args.min_count == 3
    assert args.qual_filter is None


def test_manifest_only_run_needs_no_positional():
    parser = _build_parser()
    args = parser.parse_args(["run", "-o", "out", "--manifest", "samples.tsv"])
    assert args.inputs == []
    assert str(args.manifest) == "samples.tsv"


def test_run_without_inputs_or_manifest_errors():
    import pytest

    from sktree.cli import main

    with pytest.raises(SystemExit):
        main(["run", "-o", "out"])
