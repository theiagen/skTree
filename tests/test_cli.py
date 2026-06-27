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
