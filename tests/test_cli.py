"""CLI argument wiring tests."""

from __future__ import annotations

import shutil

import pytest

from sktree.cli import _build_parser

needs_ska = pytest.mark.skipif(shutil.which("ska") is None, reason="ska binary not installed")


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


def test_debug_flag_parses():
    parser = _build_parser()
    args = parser.parse_args(["run", "a.fa", "b.fa", "-o", "out", "--debug"])
    assert args.debug is True


def test_debug_defaults_false():
    parser = _build_parser()
    args = parser.parse_args(["run", "a.fa", "b.fa", "-o", "out"])
    assert args.debug is False


def test_failed_run_writes_log_and_points_to_it(tmp_path, capsys):
    # A run that fails (here: only one input) must still leave a log file behind
    # and tell the user where to find it, so the failure can be investigated.
    from sktree.cli import main

    out = tmp_path / "out"
    rc = main(["run", str(tmp_path / "only.fasta"), "-o", str(out)])

    assert rc == 1
    log_file = out / "sktree.log"
    assert log_file.exists()
    assert "two input" in log_file.read_text()
    # the user is told where the log lives
    assert "sktree.log" in capsys.readouterr().err


@needs_ska
def test_successful_run_logs_ska_commands(synthetic_genomes, tmp_path):
    # End-to-end: a clean run records the exact ska commands it executed, so the
    # log is a faithful trace even when nothing went wrong.
    from sktree.cli import main

    out = tmp_path / "out"
    rc = main(["run", *[str(p) for p in synthetic_genomes], "-o", str(out), "-k", "31"])
    assert rc == 0

    log_text = (out / "sktree.log").read_text()
    assert "running: " in log_text  # DEBUG command trace is always captured
    assert "ska" in log_text and "build" in log_text
    assert "Run completed successfully" in log_text
