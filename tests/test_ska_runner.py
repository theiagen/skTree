"""Unit tests for the SKA2 subprocess wrapper.

These tests mock ``subprocess.run`` so they assert the *contract* (the argv we
hand to the ``ska`` binary and how we parse results) without needing the binary
installed. Real-binary checks live in the integration suite.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from sktree.engine.ska import (
    SkaError,
    SkaNotFoundError,
    SkaRunner,
    clean_sample_name,
)


@pytest.mark.parametrize(
    "path,expected",
    [
        ("/data/SAMD00052601.fa.gz", "SAMD00052601"),
        ("sample1.fasta", "sample1"),
        ("reads.fastq.gz", "reads"),
        ("/x/y/strain.fna", "strain"),
        ("E.coli.fa", "E.coli"),
        ("noext", "noext"),
    ],
)
def test_clean_sample_name(path, expected):
    assert clean_sample_name(path) == expected


@pytest.fixture
def fake_run(monkeypatch):
    """Capture the argv passed to subprocess.run and return a canned result."""
    calls = []

    def _fake(argv, **kwargs):
        calls.append((argv, kwargs))
        return subprocess.CompletedProcess(argv, 0, stdout="OUT", stderr="ERR")

    monkeypatch.setattr("sktree.engine.ska.shutil.which", lambda _: "/usr/bin/ska")
    monkeypatch.setattr("sktree.engine.ska.subprocess.run", _fake)
    return calls


def test_missing_binary_raises():
    with pytest.raises(SkaNotFoundError):
        SkaRunner(binary="definitely-not-a-real-binary-xyz")


def test_build_argv(fake_run, tmp_path, monkeypatch):
    runner = SkaRunner()
    files = [tmp_path / "a.fasta", tmp_path / "b.fasta"]
    captured: dict[str, str] = {}

    def _fake(argv, **kwargs):
        fake_run.append((argv, kwargs))
        # Read the file-list now; build() unlinks it in its finally clause.
        captured["filelist"] = Path(argv[argv.index("-f") + 1]).read_text()
        return subprocess.CompletedProcess(argv, 0, stdout="OUT", stderr="ERR")

    monkeypatch.setattr("sktree.engine.ska.subprocess.run", _fake)
    out = runner.build(files, out_prefix=tmp_path / "all", k=15)
    argv = fake_run[0][0]
    assert argv[:1] == ["/usr/bin/ska"]
    assert argv[1] == "build"
    assert "-o" in argv and str(tmp_path / "all") in argv
    assert "-k" in argv and "15" in argv
    # input files are passed via a clean-named file-list, not positionally.
    assert "-f" in argv
    assert f"a\t{files[0]}" in captured["filelist"]
    assert f"b\t{files[1]}" in captured["filelist"]
    assert out == tmp_path / "all.skf"


def test_build_rejects_even_k(fake_run, tmp_path):
    runner = SkaRunner()
    with pytest.raises(ValueError, match="odd"):
        runner.build([tmp_path / "a.fasta"], out_prefix=tmp_path / "x", k=16)


def test_build_single_strand_and_min_count(fake_run, tmp_path):
    runner = SkaRunner()
    runner.build(
        [tmp_path / "a.fastq"], out_prefix=tmp_path / "x", k=21,
        single_strand=True, min_count=4,
    )
    argv = fake_run[0][0]
    assert "--single-strand" in argv
    assert "--min-count" in argv and "4" in argv


def test_align_to_file(fake_run, tmp_path):
    runner = SkaRunner()
    skf = tmp_path / "all.skf"
    out = runner.align(skf, output=tmp_path / "aln.fasta", min_freq=0.8)
    argv = fake_run[0][0]
    assert argv[1] == "align"
    assert "-o" in argv and str(tmp_path / "aln.fasta") in argv
    assert "--min-freq" in argv and "0.8" in argv
    assert str(skf) in argv
    assert out == tmp_path / "aln.fasta"


def test_align_to_stdout_returns_text(fake_run, tmp_path):
    runner = SkaRunner()
    text = runner.align(tmp_path / "all.skf", output=None)
    # when no output path, stdout text is returned
    assert text == "OUT"
    argv = fake_run[0][0]
    assert "-o" not in argv


def test_distance_argv(fake_run, tmp_path):
    runner = SkaRunner()
    runner.distance(tmp_path / "all.skf", output=tmp_path / "d.tsv")
    argv = fake_run[0][0]
    assert argv[1] == "distance"
    assert str(tmp_path / "all.skf") in argv
    assert "-o" in argv and str(tmp_path / "d.tsv") in argv


def test_merge_argv(fake_run, tmp_path):
    runner = SkaRunner()
    skfs = [tmp_path / "a.skf", tmp_path / "b.skf"]
    out = runner.merge(skfs, out_prefix=tmp_path / "m")
    argv = fake_run[0][0]
    assert argv[1] == "merge"
    assert "-o" in argv and str(tmp_path / "m") in argv
    assert all(str(s) in argv for s in skfs)
    assert out == tmp_path / "m.skf"


def test_map_argv(fake_run, tmp_path):
    runner = SkaRunner()
    ref = tmp_path / "ref.fasta"
    skf = tmp_path / "all.skf"
    out = runner.map(ref, skf, output=tmp_path / "v.vcf", out_format="vcf")
    argv = fake_run[0][0]
    assert argv[1] == "map"
    assert "-f" in argv and "vcf" in argv
    # reference must come before the input on the ska command line
    assert argv.index(str(ref)) < argv.index(str(skf))
    assert "-o" in argv and str(tmp_path / "v.vcf") in argv
    assert out == tmp_path / "v.vcf"


def test_map_aln_argv(fake_run, tmp_path):
    runner = SkaRunner()
    ref = tmp_path / "ref.fasta"
    skf = tmp_path / "all.skf"
    out = runner.map(ref, skf, output=tmp_path / "ref.aln", out_format="aln", threads=4)
    argv = fake_run[0][0]
    assert argv[1] == "map"
    assert "-f" in argv and "aln" in argv
    assert "--threads" in argv and "4" in argv
    # reference must come before the input on the ska command line
    assert argv.index(str(ref)) < argv.index(str(skf))
    assert out == tmp_path / "ref.aln"


def test_map_rejects_bad_format(fake_run, tmp_path):
    runner = SkaRunner()
    with pytest.raises(ValueError, match="vcf"):
        runner.map(tmp_path / "ref.fasta", tmp_path / "all.skf", out_format="bogus")


def test_nonzero_exit_raises(monkeypatch, tmp_path):
    monkeypatch.setattr("sktree.engine.ska.shutil.which", lambda _: "/usr/bin/ska")

    def _boom(argv, **kwargs):
        return subprocess.CompletedProcess(argv, 2, stdout="", stderr="kaboom")

    monkeypatch.setattr("sktree.engine.ska.subprocess.run", _boom)
    runner = SkaRunner()
    with pytest.raises(SkaError, match="kaboom"):
        runner.distance(tmp_path / "all.skf")
