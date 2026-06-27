# tests/test_cluster.py
"""Unit tests for the fastbaps subprocess wrapper.

Mocks ``subprocess.run`` to assert the argv contract and CSV parsing without
needing the real ``fastbaps`` binary.
"""

from __future__ import annotations

import subprocess

import pytest

from sktree.engine.cluster import FastbapsError, FastbapsNotAvailable, FastbapsRunner

CSV = "Isolates,Clusters\ns1,1\ns2,1\ns3,2\n"


@pytest.fixture
def fake_run(monkeypatch):
    calls = []

    def _fake(argv, **kwargs):
        calls.append((argv, kwargs))
        # The wrapper passes -o <path>; emulate fastbaps writing that file.
        out = argv[argv.index("-o") + 1]
        from pathlib import Path

        Path(out).write_text(CSV)
        return subprocess.CompletedProcess(argv, 0, stdout="done", stderr="")

    monkeypatch.setattr("sktree.engine.cluster.shutil.which", lambda _: "/usr/bin/fastbaps")
    monkeypatch.setattr("sktree.engine.cluster.subprocess.run", _fake)
    return calls


def test_missing_binary_raises(monkeypatch):
    monkeypatch.setattr("sktree.engine.cluster.shutil.which", lambda _: None)
    with pytest.raises(FastbapsNotAvailable):
        FastbapsRunner()


def test_cluster_argv_and_parse(fake_run, tmp_path):
    runner = FastbapsRunner()
    aln = tmp_path / "core_snps.fasta"
    out = tmp_path / "clusters.csv"
    result = runner.cluster(aln, out, prior="baps", threads=4)
    argv = fake_run[0][0]
    assert argv[:1] == ["/usr/bin/fastbaps"]
    assert "-i" in argv and str(aln) in argv
    assert "-o" in argv and str(out) in argv
    assert "--prior" in argv and "baps" in argv
    assert "--threads" in argv and "4" in argv
    assert result == {"s1": 1, "s2": 1, "s3": 2}


def test_cluster_omits_threads_when_none(fake_run, tmp_path):
    runner = FastbapsRunner()
    runner.cluster(tmp_path / "a.fasta", tmp_path / "c.csv")
    argv = fake_run[0][0]
    assert "--threads" not in argv


def test_nonzero_exit_raises(monkeypatch, tmp_path):
    monkeypatch.setattr("sktree.engine.cluster.shutil.which", lambda _: "/usr/bin/fastbaps")

    def _boom(argv, **kwargs):
        return subprocess.CompletedProcess(argv, 1, stdout="", stderr="kaboom")

    monkeypatch.setattr("sktree.engine.cluster.subprocess.run", _boom)
    runner = FastbapsRunner()
    with pytest.raises(FastbapsError, match="kaboom"):
        runner.cluster(tmp_path / "a.fasta", tmp_path / "c.csv")
