"""Tests for the optional ML engine wrapper (argv construction + skip contract).

Neither IQ-TREE nor RAxML-NG is assumed installed, so engine detection and the
subprocess call are mocked. These tests pin the command line we generate and the
"no engine -> raise" behaviour the pipeline relies on to skip ML cleanly.
"""

from __future__ import annotations

from pathlib import Path
from unittest import mock

import pytest

from sktree.engine import ml
from sktree.engine.ml import MlError, MlNotAvailable, MlTreeBuilder


def _fake_which(mapping):
    return lambda name: mapping.get(name)


def test_no_engine_raises():
    with mock.patch.object(ml.shutil, "which", _fake_which({})):
        with pytest.raises(MlNotAvailable):
            MlTreeBuilder()


def test_auto_prefers_iqtree():
    both = {"iqtree2": "/usr/bin/iqtree2", "raxml-ng": "/usr/bin/raxml-ng"}
    with mock.patch.object(ml.shutil, "which", _fake_which(both)):
        b = MlTreeBuilder()
    assert b.engine == "iqtree"
    assert b.binary == "/usr/bin/iqtree2"


def test_auto_falls_back_to_raxml():
    with mock.patch.object(ml.shutil, "which", _fake_which({"raxml-ng": "/usr/bin/raxml-ng"})):
        b = MlTreeBuilder()
    assert b.engine == "raxml"


def test_force_iqtree_absent_raises():
    with mock.patch.object(ml.shutil, "which", _fake_which({"raxml-ng": "/x"})):
        with pytest.raises(MlNotAvailable):
            MlTreeBuilder(engine="iqtree")


def test_iqtree_argv_and_treefile():
    with mock.patch.object(ml.shutil, "which", _fake_which({"iqtree2": "/usr/bin/iqtree2"})):
        b = MlTreeBuilder()
    argv, treefile = b._argv(Path("aln.fasta"), Path("/out/ml"), "GTR+G", 1000, 4)
    assert argv[0] == "/usr/bin/iqtree2"
    assert "-s" in argv and "aln.fasta" in argv
    assert "-m" in argv and "GTR+G" in argv
    assert "-B" in argv and "1000" in argv
    assert "-nt" in argv and "4" in argv
    assert treefile == Path("/out/ml.treefile")


def test_iqtree_default_omits_model_for_modelfinder():
    # model=None -> no -m flag, so IQ-TREE runs its default ModelFinder Plus.
    with mock.patch.object(ml.shutil, "which", _fake_which({"iqtree2": "/i"})):
        b = MlTreeBuilder()
    argv, _ = b._argv(Path("aln.fasta"), Path("/out/ml"), None, None, None)
    assert "-m" not in argv
    assert "-s" in argv and "--prefix" in argv
    assert "-nt" in argv and "AUTO" in argv  # threads default
    assert "-B" not in argv  # no bootstrap by default


def test_raxml_default_model_fallback():
    # raxml-ng has no auto-search default, so model=None falls back to GTR+G.
    with mock.patch.object(ml.shutil, "which", _fake_which({"raxml-ng": "/r"})):
        b = MlTreeBuilder()
    argv, _ = b._argv(Path("aln.fasta"), Path("/out/ml"), None, None, None)
    assert "--model" in argv and "GTR+G" in argv


def test_iqtree_bootstrap_floored_to_1000():
    with mock.patch.object(ml.shutil, "which", _fake_which({"iqtree2": "/i"})):
        b = MlTreeBuilder()
    argv, _ = b._argv(Path("a"), Path("p"), "GTR+G", 100, None)
    assert "1000" in argv  # IQ-TREE ufboot requires >= 1000 replicates
    assert "AUTO" in argv  # no threads -> -nt AUTO


def test_raxml_argv_and_treefile():
    with mock.patch.object(ml.shutil, "which", _fake_which({"raxml-ng": "/usr/bin/raxml-ng"})):
        b = MlTreeBuilder()
    argv, treefile = b._argv(Path("aln.fasta"), Path("/out/ml"), "GTR+G", 200, 2)
    assert argv[0] == "/usr/bin/raxml-ng"
    assert "--all" in argv and "--msa" in argv
    assert "--bs-trees" in argv and "200" in argv
    assert "--threads" in argv and "2" in argv
    assert treefile == Path("/out/ml.raxml.bestTree")


def test_build_returns_treefile_on_success():
    with mock.patch.object(ml.shutil, "which", _fake_which({"iqtree2": "/i"})):
        b = MlTreeBuilder()
    completed = mock.Mock(returncode=0, stdout="", stderr="")
    with mock.patch.object(ml.subprocess, "run", return_value=completed) as run:
        out = b.build(Path("aln.fasta"), Path("/out/ml"))
    run.assert_called_once()
    assert out == Path("/out/ml.treefile")


def test_build_raises_on_nonzero():
    with mock.patch.object(ml.shutil, "which", _fake_which({"iqtree2": "/i"})):
        b = MlTreeBuilder()
    completed = mock.Mock(returncode=2, stdout="", stderr="boom")
    with mock.patch.object(ml.subprocess, "run", return_value=completed):
        with pytest.raises(MlError, match="boom"):
            b.build(Path("aln.fasta"), Path("/out/ml"))
