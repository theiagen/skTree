"""Unit tests for input resolution (positional auto-pairing + manifest).

These exercise the boundary that turns a user's CLI inputs into the list of
``Sample`` objects ``ska build`` consumes. No ``ska`` binary is needed.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from sktree.inputs import Sample, clean_sample_name, resolve_inputs


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


# -- positional: assemblies -------------------------------------------------


def test_positional_assemblies_one_sample_each():
    samples = resolve_inputs([Path("a.fasta"), Path("b.fa.gz")], None)
    assert samples == [
        Sample("a", (Path("a.fasta"),), is_reads=False),
        Sample("b", (Path("b.fa.gz"),), is_reads=False),
    ]


# -- positional: auto-pairing reads -----------------------------------------


@pytest.mark.parametrize(
    "fwd,rev,name",
    [
        ("B_R1.fastq.gz", "B_R2.fastq.gz", "B"),
        ("B_1.fq.gz", "B_2.fq.gz", "B"),
        ("SampleA_R1_001.fastq.gz", "SampleA_R2_001.fastq.gz", "SampleA"),
        ("x.R1.fastq", "x.R2.fastq", "x"),
    ],
)
def test_positional_auto_pairs_reads(fwd, rev, name):
    samples = resolve_inputs([Path(rev), Path(fwd)], None)  # order shouldn't matter
    assert len(samples) == 1
    s = samples[0]
    assert s.name == name
    assert s.is_reads is True
    # forward read is always first in the tuple
    assert s.files == (Path(fwd), Path(rev))


def test_lone_fastq_without_mate_errors():
    with pytest.raises(ValueError, match="mate"):
        resolve_inputs([Path("only_R1.fastq.gz")], None)


def test_fastq_with_unrecognized_naming_errors():
    with pytest.raises(ValueError, match="mate"):
        resolve_inputs([Path("weird.fastq.gz"), Path("other.fastq.gz")], None)


def test_mixed_positional_assembly_and_reads():
    samples = resolve_inputs(
        [Path("asm.fasta"), Path("S_R1.fastq.gz"), Path("S_R2.fastq.gz")], None
    )
    names = {s.name: s for s in samples}
    assert set(names) == {"asm", "S"}
    assert names["asm"].is_reads is False
    assert names["S"].is_reads is True
    assert names["S"].files == (Path("S_R1.fastq.gz"), Path("S_R2.fastq.gz"))


# -- manifest ----------------------------------------------------------------


def test_manifest_two_and_three_column(tmp_path):
    manifest = tmp_path / "samples.tsv"
    manifest.write_text(
        "# a comment\n"
        "\n"
        "strainA\tstrainA.fasta\n"
        "strainB\treads/B_1.fq.gz\treads/B_2.fq.gz\n"
    )
    samples = resolve_inputs([], manifest)
    assert samples == [
        Sample("strainA", (Path("strainA.fasta"),), is_reads=False),
        Sample("strainB", (Path("reads/B_1.fq.gz"), Path("reads/B_2.fq.gz")), is_reads=True),
    ]


def test_manifest_bad_column_count_errors(tmp_path):
    manifest = tmp_path / "bad.tsv"
    manifest.write_text("onlyname\n")
    with pytest.raises(ValueError, match="column"):
        resolve_inputs([], manifest)


def test_manifest_and_positional_combine(tmp_path):
    manifest = tmp_path / "samples.tsv"
    manifest.write_text("m1\tm1.fasta\n")
    samples = resolve_inputs([Path("p1.fasta")], manifest)
    assert {s.name for s in samples} == {"m1", "p1"}


# -- name de-duplication -----------------------------------------------------


def test_duplicate_names_disambiguated():
    samples = resolve_inputs([Path("dir1/x.fasta"), Path("dir2/x.fasta")], None)
    names = [s.name for s in samples]
    assert names[0] == "x"
    assert names[1] != "x"
    assert len(set(names)) == 2
