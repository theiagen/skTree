"""Tests for SNP-alignment parsing and core/majority partitioning."""

from __future__ import annotations

import numpy as np
import pytest

from sktree.snps import parse_alignment

# A hand-built alignment with 4 samples and 5 loci.
#   locus 0: A A A A  -> constant (not a SNP)
#   locus 1: A C A C  -> variable, present in all   -> CORE SNP
#   locus 2: A C - C  -> variable, present in 3/4    -> majority (>=0.5) not core
#   locus 3: G G - -   -> present in 2/4, single allele -> not variable
#   locus 4: A C T -   -> variable, present in 3/4, multi-allelic
ALN = """>s1
AAAGA
>s2
ACCGC
>s3
AA--T
>s4
ACC--
"""


def test_parse_basic_shape():
    m = parse_alignment(ALN)
    assert m.sample_names == ["s1", "s2", "s3", "s4"]
    assert m.n_samples == 4
    assert m.n_loci == 5
    assert m.matrix.shape == (4, 5)


def test_parse_from_path(tmp_path):
    p = tmp_path / "aln.fasta"
    p.write_text(ALN)
    m = parse_alignment(p)
    assert m.n_samples == 4


def test_present_fraction_counts_gaps_and_n_as_missing():
    m = parse_alignment(ALN)
    frac = m.present_fraction()
    # locus 2 (index 2): A C - C -> 3/4 present
    assert frac[2] == pytest.approx(0.75)
    # locus 3 (index 3): G G - - -> 2/4 present
    assert frac[3] == pytest.approx(0.5)


def test_variable_mask():
    m = parse_alignment(ALN)
    var = m.is_variable()
    assert list(var) == [False, True, True, False, True]


def test_classify_core_and_majority():
    m = parse_alignment(ALN)
    cls = m.classify(majority_threshold=0.5)
    # core = variable AND present in all samples -> only locus 1
    assert list(np.where(cls.core)[0]) == [1]
    # majority = variable AND present fraction >= 0.5 -> loci 1, 2, 4
    assert list(np.where(cls.majority)[0]) == [1, 2, 4]
    # all SNPs = variable -> 1, 2, 4
    assert list(np.where(cls.snp)[0]) == [1, 2, 4]


def test_core_alignment_subset():
    m = parse_alignment(ALN)
    core = m.core_alignment()
    assert core.n_loci == 1
    # the single core locus is column 1: A C A C
    assert "".join(core.matrix[:, 0]) == "ACAC"
    assert core.sample_names == m.sample_names


def test_missing_char_treated_as_missing():
    aln = ">a\nACN\n>b\nACT\n"
    m = parse_alignment(aln)
    frac = m.present_fraction()
    assert frac[2] == pytest.approx(0.5)  # N is missing for sample a


def test_n_loci_zero_alignment():
    aln = ">a\nAAA\n>b\nAAA\n"  # all constant
    m = parse_alignment(aln)
    cls = m.classify()
    assert cls.snp.sum() == 0
    assert m.core_alignment().n_loci == 0


def test_unequal_lengths_raise():
    with pytest.raises(ValueError, match="length"):
        parse_alignment(">a\nACGT\n>b\nACG\n")


def test_empty_alignment_raises():
    with pytest.raises(ValueError, match="no sequences"):
        parse_alignment("")
