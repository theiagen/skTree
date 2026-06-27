"""Tests for optimal odd-k selection (Kchooser-style)."""

from __future__ import annotations

import random

import pytest

from sktree.kselect import canonical_kmer, kmer_uniqueness, select_k


def test_canonical_kmer_picks_lexicographic_min_of_strand():
    # AAA vs its rc TTT -> AAA is smaller
    assert canonical_kmer("AAA") == "AAA"
    # TTT -> rc AAA
    assert canonical_kmer("TTT") == "AAA"
    # a palindromic-ish case
    assert canonical_kmer("ACG") == canonical_kmer("CGT")


def test_uniqueness_all_distinct_is_one():
    rng = random.Random(0)
    seq = "".join(rng.choice("ACGT") for _ in range(2000))
    # at large k, a random sequence has essentially all-unique k-mers
    assert kmer_uniqueness(seq, 31) == pytest.approx(1.0)


def test_uniqueness_low_for_repetitive_small_k():
    seq = "ACGT" * 200  # highly repetitive
    u_small = kmer_uniqueness(seq, 3)
    assert u_small < 0.2


def test_uniqueness_increases_with_k():
    rng = random.Random(2)
    seq = "".join(rng.choice("ACGT") for _ in range(1500))
    assert kmer_uniqueness(seq, 7) <= kmer_uniqueness(seq, 21) + 1e-9


def test_select_k_returns_odd_in_bounds(tmp_path):
    rng = random.Random(3)
    seq = "".join(rng.choice("ACGT") for _ in range(3000))
    p = tmp_path / "g.fasta"
    p.write_text(f">g\n{seq}\n")
    sel = select_k([p], k_min=9, k_max=31, target=0.99)
    assert sel.k % 2 == 1
    assert 9 <= sel.k <= 31
    assert sel.uniqueness >= 0.99 or sel.k == 31


def test_select_k_skips_non_acgt(tmp_path):
    p = tmp_path / "g.fasta"
    p.write_text(">g\n" + "ACGTN" * 100 + "\n")
    sel = select_k([p], k_min=9, k_max=21)
    assert sel.k % 2 == 1
