"""Tests for distance computation and neighbor-joining tree inference."""

from __future__ import annotations

import dendropy
import numpy as np

from sktree.phylo import neighbor_joining, snp_distance_matrix
from sktree.snps import parse_alignment

ALN = """>s1
AAAGA
>s2
ACCGC
>s3
AA--T
>s4
ACC--
"""


def test_snp_distance_matrix_values():
    m = parse_alignment(ALN)
    names, d = snp_distance_matrix(m)
    assert names == ["s1", "s2", "s3", "s4"]
    idx = {n: i for i, n in enumerate(names)}

    def dist(a, b):
        return d[idx[a], idx[b]]

    # hand-computed over co-present A/C/G/T positions only
    assert dist("s1", "s2") == 3
    assert dist("s1", "s3") == 1
    assert dist("s1", "s4") == 2
    assert dist("s2", "s3") == 2
    assert dist("s2", "s4") == 0
    assert dist("s3", "s4") == 1
    # symmetry + zero diagonal
    assert np.array_equal(d, d.T)
    assert np.all(np.diag(d) == 0)


def test_neighbor_joining_returns_tree_with_all_taxa():
    m = parse_alignment(ALN)
    tree = neighbor_joining(m)
    labels = {leaf.taxon.label for leaf in tree.leaf_node_iter()}
    assert labels == {"s1", "s2", "s3", "s4"}
    newick = tree.as_string(schema="newick")
    # round-trips as valid newick
    reparsed = dendropy.Tree.get(data=newick, schema="newick")
    assert len(reparsed.leaf_nodes()) == 4


def test_two_sample_tree():
    m = parse_alignment(">a\nAC\n>b\nAT\n")
    tree = neighbor_joining(m)
    labels = {leaf.taxon.label for leaf in tree.leaf_node_iter()}
    assert labels == {"a", "b"}
