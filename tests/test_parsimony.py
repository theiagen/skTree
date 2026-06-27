"""Tests for pure-Python maximum-parsimony tree inference (Fitch + search)."""

from __future__ import annotations

import dendropy

from sktree.parsimony import parsimony_tree, score_tree
from sktree.snps import parse_alignment


def _splits(tree: dendropy.Tree) -> set[frozenset[str]]:
    tree.encode_bipartitions()
    labels = {leaf.taxon.label for leaf in tree.leaf_node_iter()}
    out = set()
    for bip in tree.bipartition_encoding:
        side = frozenset(t.label for t in bip.leafset_taxa(tree.taxon_namespace))
        if 1 < len(side) < len(labels):
            out.add(side)
            out.add(frozenset(labels - side))
    return out


def test_fitch_score_known_value():
    # topology ((a,b),(c,d)); single site a=A b=A c=T d=T -> 1 change
    aln = ">a\nA\n>b\nA\n>c\nT\n>d\nT\n"
    m = parse_alignment(aln)
    tns = dendropy.TaxonNamespace()
    tree = dendropy.Tree.get(
        data="((a,b),(c,d));", schema="newick", taxon_namespace=tns
    )
    assert score_tree(tree, m) == 1
    # the cross topology needs 2 changes for this site
    tree2 = dendropy.Tree.get(
        data="((a,c),(b,d));", schema="newick", taxon_namespace=tns
    )
    assert score_tree(tree2, m) == 2


def test_parsimony_groups_identical_pairs():
    # s1==s3, s2==s4 -> MP must put each identical pair together
    aln = ">s1\nAAAAA\n>s2\nTTTTT\n>s3\nAAAAA\n>s4\nTTTTT\n"
    m = parse_alignment(aln)
    tree, score = parsimony_tree(m)
    assert score == 5  # one change per site on the central edge
    splits = _splits(tree)
    assert frozenset({"s1", "s3"}) in splits


def test_parsimony_missing_data_costs_nothing():
    aln = ">a\nA\n>b\nA\n>c\n-\n"  # c missing -> no extra cost
    m = parse_alignment(aln)
    tree, score = parsimony_tree(m)
    assert score == 0


def test_parsimony_score_no_worse_than_nj():
    from sktree.phylo import neighbor_joining

    aln = ">s1\nACGTACGT\n>s2\nACGAACGT\n>s3\nTCGTACGA\n>s4\nTCGTTCGA\n>s5\nACGTACGT\n"
    m = parse_alignment(aln)
    nj = neighbor_joining(m)
    _, p_score = parsimony_tree(m)
    assert p_score <= score_tree(nj, m)
