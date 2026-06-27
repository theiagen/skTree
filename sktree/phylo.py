"""Phylogenetic inference from a SNP alignment.

The MVP supports neighbor-joining (pure Python via DendroPy). Distances are SNP
counts computed only over positions where *both* samples have an unambiguous
base, which is the standard pairwise-deletion scheme for gappy SNP alignments.
"""

from __future__ import annotations

import io

import dendropy
import numpy as np

from .snps import VALID_BASES, SnpMatrix


def snp_distance_matrix(m: SnpMatrix) -> tuple[list[str], np.ndarray]:
    """Pairwise SNP-distance matrix (counts), pairwise-deletion of missing sites."""
    present = np.isin(m.matrix, VALID_BASES)
    n = m.n_samples
    dist = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            both = present[i] & present[j]
            diff = both & (m.matrix[i] != m.matrix[j])
            d = float(diff.sum())
            dist[i, j] = dist[j, i] = d
    return list(m.sample_names), dist


def _phylogenetic_distance_matrix(
    names: list[str], dist: np.ndarray
) -> dendropy.PhylogeneticDistanceMatrix:
    """Build a DendroPy PDM from a labelled distance matrix via its CSV reader."""
    header = "," + ",".join(names)
    lines = [header]
    for i, name in enumerate(names):
        row = ",".join(str(dist[i, j]) for j in range(len(names)))
        lines.append(f"{name},{row}")
    csv = "\n".join(lines) + "\n"
    return dendropy.PhylogeneticDistanceMatrix.from_csv(
        io.StringIO(csv),
        is_first_row_column_names=True,
        is_first_column_row_names=True,
        taxon_namespace=dendropy.TaxonNamespace(),
    )


def neighbor_joining(m: SnpMatrix) -> dendropy.Tree:
    """Infer a neighbor-joining tree from the SNP alignment."""
    names, dist = snp_distance_matrix(m)
    pdm = _phylogenetic_distance_matrix(names, dist)
    return pdm.nj_tree()


def to_newick(tree: dendropy.Tree) -> str:
    """Serialise a tree to Newick with a canonical (ladderized) rotation.

    Neighbor-joining tie-breaking can rotate sibling clades run-to-run without
    changing topology or branch lengths. Ladderizing sorts each node's children
    by descendant count, giving a stable, readable ordering so the same data
    always renders identically.
    """
    tree.ladderize(ascending=True)
    return tree.as_string(schema="newick").strip() + "\n"
