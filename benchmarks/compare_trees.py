#!/usr/bin/env python3
"""Compare two phylogenetic trees by Robinson-Foulds distance.

Both tools (skTree and kSNP4) label leaves by the input genome's filename stem,
so the taxon namespaces line up directly. We report the symmetric RF distance and
its normalised form (RF / max possible RF) — 0.0 means identical topology, 1.0
means maximally different. RF is the standard yardstick for "did two methods
recover the same tree?".
"""

from __future__ import annotations

import argparse
import sys

import dendropy
from dendropy.calculate import treecompare


def _load(path: str, taxa: dendropy.TaxonNamespace) -> dendropy.Tree:
    # Trees from different tools must share ONE taxon namespace for RF to work.
    return dendropy.Tree.get(
        path=path,
        schema="newick",
        taxon_namespace=taxa,
        preserve_underscores=True,
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("tree_a")
    ap.add_argument("tree_b")
    ap.add_argument("--label-a", default="tree_a")
    ap.add_argument("--label-b", default="tree_b")
    args = ap.parse_args()

    taxa = dendropy.TaxonNamespace()
    a = _load(args.tree_a, taxa)
    b = _load(args.tree_b, taxa)

    leaves_a = {lf.taxon.label for lf in a.leaf_node_iter()}
    leaves_b = {lf.taxon.label for lf in b.leaf_node_iter()}
    shared = leaves_a & leaves_b
    if leaves_a != leaves_b:
        print(
            f"WARNING: leaf sets differ ({args.label_a}={len(leaves_a)}, "
            f"{args.label_b}={len(leaves_b)}, shared={len(shared)})",
            file=sys.stderr,
        )

    a.encode_bipartitions()
    b.encode_bipartitions()
    rf = treecompare.symmetric_difference(a, b)
    # Max RF for two fully-resolved unrooted trees on n taxa is 2*(n-3).
    n = len(shared)
    max_rf = 2 * (n - 3) if n > 3 else 0
    norm = rf / max_rf if max_rf else 0.0

    print(f"taxa\t{n}")
    print(f"rf_distance\t{rf}")
    print(f"max_rf\t{max_rf}")
    print(f"normalised_rf\t{norm:.4f}")
    print(f"identical_topology\t{rf == 0}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
