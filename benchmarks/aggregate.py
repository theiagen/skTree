#!/usr/bin/env python3
"""Aggregate per-cohort skTree-vs-kSNP4 benchmark cells into one results table.

Reads the directory tree produced by ``bench.sh`` (one subdir per cohort, each
holding ``ksnp4.metrics.tsv``, ``sktree.metrics.tsv``, the kSNP4 ``ksnp_out/`` and
the skTree ``sktree_out/``) and emits a tidy TSV plus a Markdown summary covering
runtime, peak memory, core-SNP concordance, and tree topology agreement (RF).
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import dendropy
from dendropy.calculate import treecompare


def _read_metrics(path: Path) -> dict[str, str]:
    lines = path.read_text().splitlines()
    if len(lines) < 2:
        return {}
    keys = lines[0].split("\t")
    vals = lines[1].split("\t")
    return dict(zip(keys, vals, strict=False))


def _ksnp_core_snps(cohort: Path) -> int | None:
    f = cohort / "ksnp_out" / "COUNT_coreSNPs"
    if not f.exists():
        return None
    m = re.search(r"Number core SNPs:\s*(\d+)", f.read_text())
    return int(m.group(1)) if m else None


def _sktree_core_snps(cohort: Path) -> int | None:
    f = cohort / "sktree_out" / "summary.txt"
    if not f.exists():
        return None
    m = re.search(r"core SNPs\s*:\s*(\d+)", f.read_text())
    return int(m.group(1)) if m else None


def _rf(tree_a: Path, tree_b: Path) -> tuple[int, int] | None:
    # A tree file can be missing or empty (e.g. kSNP4 writes a 0-byte core tree
    # when it finds no core SNPs); treat any parse failure as "no comparison".
    if not (tree_a.exists() and tree_b.exists()):
        return None
    if tree_a.stat().st_size == 0 or tree_b.stat().st_size == 0:
        return None
    taxa = dendropy.TaxonNamespace()
    try:
        a = dendropy.Tree.get(path=str(tree_a), schema="newick",
                              taxon_namespace=taxa, preserve_underscores=True)
        b = dendropy.Tree.get(path=str(tree_b), schema="newick",
                              taxon_namespace=taxa, preserve_underscores=True)
    except Exception:
        return None
    a.encode_bipartitions()
    b.encode_bipartitions()
    n = len(a.leaf_nodes())
    return treecompare.symmetric_difference(a, b), 2 * (n - 3)


def _fmt_gb(b: str) -> str:
    try:
        return f"{int(b) / 1e9:.2f}"
    except (ValueError, TypeError):
        return "NA"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("bench_dir", type=Path)
    ap.add_argument("-o", "--out", type=Path, default=None, help="write Markdown here")
    args = ap.parse_args()

    rows = []
    for cohort in sorted(p for p in args.bench_dir.iterdir() if p.is_dir()):
        ks = _read_metrics(cohort / "ksnp4.metrics.tsv")
        sk = _read_metrics(cohort / "sktree.metrics.tsv")
        if not ks or not sk:
            continue
        ks_core = _ksnp_core_snps(cohort)
        sk_core = _sktree_core_snps(cohort)
        rf = _rf(cohort / "sktree_out" / "tree_nj.nwk",
                 cohort / "ksnp_out" / "tree.core_SNPs.parsimony.tre")
        rows.append({
            "cohort": cohort.name,
            "n": ks.get("n_genomes", "?"),
            "ksnp_wall": ks.get("wall_seconds", "NA"),
            "sktree_wall": sk.get("wall_seconds", "NA"),
            "ksnp_gb": _fmt_gb(ks.get("peak_rss_bytes", "")),
            "sktree_gb": _fmt_gb(sk.get("peak_rss_bytes", "")),
            "ksnp_core": ks_core,
            "sktree_core": sk_core,
            "rf": f"{rf[0]}/{rf[1]}" if rf else "NA",
        })

    def _speedup(r):
        try:
            return f"{float(r['ksnp_wall']) / float(r['sktree_wall']):.1f}x"
        except (ValueError, ZeroDivisionError):
            return "NA"

    def _memratio(r):
        try:
            return f"{float(r['ksnp_gb']) / float(r['sktree_gb']):.1f}x"
        except (ValueError, ZeroDivisionError):
            return "NA"

    header = ("| cohort | n | kSNP4 s | skTree s | speedup | kSNP4 GB | skTree GB | "
              "mem ratio | kSNP4 core | skTree core | RF (NJ vs core) |")
    sep = "|" + "---|" * 11
    lines = [header, sep]
    for r in rows:
        lines.append(
            f"| {r['cohort']} | {r['n']} | {r['ksnp_wall']} | {r['sktree_wall']} | "
            f"{_speedup(r)} | {r['ksnp_gb']} | {r['sktree_gb']} | {_memratio(r)} | "
            f"{r['ksnp_core']} | {r['sktree_core']} | {r['rf']} |"
        )
    table = "\n".join(lines)
    print(table)
    if args.out:
        args.out.write_text(table + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
