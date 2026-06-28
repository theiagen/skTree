"""Command-line interface for skTree."""

from __future__ import annotations

import argparse
import logging
import shlex
import sys
from pathlib import Path

from . import __version__
from .engine.ska import SkaError
from .logconfig import configure_logging
from .pipeline import run_pipeline

logger = logging.getLogger("sktree")


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="sktree",
        description=(
            "Reference-free, alignment-free SNP discovery and phylogenetics "
            "powered by the SKA2 split-k-mer engine."
        ),
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="run the full pipeline: inputs -> SNPs -> tree")
    run.add_argument(
        "inputs", nargs="*", type=Path,
        help="FASTA assemblies and/or paired FASTQ reads. Paired reads are "
             "auto-detected by _R1/_R2, _1/_2 or .R1/.R2 naming.",
    )
    run.add_argument(
        "--manifest", type=Path, default=None,
        help="TSV sample sheet: 'name<TAB>file' (assembly/single) or "
             "'name<TAB>fwd<TAB>rev' (paired reads) per line. Combines with "
             "positional inputs.",
    )
    run.add_argument("-o", "--outdir", type=Path, required=True, help="output directory")
    run.add_argument("-k", type=int, default=31, help="odd k-mer size (default: 31)")
    run.add_argument(
        "--auto-k", action="store_true",
        help="auto-select k (Kchooser-style) instead of using -k",
    )
    run.add_argument(
        "-m", "--min-freq", type=float, default=0.9,
        help="min fraction of samples a k-mer must appear in (default: 0.9)",
    )
    run.add_argument(
        "--majority-threshold", type=float, default=0.5,
        help="present-fraction cutoff for a 'majority' SNP (default: 0.5)",
    )
    run.add_argument(
        "--min-count", type=int, default=3,
        help="minimum k-mer count for read samples; filters sequencing error "
             "(default: 3; ignored for assemblies)",
    )
    run.add_argument(
        "--min-qual", type=int, default=None,
        help="minimum base quality for read samples (ska default: 20)",
    )
    run.add_argument(
        "--qual-filter", choices=("no-filter", "middle", "strict"), default=None,
        help="read quality-filtering strategy (ska default: strict)",
    )
    run.add_argument(
        "--parsimony", action="store_true",
        help="also infer a maximum-parsimony tree (pure Python; slow past ~30 samples)",
    )
    run.add_argument(
        "--ml", action="store_true",
        help="also infer a maximum-likelihood tree (needs IQ-TREE or RAxML-NG on PATH)",
    )
    run.add_argument(
        "--reference", type=Path, default=None,
        help="reference FASTA to annotate SNPs by gene (maps SNPs + predicts genes via pyrodigal)",
    )
    run.add_argument(
        "--map-tree", action="store_true",
        help="build reference-anchored tree(s) from the ska map pseudo-alignment "
             "(requires --reference; honors --ml/--parsimony)",
    )
    run.add_argument(
        "--cluster", action="store_true",
        help="assign population-structure clusters with fastbaps (writes clusters.csv)",
    )
    run.add_argument(
        "--html", action="store_true",
        help="write a self-contained interactive HTML report (runs fastbaps if available)",
    )
    run.add_argument("--threads", type=int, default=None, help="CPU threads for SKA")
    run.add_argument("-v", "--verbose", action="store_true", help="verbose logging (INFO)")
    run.add_argument(
        "--debug", action="store_true",
        help="very verbose console logging (DEBUG); a full DEBUG log is always "
             "written to <outdir>/sktree.log regardless of this flag",
    )
    return parser


def _log_run_header(argv: list[str] | None, args: argparse.Namespace) -> None:
    """Record what was run, so a log read in isolation explains the context."""
    invocation = sys.argv if argv is None else ["sktree", *argv]
    logger.info("skTree %s", __version__)
    logger.debug("command: %s", shlex.join(str(a) for a in invocation))
    logger.debug(
        "params: k=%s auto_k=%s min_freq=%s min_count=%s min_qual=%s "
        "qual_filter=%s threads=%s",
        args.k, args.auto_k, args.min_freq, args.min_count, args.min_qual,
        args.qual_filter, args.threads,
    )
    n_inputs = len(args.inputs)
    logger.debug(
        "inputs: %d positional file(s)%s -> %s",
        n_inputs,
        f", manifest={args.manifest}" if args.manifest else "",
        [str(p) for p in args.inputs],
    )
    try:
        from .engine.ska import SkaRunner

        logger.debug("engine: %s", SkaRunner().version())
    except SkaError as exc:
        # Don't fail here; run_pipeline will report the missing binary properly.
        logger.debug("engine: ska version unavailable (%s)", exc)


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.command == "run" and not args.inputs and args.manifest is None:
        parser.error("provide input files and/or --manifest")

    if getattr(args, "map_tree", False) and args.reference is None:
        parser.error("--map-tree requires --reference")

    if args.command == "run":
        log_path = configure_logging(
            args.outdir, verbose=args.verbose, debug=args.debug
        )
        _log_run_header(argv, args)
        try:
            result = run_pipeline(
                args.inputs,
                args.outdir,
                manifest=args.manifest,
                k=args.k,
                auto_k=args.auto_k,
                min_freq=args.min_freq,
                majority_threshold=args.majority_threshold,
                min_count=args.min_count,
                min_qual=args.min_qual,
                qual_filter=args.qual_filter,
                parsimony=args.parsimony,
                ml=args.ml,
                reference=args.reference,
                map_tree=args.map_tree,
                html=args.html,
                cluster=args.cluster,
                threads=args.threads,
            )
        except (SkaError, ValueError) as exc:
            # Expected, user-actionable failures: a one-line console message,
            # full detail (including the command and stderr) already in the log.
            logger.error("Run failed: %s", exc)
            print(f"error: {exc}", file=sys.stderr)
            print(f"See {log_path} for the full log.", file=sys.stderr)
            return 1
        except Exception:
            # Anything unexpected: capture the traceback in the log so a crash
            # is debuggable from the file alone, then surface where to look.
            logger.exception("Unexpected error during run")
            print("error: unexpected failure; see the log for details.", file=sys.stderr)
            print(f"See {log_path} for the full traceback.", file=sys.stderr)
            return 1
        print(f"Results written to {args.outdir}/")
        print(f"  alignment : {result.alignment.name}")
        print(f"  SNP matrix: {result.snp_matrix.name}")
        print(f"  NJ tree   : {result.nj_tree.name}")
        if result.parsimony_tree:
            print(f"  MP tree   : {result.parsimony_tree.name}")
        if result.ml_tree:
            print(f"  ML tree   : {result.ml_tree.name}")
        if result.ml_report:
            print(f"  ML report : {result.ml_report.name}")
        if result.ref_alignment:
            print(f"  ref align : {result.ref_alignment.name}")
        if result.map_nj_tree:
            print(f"  map tree  : {result.map_nj_tree.name}")
        if result.annotation:
            print(f"  annotation: {result.annotation.name}")
        if result.clusters:
            print(f"  clusters  : {result.clusters.name}")
        if result.html_report:
            print(f"  report    : {result.html_report.name}")
        print(f"  summary   : {result.summary.name}")
        print(f"  log       : {log_path.name}")
        logger.info("Run completed successfully.")
        return 0

    parser.error(f"unknown command {args.command!r}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
