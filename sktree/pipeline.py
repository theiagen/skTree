"""End-to-end orchestration: FASTA/FASTQ inputs -> SNPs -> tree -> reports.

This is the body of ``sktree run``. It drives the SKA2 engine for the heavy
k-mer work, then applies the Python analysis layer (SNP partitioning, NJ
tree, reports).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

from .annotate import (
    annotate_variants,
    parse_vcf,
    predict_genes,
    read_reference,
    write_annotation_tsv,
)
from .engine.cluster import FastbapsError, FastbapsNotAvailable, FastbapsRunner
from .engine.ml import MlError, MlNotAvailable, MlTreeBuilder
from .engine.ska import SkaError, SkaRunner
from .html_report import build_report_data, write_html_report
from .kselect import select_k
from .parsimony import parsimony_tree
from .phylo import neighbor_joining, to_newick
from .report import build_stats, format_summary, write_snp_matrix_tsv
from .snps import parse_alignment

logger = logging.getLogger("sktree")

# SKA's exact-flank matching loses recall beyond ~1% sequence divergence.
DIVERGENCE_WARN_THRESHOLD = 0.01

# Pure-Python parsimony search is O(n^3 * sites); warn past this many samples.
PARSIMONY_SLOW_ABOVE = 30


@dataclass
class RunResult:
    skf: Path
    alignment: Path
    snp_matrix: Path
    core_alignment: Path
    nj_tree: Path
    distance: Path
    summary: Path
    max_divergence: float | None
    parsimony_tree: Path | None = None
    ml_tree: Path | None = None
    ml_report: Path | None = None
    ref_alignment: Path | None = None
    map_nj_tree: Path | None = None
    map_parsimony_tree: Path | None = None
    map_ml_tree: Path | None = None
    map_ml_report: Path | None = None
    annotation: Path | None = None
    clusters: Path | None = None
    html_report: Path | None = None


def _output_listing(files: list[tuple[Path | None, str]]) -> list[dict]:
    """Name + description for every result file that was actually written.

    Entries whose path is ``None`` (optional steps that did not run) are
    dropped, so the report only ever lists files the user can find on disk.
    """
    return [
        {"name": path.name, "description": description}
        for path, description in files
        if path is not None
    ]


def _max_divergence(distance_tsv: Path) -> float | None:
    """Largest pairwise sequence divergence across sample pairs, if available.

    ``ska distance`` columns are: Sample1, Sample2, Distance (SNP count),
    Mismatches (k-mer proportion), Match count, Mismatch count. A single SNP
    invalidates ~k split-k-mers, so the k-mer proportion overstates divergence
    by ~k-fold. The faithful per-base substitution rate is the SNP count over
    the number of comparable split-k-mer positions (match + mismatch).
    """
    divergences: list[float] = []
    for line in distance_tsv.read_text().splitlines():
        if not line or line.startswith("Sample1"):
            continue
        fields = line.split("\t")
        if len(fields) >= 6:
            try:
                snps = float(fields[2])
                compared = float(fields[4]) + float(fields[5])
            except ValueError:
                continue
            if compared > 0:
                divergences.append(snps / compared)
    return max(divergences) if divergences else None


@dataclass
class MlOutputs:
    """An ML run's tree plus the engine's auxiliary files.

    IQ-TREE writes a whole family of files next to the tree -- a human-readable
    ``.iqtree`` report (chosen model, log-likelihood, tree), a ``.log`` run log,
    and checkpoint/distance side files. We surface the report and log explicitly
    and keep the rest grouped so the run's full ML output is collected, not just
    the bare Newick tree.
    """

    tree: Path
    report: Path | None
    log: Path | None
    extra: list[Path]


def _collect_ml_outputs(prefix: Path, treefile: Path) -> MlOutputs:
    """Gather every file the ML engine wrote under ``prefix`` beside the tree."""
    produced = sorted(p for p in prefix.parent.glob(prefix.name + ".*") if p != treefile)
    report = next((p for p in produced if p.suffix == ".iqtree"), None)
    log = next((p for p in produced if p.suffix == ".log"), None)
    extra = [p for p in produced if p not in (report, log)]
    return MlOutputs(tree=treefile, report=report, log=log, extra=extra)


def _build_ml_tree(
    alignment: Path, outdir: Path, threads: int | None, *, prefix_name: str = "tree_ml"
) -> MlOutputs | None:
    """Run ML inference on the core SNP alignment; skip gracefully if no engine.

    Uses the engine's own default parameters (IQ-TREE's ModelFinder picks the
    substitution model). Returns the tree together with the engine's auxiliary
    files, or ``None`` when no engine is available or the run fails. ``prefix_name``
    names the IQ-TREE output family so reference-anchored runs can sit beside the
    reference-free tree (e.g. ``tree_ref_ml`` vs ``tree_ml``).
    """
    try:
        builder = MlTreeBuilder()
    except MlNotAvailable as exc:
        logger.warning("Skipping ML tree: %s", exc)
        return None
    logger.info("Inferring maximum-likelihood tree with %s (engine defaults)", builder.engine)
    prefix = outdir / prefix_name
    try:
        treefile = builder.build(alignment, prefix, threads=threads)
    except MlError as exc:
        logger.warning("ML inference failed: %s", exc)
        return None
    outputs = _collect_ml_outputs(prefix, treefile)
    n_collected = 1 + len(outputs.extra) + (outputs.report is not None) + (outputs.log is not None)
    logger.info("Collected %d ML output file(s) under %s.*", n_collected, prefix.name)
    return outputs


def _build_clusters(
    core_alignment: Path,
    outdir: Path,
    *,
    threads: int | None,
    write_csv: bool,
) -> tuple[Path | None, dict[str, int] | None]:
    """Run fastbaps on the core alignment; skip gracefully if unavailable."""
    try:
        runner = FastbapsRunner()
    except FastbapsNotAvailable as exc:
        logger.warning("Skipping clustering: %s", exc)
        return None, None
    out = outdir / "clusters.csv"
    try:
        clusters = runner.cluster(core_alignment, out, threads=threads)
    except FastbapsError as exc:
        logger.warning("fastbaps failed: %s", exc)
        return None, None
    logger.info("fastbaps assigned %d clusters", len(set(clusters.values())))
    return (out if write_csv else None), clusters


def _annotation_counts(annotation_tsv: Path) -> dict[str, int]:
    coding = nonsyn = intergenic = 0
    lines = annotation_tsv.read_text().splitlines()
    header = lines[0].split("\t") if lines else []
    try:
        region_i = header.index("region")
        effect_i = header.index("effect")
    except ValueError:
        return {"coding": 0, "nonsynonymous": 0, "intergenic": 0}
    for line in lines[1:]:
        fields = line.split("\t")
        if len(fields) <= max(region_i, effect_i):
            continue
        if fields[region_i] == "coding":
            coding += 1
            if fields[effect_i] == "nonsynonymous":
                nonsyn += 1
        else:
            intergenic += 1
    return {"coding": coding, "nonsynonymous": nonsyn, "intergenic": intergenic}


def _annotate_snps(
    runner: SkaRunner,
    skf: Path,
    reference: Path,
    outdir: Path,
    threads: int | None,
) -> Path:
    """Map SNPs onto a reference, predict genes, and write the annotation TSV."""
    logger.info("Annotating SNPs against reference %s", reference.name)
    vcf = outdir / "reference_snps.vcf"
    runner.map(reference, skf, output=vcf, out_format="vcf", threads=threads)
    variants = parse_vcf(vcf)
    genes_by_seq = predict_genes(reference)
    ref_seqs = read_reference(reference)
    annotations = annotate_variants(variants, genes_by_seq, ref_seqs)
    out = outdir / "snp_annotation.tsv"
    write_annotation_tsv(out, annotations)
    coding = sum(1 for a in annotations if a.region == "coding")
    nonsyn = sum(1 for a in annotations if a.effect == "nonsynonymous")
    logger.info(
        "Annotated %d variants: %d coding (%d non-synonymous), %d intergenic",
        len(annotations), coding, nonsyn, len(annotations) - coding,
    )
    return out


def _reference_align(
    runner: SkaRunner,
    skf: Path,
    reference: Path,
    outdir: Path,
    threads: int | None,
) -> Path | None:
    """Write the reference-coordinate pseudo-alignment via ``ska map -f aln``.

    Unlike the reference-free ``align`` output, this is full reference length
    with one row per sample at reference coordinates -- the basis for
    reference-anchored trees. Returns ``None`` (and logs) if mapping fails, so
    the rest of the run still completes.
    """
    out = outdir / "ref_aligned.fasta"
    logger.info("Mapping split k-mers onto reference %s (pseudo-alignment)", reference.name)
    try:
        runner.map(reference, skf, output=out, out_format="aln", threads=threads)
    except SkaError as exc:
        logger.warning("Reference mapping failed: %s", exc)
        return None
    return out


@dataclass
class MapOutputs:
    """Reference-anchored trees built from the ``ska map`` pseudo-alignment."""

    nj_tree: Path
    parsimony_tree: Path | None
    ml_tree: Path | None
    ml_report: Path | None


def _build_map_trees(
    ref_alignment: Path,
    outdir: Path,
    *,
    parsimony: bool,
    ml: bool,
    threads: int | None,
) -> MapOutputs:
    """Build reference-anchored trees from the pseudo-alignment.

    NJ and parsimony run on the variable subset (the reference-coordinate
    equivalent of SKA's reference-free SNP alignment); ML runs on the core-SNP
    subset file, mirroring the reference-free path. ``classify`` scans the full
    reference-length matrix once -- the performance watch-point on large
    references (see FOR-DEVELOPERS.md).
    """
    full = parse_alignment(ref_alignment)
    cls = full.classify()
    variable = full.subset(cls.variable)

    nj = outdir / "tree_ref_nj.nwk"
    nj.write_text(to_newick(neighbor_joining(variable)))

    parsimony_path: Path | None = None
    if parsimony:
        logger.info("Inferring reference-anchored maximum-parsimony tree")
        tree, _ = parsimony_tree(variable)
        parsimony_path = outdir / "tree_ref_parsimony.nwk"
        parsimony_path.write_text(to_newick(tree))

    ml_tree: Path | None = None
    ml_report: Path | None = None
    if ml:
        core = full.subset(cls.core)
        core_fasta = outdir / "ref_core_snps.fasta"
        core_fasta.write_text(core.to_fasta())
        ml_out = _build_ml_tree(core_fasta, outdir, threads, prefix_name="tree_ref_ml")
        if ml_out is not None:
            ml_tree = ml_out.tree
            ml_report = ml_out.report
    return MapOutputs(
        nj_tree=nj, parsimony_tree=parsimony_path, ml_tree=ml_tree, ml_report=ml_report
    )


def run_pipeline(
    inputs: list[Path],
    outdir: Path,
    *,
    k: int = 31,
    auto_k: bool = False,
    min_freq: float = 0.9,
    majority_threshold: float = 0.5,
    parsimony: bool = False,
    ml: bool = False,
    reference: Path | None = None,
    map_tree: bool = False,
    html: bool = False,
    cluster: bool = False,
    threads: int | None = None,
    runner: SkaRunner | None = None,
) -> RunResult:
    if len(inputs) < 2:
        raise ValueError("need at least two input genomes to call SNPs")
    outdir.mkdir(parents=True, exist_ok=True)
    runner = runner or SkaRunner()

    if auto_k:
        selection = select_k(inputs)
        k = selection.k
        logger.info(
            "Auto-selected k=%d (uniqueness %.4f on %s)",
            k,
            selection.uniqueness,
            selection.representative.name,
        )

    logger.info("Building split-k-mer file (k=%d) from %d inputs", k, len(inputs))
    skf = runner.build(inputs, out_prefix=outdir / "combined", k=k, threads=threads)

    logger.info("Writing reference-free SNP alignment (min-freq=%s)", min_freq)
    alignment = outdir / "alignment.fasta"
    runner.align(skf, output=alignment, min_freq=min_freq, threads=threads)

    distance = outdir / "distances.tsv"
    runner.distance(skf, output=distance, threads=threads)
    max_div = _max_divergence(distance)
    if max_div is not None and max_div > DIVERGENCE_WARN_THRESHOLD:
        logger.warning(
            "Max pairwise divergence %.2f%% exceeds ~%.0f%%; SKA recall degrades "
            "on divergent samples. Treat results with caution.",
            max_div * 100,
            DIVERGENCE_WARN_THRESHOLD * 100,
        )

    m = parse_alignment(alignment)
    cls = m.classify(majority_threshold=majority_threshold)

    snp_matrix = outdir / "snp_matrix.tsv"
    write_snp_matrix_tsv(snp_matrix, m, cls)

    core = m.subset(cls.core)
    core_alignment = outdir / "core_snps.fasta"
    core_alignment.write_text(core.to_fasta())

    logger.info("Inferring neighbor-joining tree")
    nj_tree = outdir / "tree_nj.nwk"
    nj_tree.write_text(to_newick(neighbor_joining(m)))

    parsimony_path: Path | None = None
    if parsimony:
        if m.n_samples > PARSIMONY_SLOW_ABOVE:
            logger.warning(
                "Parsimony search on %d samples may be slow (pure-Python, ~O(n^3)).",
                m.n_samples,
            )
        logger.info("Inferring maximum-parsimony tree")
        tree, score = parsimony_tree(m)
        parsimony_path = outdir / "tree_parsimony.nwk"
        parsimony_path.write_text(to_newick(tree))
        logger.info("Parsimony tree score: %d steps", score)

    ml_path: Path | None = None
    ml_outputs: MlOutputs | None = None
    if ml:
        ml_outputs = _build_ml_tree(core_alignment, outdir, threads)
        if ml_outputs is not None:
            ml_path = ml_outputs.tree

    ref_alignment_path: Path | None = None
    annotation_path: Path | None = None
    if reference is not None:
        ref_alignment_path = _reference_align(runner, skf, reference, outdir, threads)
        annotation_path = _annotate_snps(runner, skf, reference, outdir, threads)

    map_outputs: MapOutputs | None = None
    if map_tree and ref_alignment_path is not None:
        logger.info("Building reference-anchored tree(s) from the pseudo-alignment")
        map_outputs = _build_map_trees(
            ref_alignment_path, outdir, parsimony=parsimony, ml=ml, threads=threads
        )

    clusters_path: Path | None = None
    clusters_map: dict[str, int] | None = None
    if cluster or html:
        clusters_path, clusters_map = _build_clusters(
            core_alignment, outdir, threads=threads, write_csv=cluster
        )

    stats = build_stats(m, cls, k=k, min_freq=min_freq, majority_threshold=majority_threshold)
    summary = outdir / "summary.txt"
    summary.write_text(format_summary(stats))
    logger.info(
        "Done: %d all / %d core / %d majority SNPs",
        stats.n_all_snps,
        stats.n_core_snps,
        stats.n_majority_snps,
    )

    html_path: Path | None = None
    if html:
        trees: dict[str, str] = {"nj": nj_tree.read_text().strip()}
        if parsimony_path is not None:
            trees["parsimony"] = parsimony_path.read_text().strip()
        if ml_path is not None:
            trees["ml"] = ml_path.read_text().strip()
        if map_outputs is not None:
            trees["ref_nj"] = map_outputs.nj_tree.read_text().strip()
            if map_outputs.parsimony_tree is not None:
                trees["ref_parsimony"] = map_outputs.parsimony_tree.read_text().strip()
            if map_outputs.ml_tree is not None:
                trees["ref_ml"] = map_outputs.ml_tree.read_text().strip()
        annotation_summary = None
        if annotation_path is not None:
            annotation_summary = _annotation_counts(annotation_path)
        outputs = _output_listing([
            (alignment, "Reference-free SNP alignment: every variable site, one row per sample"),
            (core_alignment, "Core SNP alignment: sites shared by all samples (tree input)"),
            (snp_matrix, "Loci-by-sample SNP matrix, each locus tagged core / majority / other"),
            (distance, "Pairwise SNP-distance matrix between every sample pair"),
            (nj_tree, "Neighbor-joining tree (Newick)"),
            (parsimony_path, "Maximum-parsimony tree (Newick)"),
            (ml_path, "Maximum-likelihood tree from IQ-TREE (Newick)"),
            (ml_outputs.report if ml_outputs else None,
             "IQ-TREE report: selected substitution model, log-likelihood, and tree"),
            (ml_outputs.log if ml_outputs else None, "IQ-TREE run log"),
            (ref_alignment_path,
             "Reference-anchored pseudo-alignment (ska map): one row per sample at "
             "reference coordinates"),
            (map_outputs.nj_tree if map_outputs else None,
             "Reference-anchored neighbor-joining tree (Newick)"),
            (map_outputs.parsimony_tree if map_outputs else None,
             "Reference-anchored maximum-parsimony tree (Newick)"),
            (map_outputs.ml_tree if map_outputs else None,
             "Reference-anchored maximum-likelihood tree from IQ-TREE (Newick)"),
            (map_outputs.ml_report if map_outputs else None,
             "Reference-anchored IQ-TREE report: model, log-likelihood, and tree"),
            (clusters_path, "fastbaps population-structure cluster assignments (CSV)"),
            (annotation_path, "Per-SNP annotation against the reference: gene and effect (TSV)"),
            (summary, "Plain-text run summary"),
        ])
        report_data = build_report_data(
            stats=stats,
            trees=trees,
            core_matrix=core,
            all_matrix=m.subset(cls.variable),
            clusters=clusters_map,
            annotation_summary=annotation_summary,
            max_divergence=max_div,
            divergence_warn_threshold=DIVERGENCE_WARN_THRESHOLD,
            outputs=outputs,
        )
        html_path = outdir / "report.html"
        write_html_report(html_path, report_data)
        logger.info("Wrote HTML report to %s", html_path.name)

    return RunResult(
        skf=skf,
        alignment=alignment,
        snp_matrix=snp_matrix,
        core_alignment=core_alignment,
        nj_tree=nj_tree,
        distance=distance,
        summary=summary,
        max_divergence=max_div,
        parsimony_tree=parsimony_path,
        ml_tree=ml_path,
        ml_report=ml_outputs.report if ml_outputs else None,
        ref_alignment=ref_alignment_path,
        map_nj_tree=map_outputs.nj_tree if map_outputs else None,
        map_parsimony_tree=map_outputs.parsimony_tree if map_outputs else None,
        map_ml_tree=map_outputs.ml_tree if map_outputs else None,
        map_ml_report=map_outputs.ml_report if map_outputs else None,
        annotation=annotation_path,
        clusters=clusters_path,
        html_report=html_path,
    )
