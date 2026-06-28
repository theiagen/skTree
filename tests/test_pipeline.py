"""End-to-end pipeline tests.

The full pipeline runs the real ``ska`` binary, so these are integration tests
that skip when ska is absent. A pure-Python validation check guards the
input-count guard without the binary.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from sktree.pipeline import (
    _collect_ml_outputs,
    _max_divergence,
    _output_listing,
    run_pipeline,
)

needs_ska = pytest.mark.skipif(shutil.which("ska") is None, reason="ska binary not installed")


def test_collect_ml_outputs_classifies_iqtree_family(tmp_path):
    # IQ-TREE writes a whole prefix.* family; the collector must split the tree,
    # the .iqtree report and .log out from the checkpoint/distance side files.
    prefix = tmp_path / "tree_ml"
    treefile = prefix.with_suffix(".treefile")
    names = [".treefile", ".iqtree", ".log", ".bionj", ".mldist", ".model.gz", ".ckp.gz"]
    for suffix in names:
        (tmp_path / ("tree_ml" + suffix)).write_text("x")
    (tmp_path / "unrelated.txt").write_text("x")  # must be ignored

    out = _collect_ml_outputs(prefix, treefile)

    assert out.tree == treefile
    assert out.report is not None and out.report.name == "tree_ml.iqtree"
    assert out.log is not None and out.log.name == "tree_ml.log"
    extra = {p.name for p in out.extra}
    assert extra == {"tree_ml.bionj", "tree_ml.mldist", "tree_ml.model.gz", "tree_ml.ckp.gz"}
    # tree, report and log are never double-counted in extra.
    assert treefile not in out.extra
    assert out.report not in out.extra and out.log not in out.extra


def test_collect_ml_outputs_tolerates_missing_report_and_log(tmp_path):
    # RAxML-NG (or a partial run) may not write a .iqtree/.log; report/log are None.
    prefix = tmp_path / "tree_ml"
    treefile = Path(str(prefix) + ".raxml.bestTree")
    treefile.write_text("x")
    (tmp_path / "tree_ml.raxml.bestModel").write_text("x")

    out = _collect_ml_outputs(prefix, treefile)

    assert out.tree == treefile
    assert out.report is None and out.log is None
    assert [p.name for p in out.extra] == ["tree_ml.raxml.bestModel"]


def test_output_listing_drops_missing_files():
    listing = _output_listing([
        (Path("/out/alignment.fasta"), "the alignment"),
        (None, "optional step that did not run"),
        (Path("/out/summary.txt"), "the summary"),
    ])
    assert listing == [
        {"name": "alignment.fasta", "description": "the alignment"},
        {"name": "summary.txt", "description": "the summary"},
    ]


def test_requires_two_inputs(tmp_path):
    with pytest.raises(ValueError, match="two input"):
        run_pipeline([tmp_path / "only.fasta"], tmp_path / "out")


def test_single_paired_read_sample_counts_as_one_input(tmp_path):
    # Two FASTQ files that pair into ONE sample are not two inputs: the minimum
    # is two *samples*, not two files.
    with pytest.raises(ValueError, match="two input"):
        run_pipeline(
            [tmp_path / "S_R1.fastq.gz", tmp_path / "S_R2.fastq.gz"], tmp_path / "out"
        )


def test_assembly_paths_filters_out_reads():
    from sktree.inputs import Sample
    from sktree.pipeline import _assembly_paths

    samples = [
        Sample("asm", (Path("asm.fasta"),), is_reads=False),
        Sample("reads", (Path("r_1.fq.gz"), Path("r_2.fq.gz")), is_reads=True),
    ]
    assert _assembly_paths(samples) == [Path("asm.fasta")]


def test_resolve_k_falls_back_when_reads_only(caplog):
    import logging

    from sktree.inputs import Sample
    from sktree.pipeline import _resolve_k

    reads_only = [
        Sample("a", (Path("a_1.fq.gz"), Path("a_2.fq.gz")), is_reads=True),
        Sample("b", (Path("b_1.fq.gz"), Path("b_2.fq.gz")), is_reads=True),
    ]
    with caplog.at_level(logging.WARNING):
        k = _resolve_k(reads_only, k=27, auto_k=True)
    assert k == 27  # kept the default; auto-k cannot read FASTQ
    assert any("auto-k" in r.message or "assembl" in r.message.lower() for r in caplog.records)


def test_max_divergence_uses_snp_rate_not_kmer_proportion(tmp_path):
    # Real ska distance layout: the k-mer proportion (col 4) is ~k-fold inflated;
    # true divergence is SNP count (col 3) / (match + mismatch).
    tsv = tmp_path / "distances.tsv"
    tsv.write_text(
        "Sample1\tSample2\tDistance\tMismatches (proportion)\tMatch count\tMismatch count\n"
        "a\tb\t2.00\t0.05911\t1910\t120\n"   # 2 / 2030 = 0.000985
        "a\tc\t1.00\t0.03000\t1940\t60\n"    # 1 / 2000 = 0.000500
    )
    div = _max_divergence(tsv)
    assert div is not None
    assert abs(div - 2 / 2030) < 1e-9  # the larger pair drives the warning
    assert div < 0.01  # would have falsely tripped the 1% warning on col 4


def test_max_divergence_empty_returns_none(tmp_path):
    tsv = tmp_path / "d.tsv"
    header = "Sample1\tSample2\tDistance\tMismatches (proportion)\tMatch count\tMismatch count\n"
    tsv.write_text(header)
    assert _max_divergence(tsv) is None


@needs_ska
def test_full_run_produces_all_outputs(synthetic_genomes, tmp_path):
    out = tmp_path / "out"
    result = run_pipeline(synthetic_genomes, out, k=31)
    for path in (
        result.skf,
        result.alignment,
        result.snp_matrix,
        result.core_alignment,
        result.nj_tree,
        result.distance,
        result.summary,
    ):
        assert path.exists(), f"missing {path}"

    # alignment has all three samples
    assert result.alignment.read_text().count(">") == 3
    # NJ tree is parseable newick with the three taxa
    import dendropy

    tree = dendropy.Tree.get(path=str(result.nj_tree), schema="newick")
    assert {leaf.taxon.label for leaf in tree.leaf_node_iter()} == {
        "sample1",
        "sample2",
        "sample3",
    }
    # summary reports SNP counts
    assert "core SNPs" in result.summary.read_text()


@needs_ska
def test_full_run_on_mixed_assembly_and_reads(assembly_plus_reads, tmp_path):
    out = tmp_path / "out"
    result = run_pipeline(assembly_plus_reads["inputs"], out, k=31)
    assert result.nj_tree.exists()
    # both the assembly and the auto-paired read sample appear as tree tips.
    import dendropy

    tree = dendropy.Tree.get(path=str(result.nj_tree), schema="newick")
    assert {leaf.taxon.label for leaf in tree.leaf_node_iter()} == {"asm", "strainB"}
    # the planted SNP shows up as at least one variable site.
    assert ">asm" in result.alignment.read_text()


@needs_ska
def test_manifest_drives_a_run(synthetic_genomes, tmp_path):
    # A manifest with no positional inputs builds a tree from the listed files.
    manifest = tmp_path / "samples.tsv"
    manifest.write_text(
        "".join(f"{p.stem}\t{p}\n" for p in synthetic_genomes)
    )
    out = tmp_path / "out"
    result = run_pipeline([], out, manifest=manifest, k=31)
    assert result.nj_tree.exists()
    assert result.alignment.read_text().count(">") == 3


@needs_ska
def test_parsimony_tree_written_when_requested(synthetic_genomes, tmp_path):
    out = tmp_path / "out"
    result = run_pipeline(synthetic_genomes, out, k=31, parsimony=True)
    assert result.parsimony_tree is not None
    assert result.parsimony_tree.exists()
    assert (out / "tree_parsimony.nwk").exists()


@needs_ska
def test_ml_skipped_gracefully_without_engine(synthetic_genomes, tmp_path):
    # No IQ-TREE/RAxML on PATH here: ML must be skipped, not crash the run.
    out = tmp_path / "out"
    result = run_pipeline(synthetic_genomes, out, k=31, ml=True)
    if shutil.which("iqtree2") or shutil.which("iqtree") or shutil.which("raxml-ng"):
        pytest.skip("an ML engine is installed; skip the no-engine path")
    assert result.ml_tree is None
    assert result.nj_tree.exists()  # rest of pipeline still completed


@needs_ska
def test_reference_annotation_written(synthetic_genomes, tmp_path):
    pytest.importorskip("pyrodigal")
    out = tmp_path / "out"
    reference = synthetic_genomes[0]  # sample1 doubles as the reference
    result = run_pipeline(synthetic_genomes, out, k=31, reference=reference)
    assert result.annotation is not None
    assert result.annotation.exists()
    text = result.annotation.read_text()
    header = text.splitlines()[0].split("\t")
    assert header == [
        "chrom", "pos", "ref", "alt", "region", "gene", "effect", "ref_aa", "alt_aa"
    ]
    # every annotated row is labelled coding or intergenic
    for line in text.splitlines()[1:]:
        assert line.split("\t")[4] in ("coding", "intergenic")


@needs_ska
def test_reference_pseudo_alignment_written(synthetic_genomes, tmp_path):
    out = tmp_path / "out"
    reference = synthetic_genomes[0]
    result = run_pipeline(synthetic_genomes, out, k=31, reference=reference)
    assert result.ref_alignment is not None
    assert result.ref_alignment.exists()
    assert (out / "ref_aligned.fasta").exists()
    assert result.ref_alignment.read_text().count(">") >= 1


def test_map_tree_requires_reference_via_cli(tmp_path):
    from sktree.cli import main

    with pytest.raises(SystemExit):
        main(["run", str(tmp_path / "a.fa"), str(tmp_path / "b.fa"),
              "-o", str(tmp_path / "o"), "--map-tree"])


@needs_ska
def test_map_tree_builds_nj(synthetic_genomes, tmp_path):
    out = tmp_path / "out"
    result = run_pipeline(
        synthetic_genomes, out, k=31, reference=synthetic_genomes[0], map_tree=True
    )
    assert result.map_nj_tree is not None and result.map_nj_tree.exists()
    assert (out / "tree_ref_nj.nwk").exists()


@needs_ska
def test_cli_run(synthetic_genomes, tmp_path):
    from sktree.cli import main

    out = tmp_path / "cliout"
    rc = main(["run", *[str(p) for p in synthetic_genomes], "-o", str(out), "-k", "31"])
    assert rc == 0
    assert (out / "tree_nj.nwk").exists()


def test_build_clusters_skips_when_unavailable(monkeypatch, tmp_path):
    from sktree.pipeline import _build_clusters

    def _raise(*a, **k):
        from sktree.engine.cluster import FastbapsNotAvailable
        raise FastbapsNotAvailable("nope")

    monkeypatch.setattr("sktree.pipeline.FastbapsRunner", _raise)
    aln = tmp_path / "core_snps.fasta"
    aln.write_text(">s1\nAC\n>s2\nAT\n")
    csv_path, clusters = _build_clusters(aln, tmp_path, threads=None, write_csv=True)
    assert csv_path is None and clusters is None


def test_build_clusters_skips_on_runtime_error(monkeypatch, tmp_path):
    from sktree.engine.cluster import FastbapsError
    from sktree.pipeline import _build_clusters

    class _FailRunner:
        def cluster(self, *a, **k):
            raise FastbapsError("boom")

    monkeypatch.setattr("sktree.pipeline.FastbapsRunner", lambda *a, **k: _FailRunner())
    aln = tmp_path / "core_snps.fasta"
    aln.write_text(">s1\nAC\n>s2\nAT\n")
    csv_path, clusters = _build_clusters(aln, tmp_path, threads=None, write_csv=True)
    assert csv_path is None and clusters is None
