"""Tests for reference-based SNP annotation (gene membership + codon effect).

The pure annotation logic is tested with hand-built ``Gene``/``Variant`` objects
so it needs neither pyrodigal nor ska. ``predict_genes`` is exercised separately
and skips when pyrodigal is absent.
"""

from __future__ import annotations

import pytest

from sktree.annotate import (
    Gene,
    Variant,
    annotate_variants,
    parse_vcf,
)

# ATG AAA TTT GGG TAA  ->  M K F G *   (a clean forward ORF)
PLUS_REF = "ATGAAATTTGGGTAA"
# reverse complement of PLUS_REF: same protein, encoded on the minus strand
MINUS_REF = "TTACCCAAATTTCAT"


def test_parse_vcf_basic():
    vcf = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n"
        "chrom1\t501\t.\tG\tA\t.\t.\t.\tGT\t0\t1\n"
        "chrom1\t1501\t.\tC\tA,T\t.\t.\t.\tGT\t0\t1\n"
    )
    variants = parse_vcf(vcf)
    assert variants[0] == Variant("chrom1", 501, "G", "A")
    # multiallelic ALT splits into one Variant per alternate allele
    assert Variant("chrom1", 1501, "C", "A") in variants
    assert Variant("chrom1", 1501, "C", "T") in variants
    assert len(variants) == 3


def test_intergenic_when_no_gene_covers_position():
    genes = {"chrom1": [Gene("chrom1", 1, 15, +1, "gene_1")]}
    refs = {"chrom1": PLUS_REF}
    ann = annotate_variants([Variant("chrom1", 100, "A", "C")], genes, refs)[0]
    assert ann.region == "intergenic"
    assert ann.gene_id is None
    assert ann.effect is None


def test_coding_nonsynonymous_plus_strand():
    genes = {"chrom1": [Gene("chrom1", 1, 15, +1, "gene_1")]}
    refs = {"chrom1": PLUS_REF}
    # pos 4 = first base of codon 2 (AAA=K); A->C gives CAA=Q
    ann = annotate_variants([Variant("chrom1", 4, "A", "C")], genes, refs)[0]
    assert ann.region == "coding"
    assert ann.gene_id == "gene_1"
    assert (ann.ref_aa, ann.alt_aa) == ("K", "Q")
    assert ann.effect == "nonsynonymous"


def test_coding_synonymous_plus_strand():
    genes = {"chrom1": [Gene("chrom1", 1, 15, +1, "gene_1")]}
    refs = {"chrom1": PLUS_REF}
    # pos 6 = third base of codon 2 (AAA=K); A->G gives AAG=K (silent)
    ann = annotate_variants([Variant("chrom1", 6, "A", "G")], genes, refs)[0]
    assert ann.effect == "synonymous"
    assert ann.ref_aa == ann.alt_aa == "K"


def test_coding_minus_strand_translates_reverse_complement():
    genes = {"chrom1": [Gene("chrom1", 1, 15, -1, "gene_1")]}
    refs = {"chrom1": MINUS_REF}
    # forward pos 14 is the middle base of the start codon (ATG=M after revcomp);
    # A->G turns the codon into ACG=T -> nonsynonymous
    ann = annotate_variants([Variant("chrom1", 14, "A", "G")], genes, refs)[0]
    assert ann.region == "coding"
    assert ann.ref_aa == "M"
    assert ann.alt_aa == "T"
    assert ann.effect == "nonsynonymous"


def test_predict_genes_finds_a_known_orf():
    pyrodigal = pytest.importorskip("pyrodigal")  # noqa: F841
    from sktree.annotate import predict_genes

    # a long ORF padded with flanks so metagenomic gene-finding can call it
    orf = "ATG" + "AAACGTGGCATTGAC" * 12 + "TAA"
    seq = "TT" * 40 + orf + "TT" * 40
    genes = predict_genes_from_string(predict_genes, seq)
    assert any(g.strand in (+1, -1) for g in genes)
    assert any(g.end - g.start + 1 >= 60 for g in genes)


def predict_genes_from_string(predict_genes, seq):
    """Helper: write a one-record FASTA to a temp path and predict on it."""
    import tempfile
    from pathlib import Path

    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "ref.fasta"
        p.write_text(f">chrom1\n{seq}\n")
        by_seq = predict_genes(p)
    return by_seq.get("chrom1", [])
