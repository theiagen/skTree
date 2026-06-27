"""Reference-based SNP annotation: which gene a SNP hits, and its codon effect.

This is the one place skTree leaves the reference-free world. The pipeline maps
split k-mers back onto a user-supplied reference with ``ska map -f vcf`` to get
SNP positions in reference coordinates; here we predict genes on that same
reference with **pyrodigal** and label each SNP as intergenic or coding, and for
coding SNPs whether the amino-acid change is synonymous or non-synonymous.

The pure annotation logic (``parse_vcf``, ``annotate_variants``) takes plain
dataclasses so it can be tested without pyrodigal or ska; only ``predict_genes``
imports pyrodigal, and it does so lazily so the module imports cleanly without it.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from Bio.Seq import Seq

from .snps import _read_fasta


def read_reference(path: Path) -> dict[str, str]:
    """Read a reference FASTA into a ``{seqid: uppercase_sequence}`` map."""
    names, seqs = _read_fasta(Path(path).read_text())
    return {name: seq.upper() for name, seq in zip(names, seqs, strict=True)}


@dataclass(frozen=True)
class Gene:
    """A predicted coding region in 1-based inclusive reference coordinates."""

    seqid: str
    start: int
    end: int
    strand: int  # +1 or -1
    gene_id: str


@dataclass(frozen=True)
class Variant:
    """A single-alternate-allele SNP at a 1-based reference position."""

    chrom: str
    pos: int
    ref: str
    alt: str


@dataclass(frozen=True)
class Annotation:
    variant: Variant
    region: str  # "coding" | "intergenic"
    gene_id: str | None
    effect: str | None  # "synonymous" | "nonsynonymous" | None
    ref_aa: str | None
    alt_aa: str | None


def parse_vcf(source: str | Path) -> list[Variant]:
    """Parse ``ska map -f vcf`` output into one ``Variant`` per alternate allele."""
    text = source if isinstance(source, str) and "\n" in source else Path(source).read_text()
    variants: list[Variant] = []
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 5:
            continue
        chrom, pos, _id, ref, alt = fields[:5]
        for allele in alt.split(","):
            if allele in (".", ""):
                continue
            variants.append(Variant(chrom, int(pos), ref, allele))
    return variants


def predict_genes(reference_fasta: Path) -> dict[str, list[Gene]]:
    """Predict genes on each reference contig with pyrodigal (metagenomic mode).

    Metagenomic mode needs no training pass, so it works on a single contig or a
    short reference without the length requirements of self-training.
    """
    import pyrodigal

    finder = pyrodigal.GeneFinder(meta=True)
    by_seq: dict[str, list[Gene]] = {}
    for seqid, seq in read_reference(reference_fasta).items():
        predicted = finder.find_genes(seq.encode())
        genes = [
            Gene(seqid, g.begin, g.end, g.strand, f"{seqid}_{i + 1}")
            for i, g in enumerate(predicted)
        ]
        by_seq[seqid] = genes
    return by_seq


def _containing_gene(genes: list[Gene], pos: int) -> Gene | None:
    for gene in genes:
        if gene.start <= pos <= gene.end:
            return gene
    return None


def _codon_effect(
    ref_seq: str, gene: Gene, pos: int, alt: str
) -> tuple[str | None, str | None, str | None]:
    """Return (effect, ref_aa, alt_aa) for a SNP inside ``gene``.

    The codon is always extracted in forward coordinates first; for a minus-strand
    gene we reverse-complement the codon before translating, so the substitution
    happens in the correct frame.
    """
    if gene.strand >= 0:
        codon_no = (pos - gene.start) // 3
        c0 = (gene.start - 1) + codon_no * 3  # 0-based start of the codon window
    else:
        codon_no = (gene.end - pos) // 3
        c0 = (gene.end - 1) - codon_no * 3 - 2

    if c0 < 0 or c0 + 3 > len(ref_seq):
        return None, None, None

    ref_codon_fwd = ref_seq[c0 : c0 + 3]
    alt_fwd = ref_seq[: pos - 1] + alt + ref_seq[pos:]
    alt_codon_fwd = alt_fwd[c0 : c0 + 3]

    if gene.strand < 0:
        ref_codon = str(Seq(ref_codon_fwd).reverse_complement())
        alt_codon = str(Seq(alt_codon_fwd).reverse_complement())
    else:
        ref_codon, alt_codon = ref_codon_fwd, alt_codon_fwd

    ref_aa = str(Seq(ref_codon).translate())
    alt_aa = str(Seq(alt_codon).translate())
    effect = "synonymous" if ref_aa == alt_aa else "nonsynonymous"
    return effect, ref_aa, alt_aa


def annotate_variants(
    variants: list[Variant],
    genes_by_seq: dict[str, list[Gene]],
    ref_seqs: dict[str, str],
) -> list[Annotation]:
    """Label each variant intergenic/coding and (for coding) its codon effect."""
    annotations: list[Annotation] = []
    for v in variants:
        gene = _containing_gene(genes_by_seq.get(v.chrom, []), v.pos)
        if gene is None:
            annotations.append(Annotation(v, "intergenic", None, None, None, None))
            continue
        effect, ref_aa, alt_aa = _codon_effect(ref_seqs[v.chrom], gene, v.pos, v.alt)
        annotations.append(Annotation(v, "coding", gene.gene_id, effect, ref_aa, alt_aa))
    return annotations


def write_annotation_tsv(path: Path, annotations: list[Annotation]) -> None:
    header = "chrom\tpos\tref\talt\tregion\tgene\teffect\tref_aa\talt_aa\n"
    lines = [header]
    for a in annotations:
        v = a.variant
        lines.append(
            "\t".join(
                str(x) if x is not None else ""
                for x in (
                    v.chrom, v.pos, v.ref, v.alt,
                    a.region, a.gene_id, a.effect, a.ref_aa, a.alt_aa,
                )
            )
            + "\n"
        )
    path.write_text("".join(lines))
