"""Subprocess wrapper around the SKA2 (``ska``) command-line tool.

skTree offloads all k-mer-heavy work to the Rust ``ska`` binary and builds the
analysis layer in Python on top. This module is the single place that knows
how to invoke ``ska`` and turn its results into Python values.

The minimal pipeline is ``build`` (one .skf holding all samples) then ``align``
(a reference-free SNP alignment). ``merge`` exists to combine separately-built
.skf files; ``distance`` yields pairwise SNP distances.
"""

from __future__ import annotations

import shutil
import subprocess
from collections.abc import Sequence
from pathlib import Path

PathLike = str | Path

_COMPRESSION_SUFFIXES = {".gz", ".bz2", ".xz", ".zst", ".zip"}
_SEQUENCE_SUFFIXES = {".fasta", ".fa", ".fna", ".ffn", ".faa", ".fastq", ".fq"}


def clean_sample_name(path: PathLike) -> str:
    """Derive a tidy sample label from an input file path.

    Strips the directory and any trailing compression and FASTA/FASTQ
    extensions, so ``/data/SAMD00052601.fa.gz`` becomes ``SAMD00052601``.
    ``ska`` otherwise names samples by the verbatim path it is handed, which
    makes for unreadable tree tips and alignment rows.
    """
    name = Path(path).name
    base = name
    for suffix in reversed(Path(name).suffixes):
        if suffix.lower() in _COMPRESSION_SUFFIXES or suffix.lower() in _SEQUENCE_SUFFIXES:
            base = base[: -len(suffix)]
        else:
            break
    return base or name


def _build_file_list(seq_files: Sequence[PathLike]) -> list[tuple[str, str]]:
    """Map each input file to a unique clean sample name (``name``, ``path``)."""
    seen: dict[str, int] = {}
    rows: list[tuple[str, str]] = []
    for f in seq_files:
        name = clean_sample_name(f)
        if name in seen:
            seen[name] += 1
            name = f"{name}_{seen[name]}"
        else:
            seen[name] = 1
        rows.append((name, str(f)))
    return rows


class SkaError(RuntimeError):
    """Raised when an invocation of ``ska`` exits non-zero."""


class SkaNotFoundError(SkaError):
    """Raised when the ``ska`` binary cannot be located on PATH."""


class SkaRunner:
    """Thin, testable wrapper over the ``ska`` CLI.

    Parameters
    ----------
    binary:
        Name or path of the ska executable. Resolved via :func:`shutil.which`;
        a :class:`SkaNotFoundError` is raised if it cannot be found.
    """

    def __init__(self, binary: str = "ska") -> None:
        resolved = shutil.which(binary)
        if resolved is None:
            raise SkaNotFoundError(
                f"Could not find the '{binary}' executable on PATH. "
                "Install SKA2 with `cargo install ska` or "
                "`conda install -c bioconda ska2`."
            )
        self.binary = resolved

    # -- internal -------------------------------------------------------------

    def _run(self, args: Sequence[str], capture: bool = True) -> str:
        """Run ``ska <args>``; return stdout. Raise :class:`SkaError` on failure."""
        argv = [self.binary, *args]
        proc = subprocess.run(
            argv,
            capture_output=capture,
            text=True,
        )
        if proc.returncode != 0:
            raise SkaError(
                f"ska exited with code {proc.returncode}: {proc.stderr or proc.stdout}"
            )
        return proc.stdout

    # -- subcommands ----------------------------------------------------------

    def version(self) -> str:
        return self._run(["--version"]).strip()

    def build(
        self,
        seq_files: Sequence[PathLike],
        out_prefix: PathLike,
        k: int = 31,
        *,
        single_strand: bool = False,
        min_count: int | None = None,
        min_qual: int | None = None,
        threads: int | None = None,
    ) -> Path:
        """Create a split-k-mer file from one or more FASTA/FASTQ inputs.

        Returns the path to the produced ``<out_prefix>.skf``.
        """
        if k % 2 == 0:
            raise ValueError(f"k-mer size must be odd (got {k}); split k-mers need a center base")
        out_prefix = Path(out_prefix)
        args: list[str] = ["build", "-o", str(out_prefix), "-k", str(k)]
        if single_strand:
            args.append("--single-strand")
        if min_count is not None:
            args += ["--min-count", str(min_count)]
        if min_qual is not None:
            args += ["--min-qual", str(min_qual)]
        if threads is not None:
            args += ["--threads", str(threads)]
        # A file-list (name<TAB>path per line) lets us hand ska clean sample
        # labels instead of letting it name samples after the verbatim path.
        file_list = out_prefix.with_suffix(".filelist.tsv")
        file_list.write_text(
            "".join(f"{name}\t{path}\n" for name, path in _build_file_list(seq_files)),
            encoding="utf-8",
        )
        args += ["-f", str(file_list)]
        try:
            self._run(args)
        finally:
            file_list.unlink(missing_ok=True)
        return out_prefix.with_suffix(".skf")

    def align(
        self,
        inputs: PathLike | Sequence[PathLike],
        output: PathLike | None = None,
        *,
        min_freq: float = 0.9,
        const_filter: str = "no-const",
        ambig_mask: bool = False,
        threads: int | None = None,
    ) -> Path | str:
        """Write an unordered reference-free SNP alignment.

        ``inputs`` may be a single ``.skf`` path or a list of FASTA files. With
        ``output`` set, returns that path; otherwise returns the alignment text.
        """
        if isinstance(inputs, (str, Path)):
            inputs = [inputs]
        args: list[str] = ["align", "--min-freq", str(min_freq), "--filter", const_filter]
        if ambig_mask:
            args.append("--ambig-mask")
        if threads is not None:
            args += ["--threads", str(threads)]
        if output is not None:
            args += ["-o", str(output)]
        args += [str(p) for p in inputs]
        out = self._run(args)
        return Path(output) if output is not None else out

    def distance(
        self,
        skf: PathLike,
        output: PathLike | None = None,
        *,
        min_freq: float = 0.0,
        allow_ambiguous: bool = False,
        threads: int | None = None,
    ) -> Path | str:
        """Compute pairwise SNP distances from a ``.skf`` file."""
        args: list[str] = ["distance", "--min-freq", str(min_freq)]
        if allow_ambiguous:
            args.append("--allow-ambiguous")
        if threads is not None:
            args += ["--threads", str(threads)]
        if output is not None:
            args += ["-o", str(output)]
        args.append(str(skf))
        out = self._run(args)
        return Path(output) if output is not None else out

    def map(
        self,
        reference: PathLike,
        inputs: PathLike | Sequence[PathLike],
        output: PathLike | None = None,
        *,
        out_format: str = "vcf",
        ambig_mask: bool = False,
        repeat_mask: bool = False,
        threads: int | None = None,
    ) -> Path | str:
        """Map split k-mers onto ``reference`` to get reference-coordinate output.

        Unlike ``align`` (reference-free), ``map`` orders variants against a
        reference so positions can be annotated. ``out_format`` is ``"vcf"`` or
        ``"aln"``. With ``output`` set, returns that path; else the text.
        """
        if out_format not in ("vcf", "aln"):
            raise ValueError(f"out_format must be 'vcf' or 'aln' (got {out_format!r})")
        if isinstance(inputs, (str, Path)):
            inputs = [inputs]
        args: list[str] = ["map", "-f", out_format]
        if ambig_mask:
            args.append("--ambig-mask")
        if repeat_mask:
            args.append("--repeat-mask")
        if threads is not None:
            args += ["--threads", str(threads)]
        if output is not None:
            args += ["-o", str(output)]
        args.append(str(reference))
        args += [str(p) for p in inputs]
        out = self._run(args)
        return Path(output) if output is not None else out

    def merge(self, skf_files: Sequence[PathLike], out_prefix: PathLike) -> Path:
        """Combine multiple ``.skf`` files into one ``<out_prefix>.skf``."""
        out_prefix = Path(out_prefix)
        args: list[str] = ["merge", "-o", str(out_prefix)]
        args += [str(f) for f in skf_files]
        self._run(args)
        return out_prefix.with_suffix(".skf")
