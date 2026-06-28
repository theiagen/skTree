"""Subprocess wrapper around the SKA2 (``ska``) command-line tool.

skTree offloads all k-mer-heavy work to the Rust ``ska`` binary and builds the
analysis layer in Python on top. This module is the single place that knows
how to invoke ``ska`` and turn its results into Python values.

The minimal pipeline is ``build`` (one .skf holding all samples) then ``align``
(a reference-free SNP alignment). ``merge`` exists to combine separately-built
.skf files; ``distance`` yields pairwise SNP distances.
"""

from __future__ import annotations

import logging
import shlex
import shutil
import subprocess
from collections.abc import Sequence
from pathlib import Path

from ..inputs import Sample

PathLike = str | Path

logger = logging.getLogger("sktree.ska")


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
        """Run ``ska <args>``; return stdout. Raise :class:`SkaError` on failure.

        The exact command is logged at ``DEBUG`` before running and, on failure,
        echoed into both the log (with stderr at ``ERROR``) and the raised
        :class:`SkaError`, so a broken run can be reproduced from the log alone.
        """
        argv = [self.binary, *args]
        command = shlex.join(argv)
        logger.debug("running: %s", command)
        proc = subprocess.run(
            argv,
            capture_output=capture,
            text=True,
        )
        if proc.returncode != 0:
            stderr = proc.stderr or proc.stdout or "(no output captured)"
            logger.error("ska failed (exit %d): %s", proc.returncode, command)
            logger.error("ska stderr: %s", stderr.strip())
            raise SkaError(
                f"ska exited with code {proc.returncode}\n"
                f"  command: {command}\n"
                f"  stderr : {stderr.strip()}"
            )
        return proc.stdout

    # -- subcommands ----------------------------------------------------------

    def version(self) -> str:
        return self._run(["--version"]).strip()

    def build(
        self,
        samples: Sequence[Sample],
        out_prefix: PathLike,
        k: int = 31,
        *,
        single_strand: bool = False,
        min_count: int | None = None,
        min_qual: int | None = None,
        qual_filter: str | None = None,
        threads: int | None = None,
    ) -> Path:
        """Create a split-k-mer file from assembly and/or paired-read samples.

        Each sample becomes one file-list line: ``name<TAB>file`` for an
        assembly/single file, or ``name<TAB>fwd<TAB>rev`` for paired reads. The
        read error filters (``--min-count``/``--min-qual``/``--qual-filter``) are
        only forwarded when at least one sample is reads, so pure-assembly builds
        are invoked exactly as before. Returns the produced ``<out_prefix>.skf``.
        """
        if k % 2 == 0:
            raise ValueError(f"k-mer size must be odd (got {k}); split k-mers need a center base")
        out_prefix = Path(out_prefix)
        has_reads = any(s.is_reads for s in samples)
        args: list[str] = ["build", "-o", str(out_prefix), "-k", str(k)]
        if single_strand:
            args.append("--single-strand")
        if has_reads:
            if min_count is not None:
                args += ["--min-count", str(min_count)]
            if min_qual is not None:
                args += ["--min-qual", str(min_qual)]
            if qual_filter is not None:
                args += ["--qual-filter", qual_filter]
        if threads is not None:
            args += ["--threads", str(threads)]
        # A file-list (name<TAB>path[<TAB>path2] per line) lets us hand ska clean
        # sample labels instead of letting it name samples after verbatim paths,
        # and is the only way to declare paired-end reads.
        file_list = out_prefix.with_suffix(".filelist.tsv")
        file_list.write_text(
            "".join(
                "\t".join([s.name, *(str(f) for f in s.files)]) + "\n" for s in samples
            ),
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
