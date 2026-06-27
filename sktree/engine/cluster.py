# sktree/engine/cluster.py
"""Optional population-structure clustering by wrapping the ``fastbaps`` CLI.

fastbaps (Bayesian hierarchical clustering of bacterial populations) is an
external Rust/Python tool. skTree shells out to it, mirroring the SKA2 and
ML-engine wrappers, and raises :class:`FastbapsNotAvailable` so callers can skip
clustering cleanly when it is not installed.
"""

from __future__ import annotations

import csv
import shutil
import subprocess
from pathlib import Path


class FastbapsNotAvailable(RuntimeError):
    """Raised when the ``fastbaps`` binary cannot be located on PATH."""


class FastbapsError(RuntimeError):
    """Raised when ``fastbaps`` exits non-zero or writes unusable output."""


class FastbapsRunner:
    """Thin, testable wrapper over the ``fastbaps`` CLI.

    Parameters
    ----------
    binary:
        Name or path of the fastbaps executable. Resolved via
        :func:`shutil.which`; a :class:`FastbapsNotAvailable` is raised if it
        cannot be found.
    """

    def __init__(self, binary: str = "fastbaps") -> None:
        resolved = shutil.which(binary)
        if resolved is None:
            raise FastbapsNotAvailable(
                f"Could not find the '{binary}' executable on PATH. "
                "Install it with `pip install fastbaps` "
                "(https://github.com/thanhleviet/fastbaps-py)."
            )
        self.binary = resolved

    def cluster(
        self,
        alignment: Path,
        output: Path,
        *,
        prior: str = "baps",
        threads: int | None = None,
    ) -> dict[str, int]:
        """Run fastbaps on a FASTA alignment; return ``{sample_name: cluster_id}``."""
        argv = [self.binary, "-i", str(alignment), "-o", str(output), "--prior", prior]
        if threads is not None:
            argv += ["--threads", str(threads)]
        proc = subprocess.run(argv, capture_output=True, text=True)
        if proc.returncode != 0:
            raise FastbapsError(
                f"fastbaps exited with code {proc.returncode}: {proc.stderr or proc.stdout}"
            )
        return self._parse(output)

    @staticmethod
    def _parse(csv_path: Path) -> dict[str, int]:
        if not csv_path.exists():
            raise FastbapsError(f"fastbaps produced no output file at {csv_path}")
        clusters: dict[str, int] = {}
        with csv_path.open(newline="") as fh:
            reader = csv.DictReader(fh)
            if reader.fieldnames != ["Isolates", "Clusters"]:
                raise FastbapsError(
                    f"unexpected fastbaps CSV header: {reader.fieldnames}"
                )
            for row in reader:
                clusters[row["Isolates"]] = int(row["Clusters"])
        if not clusters:
            raise FastbapsError("fastbaps output contained no cluster assignments")
        return clusters
