"""Optional maximum-likelihood tree inference by wrapping an external engine.

ML on a SNP alignment needs a real heuristic search engine; skTree shells out to
**IQ-TREE** (preferred) or **RAxML-NG** if either is on PATH, and otherwise
raises :class:`MlNotAvailable` so callers can skip ML cleanly.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

# (binary names to probe, in preference order)
IQTREE_NAMES = ("iqtree2", "iqtree")
RAXML_NAMES = ("raxml-ng",)


class MlNotAvailable(RuntimeError):
    """Raised when no supported ML engine is found on PATH."""


class MlError(RuntimeError):
    """Raised when the ML engine exits non-zero."""


def _which(names: tuple[str, ...]) -> str | None:
    for name in names:
        found = shutil.which(name)
        if found:
            return found
    return None


class MlTreeBuilder:
    """Detects and drives an external ML engine.

    Parameters
    ----------
    engine:
        ``"auto"`` (default) probes IQ-TREE then RAxML-NG. ``"iqtree"`` or
        ``"raxml"`` force a specific engine.
    """

    def __init__(self, engine: str = "auto") -> None:
        if engine in ("auto", "iqtree"):
            iq = _which(IQTREE_NAMES)
            if iq:
                self.engine = "iqtree"
                self.binary = iq
                return
            if engine == "iqtree":
                raise MlNotAvailable("iqtree not found on PATH")
        if engine in ("auto", "raxml"):
            rax = _which(RAXML_NAMES)
            if rax:
                self.engine = "raxml"
                self.binary = rax
                return
        raise MlNotAvailable(
            "No ML engine found. Install IQ-TREE (`iqtree2`) or RAxML-NG (`raxml-ng`)."
        )

    def _argv(
        self,
        alignment: Path,
        prefix: Path,
        model: str | None,
        bootstrap: int | None,
        threads: int | None,
    ) -> tuple[list[str], Path]:
        if self.engine == "iqtree":
            argv = [self.binary, "-s", str(alignment), "--prefix", str(prefix)]
            # No model -> let IQ-TREE run its default ModelFinder Plus (MFP),
            # which selects the best-fitting model itself. Only pin -m when the
            # caller explicitly overrides that default.
            if model:
                argv += ["-m", model]
            argv += ["-nt", str(threads) if threads else "AUTO"]
            if bootstrap:
                argv += ["-B", str(max(bootstrap, 1000))]
            return argv, prefix.with_suffix(".treefile")
        # raxml-ng has no auto-default search, so fall back to a sensible model.
        argv = [self.binary, "--all", "--msa", str(alignment), "--model", model or "GTR+G",
                "--prefix", str(prefix)]
        if threads:
            argv += ["--threads", str(threads)]
        if bootstrap:
            argv += ["--bs-trees", str(bootstrap)]
        return argv, Path(str(prefix) + ".raxml.bestTree")

    def build(
        self,
        alignment: Path,
        prefix: Path,
        *,
        model: str | None = None,
        bootstrap: int | None = None,
        threads: int | None = None,
    ) -> Path:
        """Run the ML engine; return the path to the best tree (Newick).

        ``model`` defaults to ``None``, which runs the engine's own default
        model selection (ModelFinder for IQ-TREE). Pass an explicit model
        string to override it.
        """
        argv, treefile = self._argv(alignment, prefix, model, bootstrap, threads)
        proc = subprocess.run(argv, capture_output=True, text=True)
        if proc.returncode != 0:
            raise MlError(f"{self.engine} exited {proc.returncode}: {proc.stderr or proc.stdout}")
        return treefile
