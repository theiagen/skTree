"""Central logging setup for skTree.

A skTree run shells out to ``ska`` and several optional engines, so when
something fails the useful evidence (the exact command, its stderr, which step
we were on) is easy to lose if it only ever reached the terminal. This module
wires up two sinks:

* a **console** handler whose level follows ``-v`` / ``--debug``, for live use;
* a **file** handler that always records ``DEBUG`` detail to
  ``<outdir>/sktree.log``, so every run leaves a full post-mortem trail on disk
  regardless of how quiet the console was.

Everything logs through the ``"sktree"`` logger (and its ``sktree.*`` children),
so configuring it once here governs the whole package.
"""

from __future__ import annotations

import logging
from pathlib import Path

LOG_FILENAME = "sktree.log"
_ROOT_NAME = "sktree"

_CONSOLE_FORMAT = "%(levelname)s: %(message)s"
_FILE_FORMAT = "%(asctime)s %(levelname)s [%(name)s] %(message)s"
_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"


def configure_logging(outdir: Path, *, verbose: bool = False, debug: bool = False) -> Path:
    """Set up console + persistent file logging for a run; return the log path.

    The file at ``<outdir>/sktree.log`` always captures ``DEBUG`` detail so a
    failed run can be debugged after the fact. The console level is ``DEBUG``
    with ``debug``, ``INFO`` with ``verbose``, otherwise ``WARNING``.

    Safe to call more than once in a process: existing handlers are torn down
    first so messages are never duplicated.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    log_path = outdir / LOG_FILENAME

    logger = logging.getLogger(_ROOT_NAME)
    for handler in list(logger.handlers):
        logger.removeHandler(handler)
        handler.close()

    # The logger passes everything through; each handler applies its own floor.
    logger.setLevel(logging.DEBUG)
    logger.propagate = False

    console_level = logging.DEBUG if debug else logging.INFO if verbose else logging.WARNING
    console = logging.StreamHandler()
    console.setLevel(console_level)
    console.setFormatter(logging.Formatter(_CONSOLE_FORMAT))
    logger.addHandler(console)

    file_handler = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(logging.Formatter(_FILE_FORMAT, datefmt=_DATE_FORMAT))
    logger.addHandler(file_handler)

    return log_path
