"""Tests for the logging configuration.

The whole point of these tests is the post-mortem story: after any run, a
full-detail log must exist on disk so a failure can be debugged without having
re-run with the right flags.
"""

from __future__ import annotations

import logging
import re

from sktree.logconfig import configure_logging


def test_configure_logging_creates_log_file_in_outdir(tmp_path):
    log_path = configure_logging(tmp_path, verbose=False, debug=False)
    assert log_path == tmp_path / "sktree.log"

    logging.getLogger("sktree").warning("hello from the run")

    assert log_path.exists()
    assert "hello from the run" in log_path.read_text()


def test_file_log_captures_debug_even_when_console_is_quiet(tmp_path):
    # Default verbosity keeps the console at WARNING, but the *file* must still
    # record DEBUG detail -- that is what makes a quiet run debuggable later.
    log_path = configure_logging(tmp_path, verbose=False, debug=False)

    logging.getLogger("sktree.ska").debug("ran: ska build -k 31")

    assert "ran: ska build -k 31" in log_path.read_text()


def test_file_log_lines_carry_a_timestamp(tmp_path):
    log_path = configure_logging(tmp_path, verbose=False, debug=False)
    logging.getLogger("sktree").info("step started")

    first_line = log_path.read_text().splitlines()[0]
    # ISO-ish "2026-06-28 12:00:00" prefix
    assert re.match(r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}", first_line)


def test_configure_logging_is_idempotent(tmp_path):
    # Re-configuring (e.g. a second run in one process) must not duplicate
    # handlers, or every message would be written twice.
    configure_logging(tmp_path, verbose=False, debug=False)
    log_path = configure_logging(tmp_path, verbose=False, debug=False)

    logging.getLogger("sktree").warning("only once")

    assert log_path.read_text().count("only once") == 1
