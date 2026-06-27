"""Assemble run artifacts into a self-contained interactive HTML report.

This module is split into two halves: data assembly (`ReportData`, the embedded
JSON payload) lives here; the rendering half (templating + asset inlining) is
added in the same module by a later task. It is pure given its inputs.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from importlib import resources
from pathlib import Path
from string import Template

import numpy as np

from .report import RunStats
from .snps import VALID_BASES, SnpMatrix

MAX_ALL_SNP_COLUMNS = 20_000


@dataclass
class ReportData:
    stats: RunStats
    sample_names: list[str]
    trees: dict[str, str]
    core_matrix: SnpMatrix
    all_matrix: SnpMatrix | None = None
    all_truncated_to: int | None = None
    clusters: dict[str, int] | None = None
    annotation_summary: dict | None = None
    max_divergence: float | None = None
    divergence_warn_threshold: float = 0.01
    outputs: list[dict] | None = None


def select_informative_columns(m: SnpMatrix, cap: int) -> np.ndarray:
    """Indices of up to ``cap`` most-informative columns (sorted ascending).

    Informativeness = minor-allele presence: how many samples carry a base
    other than the column's most common unambiguous base. Ties broken by column
    order; the returned indices are sorted so the alignment stays left-to-right.
    """
    present = np.isin(m.matrix, VALID_BASES)
    scores = np.zeros(m.n_loci, dtype=np.int64)
    for j in range(m.n_loci):
        col = m.matrix[present[:, j], j]
        if col.size == 0:
            continue
        _, counts = np.unique(col, return_counts=True)
        scores[j] = int(col.size - counts.max())
    if m.n_loci <= cap:
        return np.arange(m.n_loci)
    top = np.argpartition(scores, m.n_loci - cap)[m.n_loci - cap:]
    return np.sort(top)


def build_report_data(
    *,
    stats: RunStats,
    trees: dict[str, str],
    core_matrix: SnpMatrix,
    all_matrix: SnpMatrix | None = None,
    clusters: dict[str, int] | None = None,
    annotation_summary: dict | None = None,
    max_divergence: float | None = None,
    divergence_warn_threshold: float = 0.01,
    outputs: list[dict] | None = None,
) -> ReportData:
    truncated_to: int | None = None
    if all_matrix is not None and all_matrix.n_loci > MAX_ALL_SNP_COLUMNS:
        idx = select_informative_columns(all_matrix, MAX_ALL_SNP_COLUMNS)
        all_matrix = all_matrix.subset(_mask_from_indices(all_matrix.n_loci, idx))
        truncated_to = MAX_ALL_SNP_COLUMNS
    return ReportData(
        stats=stats,
        sample_names=list(core_matrix.sample_names),
        trees=dict(trees),
        core_matrix=core_matrix,
        all_matrix=all_matrix,
        all_truncated_to=truncated_to,
        clusters=clusters,
        annotation_summary=annotation_summary,
        max_divergence=max_divergence,
        divergence_warn_threshold=divergence_warn_threshold,
        outputs=outputs,
    )


def _mask_from_indices(n: int, idx: np.ndarray) -> np.ndarray:
    mask = np.zeros(n, dtype=bool)
    mask[idx] = True
    return mask


def _matrix_block(m: SnpMatrix | None) -> dict | None:
    if m is None:
        return None
    rows = {name: "".join(row) for name, row in zip(m.sample_names, m.matrix, strict=True)}
    return {"columns": int(m.n_loci), "rows": rows}


def to_payload(data: ReportData) -> dict:
    s = data.stats
    all_block = _matrix_block(data.all_matrix)
    if all_block is not None:
        all_block["truncatedTo"] = data.all_truncated_to
    return {
        "samples": list(data.sample_names),
        "trees": dict(data.trees),
        "alignment": {"core": _matrix_block(data.core_matrix), "all": all_block},
        "clusters": data.clusters,
        "stats": {
            "nSamples": s.n_samples, "k": s.k, "minFreq": s.min_freq,
            "majorityThreshold": s.majority_threshold, "nAllSnps": s.n_all_snps,
            "nCoreSnps": s.n_core_snps, "nMajoritySnps": s.n_majority_snps,
        },
        "annotation": data.annotation_summary,
        "divergence": {
            "max": data.max_divergence,
            "warnThreshold": data.divergence_warn_threshold,
        },
        "outputs": data.outputs,
    }


def _asset(name: str) -> str:
    return (resources.files("sktree") / "report_assets" / name).read_text(encoding="utf-8")


def render_html(data: ReportData) -> str:
    """Return a single self-contained HTML document for ``data``."""
    template = Template(_asset("template.html"))
    data_json = json.dumps(to_payload(data))
    data_json = data_json.replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
    return template.safe_substitute(
        STYLES=_asset("styles.css"),
        TPV_JS=_asset("theiaphyloviewer.iife.js"),
        APP_JS=_asset("app.js"),
        DATA_JSON=data_json,
    )


def write_html_report(path: Path, data: ReportData) -> None:
    path.write_text(render_html(data), encoding="utf-8")
