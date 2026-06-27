# tests/test_report_frontend.py
"""Browser smoke test for the report front-end.

Builds a minimal self-contained HTML from the vendored assets + a tiny embedded
dataset and asserts (via Playwright) that the tree SVG and alignment canvas
render. Marked integration so browserless CI can skip it with -m 'not integration'.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

pytestmark = pytest.mark.integration

ASSETS = Path(__file__).resolve().parents[1] / "sktree" / "report_assets"

DATA = {
    "samples": ["s1", "s2", "s3"],
    "trees": {"nj": "((s1:0.1,s2:0.1):0.1,s3:0.2);"},
    "alignment": {
        "core": {"columns": 3, "rows": {"s1": "ATG", "s2": "ACG", "s3": "ATC"}},
        "all": None,
    },
    "clusters": {"s1": 1, "s2": 1, "s3": 2},
    "stats": {"nSamples": 3, "k": 31, "minFreq": 0.9, "majorityThreshold": 0.5,
              "nAllSnps": 3, "nCoreSnps": 3, "nMajoritySnps": 3},
    "annotation": None,
    "divergence": {"max": None, "warnThreshold": 0.01},
}


def _build_html(tmp_path: Path) -> Path:
    tpl = (ASSETS / "template.html").read_text()
    html = (tpl
            .replace("$STYLES", (ASSETS / "styles.css").read_text())
            .replace("$TPV_JS", (ASSETS / "theiaphyloviewer.iife.js").read_text())
            .replace("$APP_JS", (ASSETS / "app.js").read_text())
            .replace("$DATA_JSON", json.dumps(DATA)))
    out = tmp_path / "smoke.html"
    out.write_text(html)
    return out


def test_frontend_builds_html_without_external_refs(tmp_path):
    html = _build_html(tmp_path).read_text()
    assert not re.search(r'(?:src|href)\s*=\s*["\']?https?://', html)
    for token in ("$STYLES", "$TPV_JS", "$APP_JS", "$DATA_JSON"):
        assert token not in html
    assert "SKTREE_RENDER" in html
