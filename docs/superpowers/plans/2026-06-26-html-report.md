# skTree HTML Report Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a self-contained interactive HTML report that shows a phylogenetic tree beside the SNP alignment, overlays fastbaps population-structure clusters, and surfaces run stats / divergence warnings / annotation.

**Architecture:** A new `engine/cluster.py` wraps the `fastbaps` CLI as a subprocess (same pattern as `engine/ska.py`/`engine/ml.py`). A new `html_report.py` assembles a `ReportData` object and renders one `.html` with all JS/CSS/data inlined. Vendored `d3` + `phylotree.js` (MIT) plus a custom `app.js` render the tree (SVG) and alignment (`<canvas>`), binding alignment rows to the tree's tip order. `pipeline.py`/`cli.py` wire it in behind `--html`/`--cluster`.

**Tech Stack:** Python 3.10+, numpy, stdlib `subprocess`/`csv`/`json`/`string.Template`/`importlib.resources`; vendored d3 v7 + phylotree.js v2.6 + hand-written `app.js`; pytest + (optional) Playwright MCP for browser smoke tests.

## Global Constraints

- `requires-python = ">=3.10"`; keep `from __future__ import annotations` in every module.
- **No new Python runtime dependency.** HTML templating uses stdlib only (`string.Template`). fastbaps is an external binary, never a pip dependency.
- Subprocess wrappers resolve the binary in `__init__` via `shutil.which` and raise a `*NotAvailable`/`*NotFound` error if absent — mirror `MlTreeBuilder`/`SkaRunner` exactly.
- ruff config: line-length 100, rules `E,F,I,UP,B,W`. Run `ruff check sktree tests` clean before each commit.
- Run tests with `pytest` (config in `pyproject.toml`). Browser tests are marked `integration` and must be deselectable with `-m "not integration"`.
- License is MIT; only vendor MIT/BSD/ISC-licensed JS. Record each vendored file's source URL + version + license in `sktree/report_assets/VENDOR.md`.
- The report MUST be a single self-contained `.html` with **no external `http(s)` `src`/`href`** — everything inlined.
- All code, comments, identifiers in English. Self-documenting code; comment only the non-obvious.

## Embedded JSON schema (shared contract)

`html_report.render_html` writes exactly this object into a `<script id="sktree-data" type="application/json">` block; `app.js` reads exactly this shape. Both tasks depend on it.

```json
{
  "samples": ["s1", "s2"],
  "trees": { "nj": "(s1:0.1,s2:0.1);", "parsimony": "...", "ml": "..." },
  "alignment": {
    "core": { "columns": 3, "rows": { "s1": "ATG", "s2": "ACG" } },
    "all":  { "columns": 5, "rows": { "s1": "ATGCA", "s2": "ACGTA" }, "truncatedTo": null }
  },
  "clusters": { "s1": 1, "s2": 2 },
  "stats": { "nSamples": 2, "k": 31, "minFreq": 0.9, "majorityThreshold": 0.5,
             "nAllSnps": 5, "nCoreSnps": 3, "nMajoritySnps": 4 },
  "annotation": { "coding": 10, "nonsynonymous": 4, "intergenic": 6 },
  "divergence": { "max": 0.012, "warnThreshold": 0.01 }
}
```

- `trees` contains only keys whose tree exists (`nj` always present).
- `alignment.all` is `null` when not computed; `truncatedTo` is an int when capped, else `null`.
- `clusters` is `null` when fastbaps did not run; `annotation` is `null` without `--reference`; `divergence.max` may be `null`.

---

### Task 1: `FastbapsRunner` subprocess wrapper

**Files:**
- Create: `sktree/engine/cluster.py`
- Test: `tests/test_cluster.py`

**Interfaces:**
- Consumes: nothing from other tasks.
- Produces:
  - `class FastbapsNotAvailable(RuntimeError)`
  - `class FastbapsError(RuntimeError)`
  - `FastbapsRunner(binary: str = "fastbaps")` — raises `FastbapsNotAvailable` if `shutil.which(binary)` is None.
  - `FastbapsRunner.cluster(alignment: Path, output: Path, *, prior: str = "baps", threads: int | None = None) -> dict[str, int]`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_cluster.py
"""Unit tests for the fastbaps subprocess wrapper.

Mocks ``subprocess.run`` to assert the argv contract and CSV parsing without
needing the real ``fastbaps`` binary.
"""

from __future__ import annotations

import subprocess

import pytest

from sktree.engine.cluster import FastbapsError, FastbapsNotAvailable, FastbapsRunner

CSV = "Isolates,Clusters\ns1,1\ns2,1\ns3,2\n"


@pytest.fixture
def fake_run(monkeypatch):
    calls = []

    def _fake(argv, **kwargs):
        calls.append((argv, kwargs))
        # The wrapper passes -o <path>; emulate fastbaps writing that file.
        out = argv[argv.index("-o") + 1]
        from pathlib import Path

        Path(out).write_text(CSV)
        return subprocess.CompletedProcess(argv, 0, stdout="done", stderr="")

    monkeypatch.setattr("sktree.engine.cluster.shutil.which", lambda _: "/usr/bin/fastbaps")
    monkeypatch.setattr("sktree.engine.cluster.subprocess.run", _fake)
    return calls


def test_missing_binary_raises(monkeypatch):
    monkeypatch.setattr("sktree.engine.cluster.shutil.which", lambda _: None)
    with pytest.raises(FastbapsNotAvailable):
        FastbapsRunner()


def test_cluster_argv_and_parse(fake_run, tmp_path):
    runner = FastbapsRunner()
    aln = tmp_path / "core_snps.fasta"
    out = tmp_path / "clusters.csv"
    result = runner.cluster(aln, out, prior="baps", threads=4)
    argv = fake_run[0][0]
    assert argv[:1] == ["/usr/bin/fastbaps"]
    assert "-i" in argv and str(aln) in argv
    assert "-o" in argv and str(out) in argv
    assert "--prior" in argv and "baps" in argv
    assert "--threads" in argv and "4" in argv
    assert result == {"s1": 1, "s2": 1, "s3": 2}


def test_cluster_omits_threads_when_none(fake_run, tmp_path):
    runner = FastbapsRunner()
    runner.cluster(tmp_path / "a.fasta", tmp_path / "c.csv")
    argv = fake_run[0][0]
    assert "--threads" not in argv


def test_nonzero_exit_raises(monkeypatch, tmp_path):
    monkeypatch.setattr("sktree.engine.cluster.shutil.which", lambda _: "/usr/bin/fastbaps")

    def _boom(argv, **kwargs):
        return subprocess.CompletedProcess(argv, 1, stdout="", stderr="kaboom")

    monkeypatch.setattr("sktree.engine.cluster.subprocess.run", _boom)
    runner = FastbapsRunner()
    with pytest.raises(FastbapsError, match="kaboom"):
        runner.cluster(tmp_path / "a.fasta", tmp_path / "c.csv")
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_cluster.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'sktree.engine.cluster'`

- [ ] **Step 3: Write the implementation**

```python
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
                "Install fastbaps-py (https://github.com/thanhleviet/fastbaps-py)."
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
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_cluster.py -v && ruff check sktree/engine/cluster.py tests/test_cluster.py`
Expected: all PASS, ruff clean

- [ ] **Step 5: Commit**

```bash
git add sktree/engine/cluster.py tests/test_cluster.py
git commit -m "feat(cluster): wrap fastbaps CLI as optional subprocess engine"
```

---

### Task 2: Vendor the front-end libraries (d3 + phylotree)

**Files:**
- Create: `sktree/report_assets/d3.v7.min.js`
- Create: `sktree/report_assets/phylotree.min.js`
- Create: `sktree/report_assets/VENDOR.md`
- Create (temporary, delete after): `/tmp/claude-1000/-home-ubuntu-skTree/8c640621-d088-4f62-ba3a-39d8679da7a1/scratchpad/vendor_check.html`

**Interfaces:**
- Produces: two vendored JS files that, loaded in order (d3 first), expose globals `d3` and `phylotree` in a browser with no network access.

- [ ] **Step 1: Download the vendored libraries**

```bash
mkdir -p sktree/report_assets
curl -fsSL https://cdn.jsdelivr.net/npm/d3@7/dist/d3.min.js -o sktree/report_assets/d3.v7.min.js
curl -fsSL https://cdn.jsdelivr.net/npm/phylotree@2.6.0/dist/phylotree.min.js -o sktree/report_assets/phylotree.min.js
ls -l sktree/report_assets/
```
Expected: both files non-empty (d3 ~280 KB, phylotree ~tens–hundreds KB).

- [ ] **Step 2: Write a standalone self-containment check harness**

```html
<!-- scratchpad/vendor_check.html : both libs inlined, NO network -->
<!doctype html><html><head><meta charset="utf-8"></head><body>
<div id="out">pending</div>
<script>/*__D3__*/</script>
<script>/*__PHYLO__*/</script>
<script>
  try {
    const t = new phylotree.phylotree("((a:1,b:1):1,c:2);");
    const leaves = t.getTips ? t.getTips() : t.get_nodes().filter(n => !n.children);
    document.getElementById("out").textContent =
      "OK d3=" + (typeof d3) + " phylotree=" + (typeof phylotree) + " tips=" + leaves.length;
  } catch (e) {
    document.getElementById("out").textContent = "ERR " + e.message;
  }
</script></body></html>
```

Build the concrete harness by inlining the two downloaded files in place of the markers:

```bash
python3 - <<'PY'
from pathlib import Path
base = Path("sktree/report_assets")
tpl = Path("/tmp/claude-1000/-home-ubuntu-skTree/8c640621-d088-4f62-ba3a-39d8679da7a1/scratchpad/vendor_check.html").read_text()
tpl = tpl.replace("/*__D3__*/", (base / "d3.v7.min.js").read_text())
tpl = tpl.replace("/*__PHYLO__*/", (base / "phylotree.min.js").read_text())
Path("/tmp/claude-1000/-home-ubuntu-skTree/8c640621-d088-4f62-ba3a-39d8679da7a1/scratchpad/vendor_check_built.html").write_text(tpl)
print("built")
PY
```

- [ ] **Step 3: Verify self-containment in a real browser (Playwright MCP)**

Open `file:///tmp/claude-1000/-home-ubuntu-skTree/8c640621-d088-4f62-ba3a-39d8679da7a1/scratchpad/vendor_check_built.html` with the Playwright MCP browser, snapshot, and read `#out`.
Expected text starts with `OK d3=object phylotree=function tips=3`.

If it shows `ERR ... phylotree is not defined` or a missing-dependency error, phylotree's dist externalizes a peer dep. Resolve by also vendoring the missing dep (commonly `d3` only; if `underscore` is referenced, add `https://cdn.jsdelivr.net/npm/underscore@1/underscore-min.js` as `underscore.min.js` and inline it before phylotree in the harness and later in the template). Re-run until the harness prints `OK ... tips=3`.

- [ ] **Step 4: Record provenance**

```markdown
<!-- sktree/report_assets/VENDOR.md -->
# Vendored front-end assets

| File | Source | Version | License |
|---|---|---|---|
| d3.v7.min.js | https://cdn.jsdelivr.net/npm/d3@7/dist/d3.min.js | 7.x | ISC |
| phylotree.min.js | https://cdn.jsdelivr.net/npm/phylotree@2.6.0/dist/phylotree.min.js | 2.6.0 | MIT |

Verified self-contained (loads offline, exposes `d3` + `phylotree` globals) via
scratchpad/vendor_check_built.html in a Playwright browser.
```

- [ ] **Step 5: Commit**

```bash
git add sktree/report_assets/d3.v7.min.js sktree/report_assets/phylotree.min.js sktree/report_assets/VENDOR.md
git commit -m "chore(report): vendor d3 and phylotree.js for offline reports"
```

---

### Task 3: Front-end (template + styles + app.js) with browser smoke test

**Files:**
- Create: `sktree/report_assets/template.html`
- Create: `sktree/report_assets/styles.css`
- Create: `sktree/report_assets/app.js`
- Create: `tests/test_report_frontend.py`

**Interfaces:**
- Consumes: vendored `d3`/`underscore`/`phylotree` globals (Task 2 — note phylotree needs `underscore` as a peer dependency; the constructor is `new phylotree.phylotree(newick)`). The embedded JSON schema (Global Constraints).
- Produces: `template.html` containing the placeholder tokens `$STYLES`, `$D3_JS`, `$UNDERSCORE_JS`, `$PHYLO_JS`, `$APP_JS`, and `$DATA_JSON` (consumed by Task 5's renderer), and an `app.js` whose `window.SKTREE_RENDER(data)` builds the tree + alignment + cluster strip. **Load order matters: d3, then underscore, then phylotree.**

- [ ] **Step 1: Write the template**

```html
<!-- sktree/report_assets/template.html -->
<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>skTree report</title>
<style>$STYLES</style>
</head>
<body>
<header>
  <h1>skTree report</h1>
  <div id="controls">
    <span id="tree-switcher"></span>
    <label><input type="checkbox" id="all-snps-toggle"> show all variable SNPs</label>
  </div>
</header>
<main>
  <section id="viz">
    <div id="tree"></div>
    <div id="strip"></div>
    <div id="aln-wrap"><canvas id="aln"></canvas></div>
  </section>
  <aside id="panels">
    <div id="stats-panel"></div>
    <div id="warn-panel"></div>
    <div id="annot-panel"></div>
    <div id="legend"></div>
  </aside>
</main>
<div id="tooltip" hidden></div>
<script id="sktree-data" type="application/json">$DATA_JSON</script>
<script>$D3_JS</script>
<script>$UNDERSCORE_JS</script>
<script>$PHYLO_JS</script>
<script>$APP_JS</script>
<script>
  window.SKTREE_RENDER(JSON.parse(document.getElementById("sktree-data").textContent));
</script>
</body>
</html>
```

- [ ] **Step 2: Write the stylesheet**

```css
/* sktree/report_assets/styles.css */
* { box-sizing: border-box; }
body { margin: 0; font: 14px/1.4 system-ui, sans-serif; color: #1a1a1a; }
header { padding: 12px 16px; border-bottom: 1px solid #ddd; }
header h1 { font-size: 18px; margin: 0 0 8px; }
#controls { display: flex; gap: 16px; align-items: center; }
#tree-switcher button { margin-right: 4px; cursor: pointer; }
#tree-switcher button.active { font-weight: 700; text-decoration: underline; }
main { display: flex; gap: 16px; padding: 16px; align-items: flex-start; }
#viz { display: flex; align-items: stretch; }
#tree { flex: 0 0 auto; }
#strip { width: 14px; flex: 0 0 auto; }
#aln-wrap { overflow-x: auto; }
#panels { flex: 0 0 280px; }
#panels > div { margin-bottom: 16px; }
#panels h2 { font-size: 14px; margin: 0 0 6px; }
table.kv { border-collapse: collapse; width: 100%; }
table.kv td { padding: 2px 6px; border-bottom: 1px solid #eee; }
table.kv td:last-child { text-align: right; font-variant-numeric: tabular-nums; }
.warn { color: #b00; }
#legend .swatch { display: inline-block; width: 12px; height: 12px; margin-right: 4px;
  vertical-align: middle; border: 1px solid #999; }
#tooltip { position: fixed; background: #222; color: #fff; padding: 4px 8px;
  border-radius: 4px; font-size: 12px; pointer-events: none; z-index: 10; }
```

- [ ] **Step 3: Write app.js**

```javascript
/* sktree/report_assets/app.js
 * Renders the skTree report: phylotree SVG tree, a <canvas> SNP alignment whose
 * rows are bound to the tree's tip order, and a cluster colour strip between
 * them. Reads the embedded data object (see schema in the plan). */
(function () {
  "use strict";

  var BASE_COLORS = { A: "#4caf50", C: "#2196f3", G: "#ff9800", T: "#e53935" };
  var GAP_COLOR = "#dddddd";
  var CLUSTER_PALETTE = d3.schemeTableau10.concat(d3.schemeSet3);

  function baseColor(b) { return BASE_COLORS[b] || GAP_COLOR; }
  function clusterColor(id) { return CLUSTER_PALETTE[(id - 1) % CLUSTER_PALETTE.length]; }

  function el(id) { return document.getElementById(id); }

  // Read each tip's vertical pixel position from the rendered phylotree SVG by
  // inspecting node-group transforms — robust across phylotree versions.
  function tipYByName(svg) {
    var map = {};
    svg.querySelectorAll("g.node").forEach(function (g) {
      var text = g.querySelector("text");
      if (!text) return;
      var name = text.textContent.trim();
      var tr = g.getAttribute("transform") || "";
      var m = /translate\(\s*([-\d.]+)[ ,]+([-\d.]+)/.exec(tr);
      if (m) map[name] = parseFloat(m[2]);
    });
    return map;
  }

  function renderTree(newick, height) {
    var host = el("tree");
    host.innerHTML = "";
    var tree = new phylotree.phylotree(newick);
    tree.render({
      container: "#tree",
      height: height,
      width: 320,
      "left-right-spacing": "fit-to-size",
      "top-bottom-spacing": "fit-to-size",
      "show-scale": false,
      "align-tips": true,
      "draw-size-bubbles": false
    });
    host.appendChild(tree.display.show());
    return host.querySelector("svg");
  }

  function drawAlignment(data, scope, order, yByName, rowHeight, totalHeight) {
    var block = data.alignment[scope];
    var canvas = el("aln");
    var cols = block.columns;
    var cellW = Math.max(1, Math.min(8, Math.floor(900 / cols)));
    canvas.width = cols * cellW;
    canvas.height = totalHeight;
    var ctx = canvas.getContext("2d");
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    order.forEach(function (name) {
      var seq = block.rows[name];
      if (seq === undefined) return;
      var y = Math.round(yByName[name] - rowHeight / 2);
      for (var j = 0; j < seq.length; j++) {
        ctx.fillStyle = baseColor(seq[j]);
        ctx.fillRect(j * cellW, y, cellW, Math.ceil(rowHeight));
      }
    });
    canvas._cellW = cellW;
    canvas._scope = scope;
  }

  function drawStrip(data, order, yByName, rowHeight, totalHeight) {
    var strip = el("strip");
    if (!data.clusters) { strip.innerHTML = ""; return; }
    var svg = d3.select(strip).html("").append("svg")
      .attr("width", 14).attr("height", totalHeight);
    order.forEach(function (name) {
      var c = data.clusters[name];
      if (c === undefined) return;
      svg.append("rect")
        .attr("x", 0).attr("y", yByName[name] - rowHeight / 2)
        .attr("width", 14).attr("height", Math.ceil(rowHeight))
        .attr("fill", clusterColor(c));
    });
  }

  function panels(data) {
    var s = data.stats;
    el("stats-panel").innerHTML =
      "<h2>Run statistics</h2><table class='kv'>" +
      row("samples", s.nSamples) + row("k", s.k) + row("min-freq", s.minFreq) +
      row("all SNPs", s.nAllSnps) + row("core SNPs", s.nCoreSnps) +
      row("majority SNPs", s.nMajoritySnps) + "</table>";

    var w = el("warn-panel");
    if (data.divergence && data.divergence.max != null &&
        data.divergence.max > data.divergence.warnThreshold) {
      w.innerHTML = "<h2>Warnings</h2><p class='warn'>Max pairwise divergence " +
        (data.divergence.max * 100).toFixed(2) + "% exceeds ~" +
        (data.divergence.warnThreshold * 100).toFixed(0) +
        "%; SKA recall degrades on divergent samples.</p>";
    } else { w.innerHTML = ""; }

    var a = el("annot-panel");
    if (data.annotation) {
      a.innerHTML = "<h2>SNP annotation</h2><table class='kv'>" +
        row("coding", data.annotation.coding) +
        row("non-synonymous", data.annotation.nonsynonymous) +
        row("intergenic", data.annotation.intergenic) + "</table>";
    } else { a.innerHTML = ""; }

    var legend = ["A", "C", "G", "T"].map(function (b) {
      return "<span class='swatch' style='background:" + baseColor(b) + "'></span>" + b;
    }).join(" ");
    legend += " <span class='swatch' style='background:" + GAP_COLOR + "'></span>gap";
    el("legend").innerHTML = "<h2>Bases</h2>" + legend;
  }

  function row(k, v) { return "<tr><td>" + k + "</td><td>" + v + "</td></tr>"; }

  function tooltip(data, scope) {
    var canvas = el("aln"), tip = el("tooltip");
    canvas.onmousemove = function (ev) {
      var rect = canvas.getBoundingClientRect();
      var col = Math.floor((ev.clientX - rect.left) / canvas._cellW);
      tip.hidden = false;
      tip.style.left = (ev.clientX + 12) + "px";
      tip.style.top = (ev.clientY + 12) + "px";
      tip.textContent = scope + " SNP #" + (col + 1) + " / " +
        data.alignment[scope].columns;
    };
    canvas.onmouseleave = function () { el("tooltip").hidden = true; };
  }

  function fullRender(data, scope) {
    var n = data.samples.length;
    var rowHeight = Math.max(8, Math.min(22, 600 / n));
    var totalHeight = Math.ceil(rowHeight * n) + 20;
    var svg = renderTree(data.trees[scope.tree], totalHeight);
    var yByName = tipYByName(svg);
    var order = Object.keys(yByName).sort(function (a, b) {
      return yByName[a] - yByName[b];
    });
    if (order.length === 0) { order = data.samples.slice(); }
    drawStrip(data, order, yByName, rowHeight, totalHeight);
    drawAlignment(data, scope.aln, order, yByName, rowHeight, totalHeight);
    tooltip(data, scope.aln);
  }

  window.SKTREE_RENDER = function (data) {
    var state = { tree: "nj", aln: "core" };
    var sw = el("tree-switcher");
    Object.keys(data.trees).forEach(function (name) {
      var b = document.createElement("button");
      b.textContent = name.toUpperCase();
      if (name === state.tree) b.className = "active";
      b.onclick = function () {
        state.tree = name;
        sw.querySelectorAll("button").forEach(function (x) { x.className = ""; });
        b.className = "active";
        fullRender(data, state);
      };
      sw.appendChild(b);
    });
    var toggle = el("all-snps-toggle");
    if (!data.alignment.all) { toggle.disabled = true; }
    toggle.onchange = function () {
      state.aln = toggle.checked ? "all" : "core";
      fullRender(data, state);
    };
    panels(data);
    fullRender(data, state);
    window.__SKTREE_READY = true;
  };
})();
```

- [ ] **Step 4: Write the browser smoke test**

```python
# tests/test_report_frontend.py
"""Browser smoke test for the report front-end.

Builds a minimal self-contained HTML from the vendored assets + a tiny embedded
dataset and asserts (via Playwright) that the tree SVG and alignment canvas
render. Marked integration so browserless CI can skip it with -m 'not integration'.
"""

from __future__ import annotations

import json
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
            .replace("$D3_JS", (ASSETS / "d3.v7.min.js").read_text())
            .replace("$UNDERSCORE_JS", (ASSETS / "underscore.min.js").read_text())
            .replace("$PHYLO_JS", (ASSETS / "phylotree.min.js").read_text())
            .replace("$APP_JS", (ASSETS / "app.js").read_text())
            .replace("$DATA_JSON", json.dumps(DATA)))
    out = tmp_path / "smoke.html"
    out.write_text(html)
    return out


def test_frontend_builds_html_without_external_refs(tmp_path):
    html = _build_html(tmp_path).read_text()
    assert "src=\"http" not in html and "href=\"http" not in html
    assert "SKTREE_RENDER" in html
```

The pure-Python assertion runs everywhere. The actual browser check is performed manually/by the executing agent with the Playwright MCP: open the built `smoke.html`, wait for `window.__SKTREE_READY === true`, then assert `#tree svg` exists and `#aln` canvas `width > 0`.

- [ ] **Step 5: Run tests + verify in browser**

Run: `pytest tests/test_report_frontend.py -v` (the non-browser assertion passes).
Then, with Playwright MCP: build `smoke.html` via `_build_html`, open it, snapshot, confirm a tree `<svg>` and the `#aln` canvas are present and the cluster strip shows 3 colored rects.
Expected: tree + alignment + strip visible; `window.__SKTREE_READY` true.

- [ ] **Step 6: Commit**

```bash
git add sktree/report_assets/template.html sktree/report_assets/styles.css sktree/report_assets/app.js tests/test_report_frontend.py
git commit -m "feat(report): add interactive front-end (tree + alignment + clusters)"
```

---

### Task 4: `ReportData` assembly + informative-column capping

**Files:**
- Create: `sktree/html_report.py`
- Test: `tests/test_html_report.py`

**Interfaces:**
- Consumes: `RunStats` (`sktree.report`), `SnpMatrix` (`sktree.snps`).
- Produces:
  - `MAX_ALL_SNP_COLUMNS = 20_000`
  - `@dataclass ReportData` with fields: `stats: RunStats`, `sample_names: list[str]`, `trees: dict[str, str]`, `core_matrix: SnpMatrix`, `all_matrix: SnpMatrix | None`, `all_truncated_to: int | None`, `clusters: dict[str, int] | None`, `annotation_summary: dict | None`, `max_divergence: float | None`, `divergence_warn_threshold: float`.
  - `select_informative_columns(m: SnpMatrix, cap: int) -> np.ndarray` (indices of up to `cap` most-informative columns, sorted ascending).
  - `build_report_data(*, stats, trees, core_matrix, all_matrix=None, clusters=None, annotation_summary=None, max_divergence=None, divergence_warn_threshold=0.01) -> ReportData` (applies the cap to `all_matrix`).
  - `to_payload(data: ReportData) -> dict` (produces the embedded JSON schema object).

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_html_report.py
"""Tests for report-data assembly and the embedded JSON payload."""

from __future__ import annotations

import numpy as np

from sktree.html_report import (
    MAX_ALL_SNP_COLUMNS,
    build_report_data,
    select_informative_columns,
    to_payload,
)
from sktree.report import build_stats
from sktree.snps import parse_alignment

ALN = ">s1\nAAAGA\n>s2\nACCGC\n>s3\nAA--T\n>s4\nACC--\n"


def _stats_and_matrix():
    m = parse_alignment(ALN)
    cls = m.classify()
    stats = build_stats(m, cls, k=31, min_freq=0.9, majority_threshold=0.5)
    return stats, m, cls


def test_select_informative_columns_caps_and_sorts():
    m = parse_alignment(">a\nACGTACGT\n>b\nAGGTTCGA\n>c\nACGAACGT\n")
    idx = select_informative_columns(m, cap=3)
    assert len(idx) == 3
    assert list(idx) == sorted(idx)


def test_build_report_data_caps_all_matrix():
    stats, m, cls = _stats_and_matrix()
    big = m.subset(np.ones(m.n_loci, dtype=bool))
    data = build_report_data(
        stats=stats, trees={"nj": "(s1,s2,(s3,s4));"},
        core_matrix=m.subset(cls.core), all_matrix=big,
    )
    # tiny here, so no truncation
    assert data.all_truncated_to is None
    assert data.all_matrix is not None


def test_to_payload_shape():
    stats, m, cls = _stats_and_matrix()
    data = build_report_data(
        stats=stats, trees={"nj": "(s1,s2,(s3,s4));"},
        core_matrix=m.subset(cls.core),
        clusters={"s1": 1, "s2": 1, "s3": 2, "s4": 2},
        annotation_summary={"coding": 3, "nonsynonymous": 1, "intergenic": 2},
        max_divergence=0.02, divergence_warn_threshold=0.01,
    )
    payload = to_payload(data)
    assert payload["samples"] == ["s1", "s2", "s3", "s4"]
    assert "nj" in payload["trees"]
    assert payload["alignment"]["core"]["columns"] == data.core_matrix.n_loci
    assert payload["alignment"]["all"] is None
    assert payload["clusters"] == {"s1": 1, "s2": 1, "s3": 2, "s4": 2}
    assert payload["annotation"]["coding"] == 3
    assert payload["divergence"]["max"] == 0.02
    assert set(payload["alignment"]["core"]["rows"]) == {"s1", "s2", "s3", "s4"}


def test_max_all_snp_columns_default():
    assert MAX_ALL_SNP_COLUMNS == 20_000
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_html_report.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'sktree.html_report'`

- [ ] **Step 3: Write the implementation**

```python
# sktree/html_report.py
"""Assemble run artifacts into a self-contained interactive HTML report.

This module is split into two halves: data assembly (`ReportData`, the embedded
JSON payload) lives here; the rendering half (templating + asset inlining) is
added in the same module by a later task. It is pure given its inputs.
"""

from __future__ import annotations

from dataclasses import dataclass, field

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
    }
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_html_report.py -v && ruff check sktree/html_report.py tests/test_html_report.py`
Expected: all PASS, ruff clean

- [ ] **Step 5: Commit**

```bash
git add sktree/html_report.py tests/test_html_report.py
git commit -m "feat(report): assemble ReportData and embedded JSON payload"
```

---

### Task 5: Render the self-contained HTML (inline template + assets + data)

**Files:**
- Modify: `sktree/html_report.py` (append rendering functions)
- Test: `tests/test_html_report.py` (append render tests)

**Interfaces:**
- Consumes: `to_payload` / `ReportData` (Task 4); `report_assets/*` (Tasks 2–3); placeholder tokens in `template.html` (Task 3).
- Produces:
  - `render_html(data: ReportData) -> str`
  - `write_html_report(path: Path, data: ReportData) -> None`

- [ ] **Step 1: Write the failing tests (append to tests/test_html_report.py)**

```python
# append to tests/test_html_report.py
from pathlib import Path

from sktree.html_report import render_html, write_html_report


def _data():
    stats, m, cls = _stats_and_matrix()
    return build_report_data(
        stats=stats, trees={"nj": "(s1,s2,(s3,s4));"},
        core_matrix=m.subset(cls.core),
        clusters={"s1": 1, "s2": 1, "s3": 2, "s4": 2},
    )


def test_render_html_is_self_contained():
    html = render_html(_data())
    assert "<canvas id=\"aln\">" in html
    assert "SKTREE_RENDER" in html
    # No external resources: proves single-file portability.
    assert "src=\"http" not in html
    assert "href=\"http" not in html
    # Embedded data + tree present.
    assert "(s1,s2,(s3,s4));" in html
    assert "\"s3\": 2" in html or "\"s3\":2" in html


def test_render_html_graceful_without_clusters_or_annotation():
    stats, m, cls = _stats_and_matrix()
    data = build_report_data(stats=stats, trees={"nj": "(s1,s2,(s3,s4));"},
                             core_matrix=m.subset(cls.core))
    html = render_html(data)
    assert "\"clusters\": null" in html or "\"clusters\":null" in html
    assert "\"annotation\": null" in html or "\"annotation\":null" in html


def test_write_html_report(tmp_path):
    out = tmp_path / "report.html"
    write_html_report(out, _data())
    assert out.exists() and out.stat().st_size > 1000
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_html_report.py -k "render or write" -v`
Expected: FAIL with `ImportError: cannot import name 'render_html'`

- [ ] **Step 3: Append the rendering implementation to sktree/html_report.py**

```python
# append to sktree/html_report.py

import json
from importlib import resources
from pathlib import Path
from string import Template


def _asset(name: str) -> str:
    return (resources.files("sktree") / "report_assets" / name).read_text(encoding="utf-8")


def render_html(data: ReportData) -> str:
    """Return a single self-contained HTML document for ``data``."""
    template = Template(_asset("template.html"))
    return template.safe_substitute(
        STYLES=_asset("styles.css"),
        D3_JS=_asset("d3.v7.min.js"),
        UNDERSCORE_JS=_asset("underscore.min.js"),
        PHYLO_JS=_asset("phylotree.min.js"),
        APP_JS=_asset("app.js"),
        DATA_JSON=json.dumps(to_payload(data)),
    )


def write_html_report(path: Path, data: ReportData) -> None:
    path.write_text(render_html(data), encoding="utf-8")
```

Note: add `import json`, `from importlib import resources`, `from pathlib import Path`, and `from string import Template` to the existing import block at the top of the module instead of leaving them inline (ruff `E402`/`I` will flag mid-file imports). Move them up when applying.

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_html_report.py -v && ruff check sktree/html_report.py tests/test_html_report.py`
Expected: all PASS, ruff clean

- [ ] **Step 5: Browser end-to-end check (Playwright MCP)**

Write a built report to a temp file via `write_html_report` using the Task-3 `DATA` shape (with a real `nj` newick), open it in the Playwright browser, wait for `window.__SKTREE_READY`, assert a `#tree svg` and a non-zero-width `#aln` canvas render and the cluster strip shows colored rects.
Expected: tree + alignment + strip visible.

- [ ] **Step 6: Commit**

```bash
git add sktree/html_report.py tests/test_html_report.py
git commit -m "feat(report): render self-contained HTML with inlined assets"
```

---

### Task 6: Wire clustering + report into the pipeline

**Files:**
- Modify: `sktree/pipeline.py`
- Test: `tests/test_pipeline.py` (append)

**Interfaces:**
- Consumes: `FastbapsRunner`/`FastbapsNotAvailable` (Task 1); `build_report_data`/`write_html_report` (Tasks 4–5); existing `parse_alignment`, `RunResult`.
- Produces: `run_pipeline(..., html: bool = False, cluster: bool = False)` and `RunResult` gains `html_report: Path | None = None` and `clusters: Path | None = None`.

- [ ] **Step 1: Write the failing test (append to tests/test_pipeline.py)**

```python
# append to tests/test_pipeline.py
from pathlib import Path

from sktree.pipeline import _build_clusters


def test_build_clusters_skips_when_unavailable(monkeypatch, tmp_path):
    def _raise(*a, **k):
        from sktree.engine.cluster import FastbapsNotAvailable
        raise FastbapsNotAvailable("nope")

    monkeypatch.setattr("sktree.pipeline.FastbapsRunner", _raise)
    aln = tmp_path / "core_snps.fasta"
    aln.write_text(">s1\nAC\n>s2\nAT\n")
    csv_path, clusters = _build_clusters(aln, tmp_path, threads=None, write_csv=True)
    assert csv_path is None and clusters is None
```

If `tests/test_pipeline.py` mocks the SKA runner for full-pipeline tests, follow that existing fixture pattern for any end-to-end additions; this unit test targets the new helper directly.

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_pipeline.py -k build_clusters -v`
Expected: FAIL with `ImportError: cannot import name '_build_clusters'`

- [ ] **Step 3: Implement the wiring in sktree/pipeline.py**

Add imports near the existing ones:

```python
from .engine.cluster import FastbapsError, FastbapsNotAvailable, FastbapsRunner
from .html_report import build_report_data, write_html_report
```

Extend `RunResult` with two fields (after `annotation`):

```python
    clusters: Path | None = None
    html_report: Path | None = None
```

Add the clustering helper:

```python
def _build_clusters(
    core_alignment: Path,
    outdir: Path,
    *,
    threads: int | None,
    write_csv: bool,
) -> tuple[Path | None, dict[str, int] | None]:
    """Run fastbaps on the core alignment; skip gracefully if unavailable."""
    try:
        runner = FastbapsRunner()
    except FastbapsNotAvailable as exc:
        logger.warning("Skipping clustering: %s", exc)
        return None, None
    out = outdir / "clusters.csv"
    try:
        clusters = runner.cluster(core_alignment, out, threads=threads)
    except FastbapsError as exc:
        logger.warning("fastbaps failed: %s", exc)
        return None, None
    logger.info("fastbaps assigned %d clusters", len(set(clusters.values())))
    return (out if write_csv else None), clusters
```

Add `html: bool = False, cluster: bool = False` to `run_pipeline`'s signature (after `reference`). Before building `stats`/`summary` at the end, assemble clusters + annotation summary + trees and write the report. Insert after the annotation block and before `stats = build_stats(...)`:

```python
    clusters_path: Path | None = None
    clusters_map: dict[str, int] | None = None
    if cluster or html:
        clusters_path, clusters_map = _build_clusters(
            core_alignment, outdir, threads=threads, write_csv=cluster
        )
```

After `stats = build_stats(...)` and `summary.write_text(...)`, add the report build:

```python
    html_path: Path | None = None
    if html:
        trees: dict[str, str] = {"nj": nj_tree.read_text().strip()}
        if parsimony_path is not None:
            trees["parsimony"] = parsimony_path.read_text().strip()
        if ml_path is not None:
            trees["ml"] = ml_path.read_text().strip()
        annotation_summary = None
        if annotation_path is not None:
            annotation_summary = _annotation_counts(annotation_path)
        report_data = build_report_data(
            stats=stats,
            trees=trees,
            core_matrix=core,
            all_matrix=m.subset(cls.variable),
            clusters=clusters_map,
            annotation_summary=annotation_summary,
            max_divergence=max_div,
            divergence_warn_threshold=DIVERGENCE_WARN_THRESHOLD,
        )
        html_path = outdir / "report.html"
        write_html_report(html_path, report_data)
        logger.info("Wrote HTML report to %s", html_path.name)
```

Add the annotation-count helper (reads the TSV written by `write_annotation_tsv`):

```python
def _annotation_counts(annotation_tsv: Path) -> dict[str, int]:
    coding = nonsyn = intergenic = 0
    lines = annotation_tsv.read_text().splitlines()
    header = lines[0].split("\t") if lines else []
    try:
        region_i = header.index("region")
        effect_i = header.index("effect")
    except ValueError:
        return {"coding": 0, "nonsynonymous": 0, "intergenic": 0}
    for line in lines[1:]:
        fields = line.split("\t")
        if len(fields) <= max(region_i, effect_i):
            continue
        if fields[region_i] == "coding":
            coding += 1
            if fields[effect_i] == "nonsynonymous":
                nonsyn += 1
        else:
            intergenic += 1
    return {"coding": coding, "nonsynonymous": nonsyn, "intergenic": intergenic}
```

Finally, pass the two new paths into the returned `RunResult(...)`:

```python
        clusters=clusters_path,
        html_report=html_path,
```

Before implementing `_annotation_counts`, open `sktree/annotate.py` and confirm the exact column names emitted by `write_annotation_tsv` (the header). If they differ from `region`/`effect`, use the actual names. (The annotate notes describe per-SNP `region` = intergenic/coding and `effect` = synonymous/nonsynonymous.)

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_pipeline.py -v && ruff check sktree/pipeline.py tests/test_pipeline.py`
Expected: all PASS, ruff clean

- [ ] **Step 5: Commit**

```bash
git add sktree/pipeline.py tests/test_pipeline.py
git commit -m "feat(pipeline): wire fastbaps clustering and HTML report into run"
```

---

### Task 7: CLI flags, packaging, and docs

**Files:**
- Modify: `sktree/cli.py`
- Modify: `pyproject.toml`
- Modify: `README.md`, `FOR-DEVELOPERS.md`
- Test: `tests/test_cli.py` (create if absent, else append)

**Interfaces:**
- Consumes: `run_pipeline(..., html=, cluster=)` and the new `RunResult` fields (Task 6).
- Produces: `sktree run --html` and `sktree run --cluster` flags; vendored assets shipped in the wheel.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_cli.py  (append if file exists)
"""CLI argument wiring tests."""

from __future__ import annotations

from sktree.cli import _build_parser


def test_run_parser_has_html_and_cluster_flags():
    parser = _build_parser()
    args = parser.parse_args(["run", "a.fa", "b.fa", "-o", "out", "--html", "--cluster"])
    assert args.html is True
    assert args.cluster is True


def test_html_and_cluster_default_false():
    parser = _build_parser()
    args = parser.parse_args(["run", "a.fa", "b.fa", "-o", "out"])
    assert args.html is False
    assert args.cluster is False
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_cli.py -v`
Expected: FAIL with `AttributeError: 'Namespace' object has no attribute 'html'`

- [ ] **Step 3: Add the flags and pass them through (sktree/cli.py)**

After the `--reference` argument in `_build_parser`, add:

```python
    run.add_argument(
        "--cluster", action="store_true",
        help="assign population-structure clusters with fastbaps (writes clusters.csv)",
    )
    run.add_argument(
        "--html", action="store_true",
        help="write a self-contained interactive HTML report (runs fastbaps if available)",
    )
```

In the `run_pipeline(...)` call inside `main`, add the two kwargs:

```python
                html=args.html,
                cluster=args.cluster,
```

After the existing `print` lines, add reporting of the new artifacts:

```python
        if result.clusters:
            print(f"  clusters  : {result.clusters.name}")
        if result.html_report:
            print(f"  report    : {result.html_report.name}")
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_cli.py -v`
Expected: PASS

- [ ] **Step 5: Ensure assets ship in the wheel**

Confirm hatchling includes the non-Python `report_assets/` files in the built wheel. Add to `pyproject.toml` under the wheel target:

```toml
[tool.hatch.build.targets.wheel.force-include]
"sktree/report_assets" = "sktree/report_assets"
```

Verify:

```bash
python -m build --wheel 2>/dev/null && \
  python - <<'PY'
import zipfile, glob
whl = sorted(glob.glob("dist/*.whl"))[-1]
names = zipfile.ZipFile(whl).namelist()
assert any(n.endswith("report_assets/app.js") for n in names), "app.js missing from wheel"
assert any(n.endswith("report_assets/phylotree.min.js") for n in names), "phylotree missing"
print("assets present in wheel")
PY
```
Expected: `assets present in wheel`. (If `python -m build` is unavailable, install `build` into `.venv` with `uv pip install build`.)

- [ ] **Step 6: Update docs**

In `README.md`, add `--html` / `--cluster` to the usage/flags section with a one-line example:

```bash
sktree run *.fasta -o out --html        # interactive report at out/report.html
```

In `FOR-DEVELOPERS.md`, add a short "HTML report" subsection covering: the `engine/cluster.py` fastbaps wrapper (optional, graceful skip), `html_report.py` assembly + rendering, the vendored `report_assets/` (offline, MIT/ISC), the canvas-bound-to-tip-y design, and the 20k-column cap with the lesson that a portable single-file report must inline everything (no CDN) so it survives air-gapped clusters and link rot.

- [ ] **Step 7: Full suite + lint**

Run: `pytest -m "not integration" -q && ruff check sktree tests`
Expected: all PASS, ruff clean

- [ ] **Step 8: Commit**

```bash
git add sktree/cli.py pyproject.toml README.md FOR-DEVELOPERS.md tests/test_cli.py
git commit -m "feat(cli): add --html and --cluster flags; ship report assets"
```

---

## Self-Review

**Spec coverage:**
- Self-contained interactive single HTML → Tasks 3 (front-end) + 5 (render/inline). ✓
- Vendored phylotree.js + canvas alignment + thin app.js → Tasks 2, 3. ✓
- fastbaps optional subprocess CLI, graceful degradation → Task 1 + Task 6 `_build_clusters`. ✓
- All trees (NJ default) + core/all toggle + stats/warnings/annotation panels → Task 3 `app.js`, payload in Task 4, wiring in Task 6. ✓
- 20k-column cap, most-informative selection, logged → Task 4 `select_informative_columns` + Task 6 log. ✓
- `--html` auto-runs fastbaps if available; `--cluster` writes standalone CSV → Task 6 (`cluster or html`, `write_csv=cluster`) + Task 7 flags. ✓
- Vendoring committed + provenance → Task 2 `VENDOR.md`. ✓
- No external refs assertion → Tasks 3 + 5 tests. ✓
- Tests: cluster, html_report, frontend, pipeline, cli → Tasks 1, 3, 4, 5, 6, 7. ✓
- Optional Playwright smoke (gated `integration`) → Tasks 3, 5. ✓
- Out of scope items (metadata upload, multi-run, radial) → not implemented. ✓

**Placeholder scan:** No "TBD"/"handle edge cases"/"similar to" — all steps carry concrete code or exact commands. The only deferred verification (annotate.py column names, phylotree peer deps) is written as an explicit inspect-then-adapt step with a concrete fallback, not a vague placeholder.

**Type consistency:** `build_report_data`/`to_payload`/`ReportData` field names match across Tasks 4–6; payload keys match `app.js` reads in Task 3 (`alignment.core.rows`, `clusters`, `stats.nCoreSnps`, `divergence.warnThreshold`); `_build_clusters` return tuple `(Path|None, dict|None)` matches its Task-6 caller and Task-6 test; template tokens `$STYLES/$D3_JS/$UNDERSCORE_JS/$PHYLO_JS/$APP_JS/$DATA_JSON` defined in Task 3 match `render_html` in Task 5.
