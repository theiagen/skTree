# skTree HTML Report — Design

**Date:** 2026-06-26
**Status:** Approved (brainstorming) → ready for implementation plan

## 1. Goal

Add an advanced, self-contained HTML report to skTree that places an interactive
phylogenetic tree **beside** the SNP alignment, overlays population-structure
**clusters** from [fastbaps-py](https://github.com/thanhleviet/fastbaps-py), and
surfaces the run's other relevant information (stats, divergence warnings, SNP
annotation). The deliverable is a single portable `.html` file that works
offline with no server.

## 2. Decisions (locked)

| Topic | Decision |
|---|---|
| Rendering | Self-contained **interactive** single `.html`; data embedded as JSON, rendered by bundled JS. Canvas alignment for scale. |
| Front-end | **Vendored libraries**: phylotree.js (MIT, D3-based, rectangular SVG) for the tree + a custom `<canvas>` alignment track, glued by a thin `app.js`. All inlined at render time. |
| fastbaps integration | **Optional subprocess CLI** (`fastbaps -i core_snps.fasta -o clusters.csv`), mirroring the existing `engine/ska.py` / `engine/ml.py` pattern. Degrades gracefully when not installed. |
| Content | All available trees (NJ default, switcher for Parsimony/ML when present); core-SNP alignment by default with an "all variable SNPs" toggle; panels for stats, divergence warnings, annotation summary. |
| All-SNP cap | Embed core SNPs fully; cap the all-SNP view at **~20 000 columns**, selecting the most-informative columns (highest minor-allele presence); log when truncated. |
| Clustering trigger | `--html` auto-runs fastbaps **if on PATH** (graceful skip otherwise); `--cluster` additionally writes the standalone `clusters.csv`. |
| Vendoring | Commit `d3.min.js` + `phylotree.min.js` into `report_assets/` (~300 KB) so the report is genuinely offline and reproducible. |

## 3. Architecture & module boundaries

Five focused units, each with one clear purpose and a well-defined interface:

| Unit | Responsibility | Depends on |
|---|---|---|
| `sktree/engine/cluster.py` | `FastbapsRunner` — subprocess wrapper around `fastbaps`; parses `Isolates,Clusters` CSV → `dict[name, int]`. Raises `FastbapsNotAvailable` / `FastbapsError` (loud-fail like `engine/ml.py`). | `subprocess`, `csv` |
| `sktree/html_report.py` | **New.** Assemble a `ReportData` object from stats/trees/alignment/clusters/annotation, then render a single self-contained `.html` via a template with all JS/CSS inlined. Pure given its inputs. | `report_assets/`, stdlib templating |
| `sktree/report_assets/` | **New.** Vendored static front-end: `template.html`, `d3.min.js`, `phylotree.min.js`, `app.js`, `styles.css`. Committed → zero network at render time. | — |
| `sktree/report.py` | **Unchanged.** Stays the text/TSV summary writer (`RunStats`, `write_snp_matrix_tsv`, `format_summary`). Keeps that file single-purpose. | — |
| `sktree/pipeline.py` / `sktree/cli.py` | Wire clustering + report into the run; add `--html` / `--cluster` flags; extend `RunResult`. | above |

Putting HTML rendering in a **new** `html_report.py` rather than swelling
`report.py` keeps each file answerable to one question — consistent with the
existing `snps.py` / `phylo.py` / `annotate.py` separation.

## 4. `engine/cluster.py` contract

```python
class FastbapsNotAvailable(RuntimeError): ...   # binary not on PATH
class FastbapsError(RuntimeError): ...          # non-zero exit / bad output

class FastbapsRunner:
    def __init__(self, binary: str = "fastbaps"): ...
    # Raises FastbapsNotAvailable in __init__ if shutil.which(binary) is None,
    # mirroring MlTreeBuilder.
    def cluster(
        self,
        alignment: Path,
        output: Path,
        *,
        prior: str = "baps",
        threads: int | None = None,
    ) -> dict[str, int]:
        """Run fastbaps on a FASTA alignment; return {sample_name: cluster_id}."""
```

- Single `_run` choke point for subprocess, matching `SkaRunner`.
- Parses the CSV header `Isolates,Clusters`; tolerant of trailing newline.
- Cluster ids are integers as emitted by fastbaps (1-indexed, contiguous).

## 5. `html_report.py` contract

```python
@dataclass
class ReportData:
    stats: RunStats
    sample_names: list[str]
    trees: dict[str, str]            # {"nj": newick, "parsimony": ..., "ml": ...}
    core_alignment: SnpMatrix        # rows aligned to sample_names
    all_alignment: SnpMatrix | None  # variable SNPs, may be column-capped
    all_truncated_to: int | None     # column count if capped, else None
    clusters: dict[str, int] | None  # {sample: cluster_id} or None
    annotation_summary: dict | None  # counts: coding / nonsyn / intergenic
    max_divergence: float | None
    divergence_warn_threshold: float

def build_report_data(...) -> ReportData: ...
def render_html(data: ReportData) -> str: ...      # returns full HTML string
def write_html_report(path: Path, data: ReportData) -> None: ...
```

**Embedding format.** Alignment is embedded as one string per sample
(`["ATGC...", ...]`) keyed by sample name, plus the column count — compact and
trivially mapped to canvas cells. Trees embed as newick strings. Clusters embed
as a `{sample: id}` map. All data goes into one `<script type="application/json">`
block; assets are inlined verbatim into the template.

**Scale guard.** `core_alignment` embeds fully (typically hundreds–few thousand
columns). `all_alignment` is capped at `MAX_ALL_SNP_COLUMNS = 20_000`; when the
variable-SNP count exceeds the cap, select the columns with the highest
minor-allele presence and set `all_truncated_to`; the renderer shows a note and
the pipeline logs the truncation.

## 6. Front-end (`app.js`) — the tree↔alignment bind

1. **Tree** — phylotree.js renders the chosen newick as a rectangular,
   branch-length-scaled SVG. A switcher exposes only the trees present in
   `data.trees`.
2. **Alignment** — a `<canvas>` grid: one row per sample, one column per SNP,
   cells colored by base (A/C/G/T/gap legend). Canvas keeps it smooth at
   thousands of columns where DOM cells would not.
3. **The bind** — after phylotree lays out, read each tip's `y` pixel
   coordinate; draw canvas rows in that exact order and vertical position.
   Re-ladderizing the tree re-reads tip `y` and redraws the canvas, so the two
   panels never desync.
4. **Cluster strip** — a color strip in the gutter between tree and alignment,
   colored by `data.clusters`; omitted entirely when clusters are absent.
5. **Interactivity** — hover tooltip (sample · SNP position · base · cluster),
   zoom/pan, tree switcher, core/all-SNP toggle, base + cluster legends.

## 7. Data flow

```
inputs → ska → alignment.fasta → SnpMatrix → classify → core_snps.fasta
                                                  │
        ┌─────────────────────────────────────────┼───────────────┐
   NJ/[pars]/[ml] newick          fastbaps clusters.csv      [snp_annotation.tsv]
        └─────────────────────────────────────────┼───────────────┘
                                          ReportData (JSON-embedded)
                                                  │
                                  html_report.render_html() → report.html
```

## 8. Pipeline & CLI wiring

- `pipeline.run_pipeline` gains `html: bool = False` and `cluster: bool = False`.
  A `_cluster()` helper (shaped like `_annotate_snps`) runs fastbaps when
  requested/available and returns the `clusters.csv` path + map.
- After trees (and optional clustering/annotation) are built, assemble
  `ReportData` and write `report.html` when `--html`.
- `RunResult` gains `html_report: Path | None` and `clusters: Path | None`.
- `cli.py`: add `--html` (emit `report.html`; auto-run fastbaps if available)
  and `--cluster` (also write standalone `clusters.csv`).

## 9. Error handling (graceful degradation)

- `fastbaps` not on PATH → warn; render report **without** the cluster strip.
- annotation absent (no `--reference`) → omit that panel.
- only NJ present → switcher shows just NJ.
- all-SNP view over cap → truncate to most-informative columns + log.

The report always renders so long as the NJ tree exists — which the pipeline
always produces.

## 10. Testing

- `tests/test_cluster.py` — mock subprocess: CSV parsing, `FastbapsNotAvailable`
  when binary missing, `prior`/`threads` passthrough, bad-output → `FastbapsError`.
- `tests/test_report.py` (extend) — build HTML from synthetic `ReportData`:
  assert embedded newick + sample names + cluster values present; assert **no
  external `http(s)` `src`/`href`** (proves self-contained); assert graceful
  render when `clusters` / `annotation_summary` are `None`; assert all-SNP cap
  truncates and records `all_truncated_to`.
- Optional Playwright smoke test (Playwright MCP available): load `report.html`,
  assert the tree SVG and the alignment canvas both render; gated so a
  browserless CI skips it.

## 11. Out of scope (YAGNI)

- Metadata upload / arbitrary trait coloring beyond fastbaps clusters.
- Multi-run comparison dashboards.
- Server-backed or live-updating reports.
- Radial/unrooted tree layouts (rectangular only for the alignment bind).
