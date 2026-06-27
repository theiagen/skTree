# Vendored front-end assets

| File | Source | Version | License |
|---|---|---|---|
| theiaphyloviewer.iife.js | TheiaPhyloViewer (SVG renderer), bundled via Vite | see `vendor_build/` | Proprietary — used with the owner's grant |

The report has a single third-party front-end dependency: **TheiaPhyloViewer**,
an SVG phylogenetic-tree viewer with always-on pan/zoom, per-tip styling and a
scale bar. The pure-ESM source (zero runtime dependencies) is bundled into one
self-contained IIFE that exposes the `TPV` global; `app.js` instantiates
`new TPV.TheiaPhyloViewer(host, { source: <newick>, ... })`.

TheiaPhyloViewer is proprietary; it is vendored here under an explicit grant
from its owner for use in skTree. Do not redistribute the source separately.

Everything else is hand-written in `app.js`: the **zoomable SNP alignment** is a
virtualized `<canvas>` (only on-screen columns are painted, so 20k-column
matrices stay smooth), and the cluster colour strip and legend use a small
built-in categorical palette. There is no other library to vendor.

The report is self-contained: every asset is inlined into a single HTML file by
`html_report.render_html` (no CDN, no network at view time), so it survives
air-gapped clusters and link rot.

## Rebuilding the bundle

See [`vendor_build/`](vendor_build/). With the TheiaPhyloViewer source checked
out and Vite available:

```bash
cd vendor_build
# point TPV_SRC at the TheiaPhyloViewer source's src/index.js, then:
npx vite build --config vite.config.mjs
cp dist/tpv.iife.js ../theiaphyloviewer.iife.js
```

Vite tree-shakes the renderer down to ~60 KB (the editing/serialisation paths
the report never calls are dropped).

## History

Earlier revisions drew the tree with d3 (a hand-written Newick parser feeding a
rectangular phylogram) and, before that, vendored `phylotree.min.js` +
`underscore.min.js`. Both were dropped: d3 in favour of TheiaPhyloViewer (the
owner's own viewer, with richer interaction and styling), removing the 280 KB
d3 bundle entirely since the alignment canvas needs no library.
