# Rebuilding `theiaphyloviewer.iife.js`

This directory reproduces the vendored TheiaPhyloViewer bundle used by the
report. The committed `../theiaphyloviewer.iife.js` is the output; regenerate it
only when updating the viewer.

## Steps

1. Obtain the TheiaPhyloViewer SVG-renderer source (pure ESM, zero runtime
   dependencies). Note its `src/index.js` path.
2. Build with Vite (Node 20+):

   ```bash
   export TPV_SRC=/abs/path/to/theiaphyloviewer/src/index.js
   npx vite build --config vite.config.mjs
   ```

   Or place the source under `tpv-src/` (so `tpv-src/src/index.js` exists) and
   omit `TPV_SRC`.
3. Copy the result over the vendored asset:

   ```bash
   cp dist/tpv.iife.js ../theiaphyloviewer.iife.js
   ```

## What this produces

- A single IIFE that assigns the global `TPV`, with `TPV.TheiaPhyloViewer` as
  the viewer constructor.
- Vite/Rollup tree-shakes the renderer to ~60 KB (gzip ~18 KB): only the
  reachable render path is kept; editing, session round-tripping and host
  adapters are dropped.

The bundle is then inlined verbatim into every report by
`html_report.render_html`, keeping the report offline and self-contained.
