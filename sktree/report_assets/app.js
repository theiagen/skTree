/* sktree/report_assets/app.js
 * Front-end for the self-contained skTree report.
 *
 * The phylogeny is drawn by TheiaPhyloViewer (vendored IIFE bundle, exposed as
 * the `TPV` global): an SVG tree with always-on pan/zoom, per-tip cluster
 * colours and a scale bar. Beside it sits an independently *zoomable* SNP
 * alignment rendered to a <canvas>: only the columns currently in view are
 * painted, so a 20k-column matrix stays smooth. The tree owns its own camera
 * (pan/zoom) but the three panels share one vertical grid: the alignment and the
 * cluster strip adopt the tree's rendered tip y-positions, so tip N, strip cell N
 * and alignment row N sit on one horizontal axis. The panels are linked both
 * ways too: hovering an alignment row highlights the matching tip via
 * `viewer.highlightNode(name)`, and clicking a tree tip rings that sample's
 * alignment row (and its cluster-strip cell). No d3, no network. */
(function () {
  "use strict";

  var TPV_CTOR = (window.TPV && window.TPV.TheiaPhyloViewer) || null;

  var BASE_COLORS = { A: "#3fa45b", C: "#2f72c4", G: "#e08a1e", T: "#d6453c" };
  var GAP_COLOR = "#e4e4e7";
  var MATCH_COLOR = "#f3f4f6";
  // A clicked tip rings its alignment row (and cluster-strip cell). An ink line
  // with a white outer fringe reads on any base colour or the pale match cell,
  // so the ring never loses contrast against the fill underneath it.
  var SELECT_INK = "#111827";
  var SELECT_HALO = "#ffffff";
  // Tableau-10 then a short Set3-ish tail, so up to a dozen clusters stay
  // distinguishable without pulling in d3's colour schemes.
  var CLUSTER_HEX = [
    "#4e79a7", "#f28e2b", "#e15759", "#76b7b2", "#59a14f", "#edc948",
    "#b07aa1", "#ff9da7", "#9c755f", "#bab0ac", "#8cd17d", "#b6992d",
  ];

  function el(id) { return document.getElementById(id); }
  function isBase(b) { return b === "A" || b === "C" || b === "G" || b === "T"; }
  function clamp(v, lo, hi) { return v < lo ? lo : v > hi ? hi : v; }

  // Fixed vertical chrome around the tree (top/bottom padding for the root edge
  // and the deepest leaf label) added on top of the rows' combined height.
  var TREE_CHROME = 60;

  // Vertical space the visualization can occupy: from the top of the #viz card
  // to the bottom of the viewport, less a small margin. Drives the per-row
  // height so the tree + alignment fill the screen instead of leaving a void.
  function availHeight() {
    var viz = el("viz");
    var top = viz ? viz.getBoundingClientRect().top : 120;
    return Math.max(260, window.innerHeight - top - 28);
  }

  // Per-leaf row height: spread the leaves across the available height, but keep
  // a floor (legible) and a cap (so a handful of samples don't become absurdly
  // tall alignment cells). Subtract the tree's fixed chrome (TREE_CHROME, the
  // same +60 added to the tree height) so the visualization fills the viewport
  // without overshooting it.
  function rowHeightFor(n) {
    if (!n) return 22;
    return clamp((availHeight() - TREE_CHROME) / n, 14, 96);
  }
  function baseColor(b) { return BASE_COLORS[b] || GAP_COLOR; }
  function clusterHex(id) { return CLUSTER_HEX[(id - 1) % CLUSTER_HEX.length]; }

  function hexToRgba(hex) {
    var h = hex.replace("#", "");
    return [
      parseInt(h.slice(0, 2), 16),
      parseInt(h.slice(2, 4), 16),
      parseInt(h.slice(4, 6), 16),
      255,
    ];
  }

  // -- Newick (leaf order only) ----------------------------------------------

  // Recursive-descent parser, used solely to read the leaf order so the
  // alignment rows line up top-to-bottom with the tree's tips (TheiaPhyloViewer
  // lays leaves out in source order, which is what this walk reproduces).
  function orderedLeafNames(text) {
    var s = text.replace(/\[[^\]]*\]/g, "").trim();
    if (s.charAt(s.length - 1) === ";") s = s.slice(0, -1);
    var pos = 0;
    var out = [];

    function node() {
      var isLeaf = true;
      if (s.charAt(pos) === "(") {
        isLeaf = false;
        pos++;
        for (;;) {
          node();
          var sep = s.charAt(pos++);
          if (sep === ")") break;
          if (sep !== ",") throw new Error("malformed Newick near " + pos);
        }
      }
      var name = label();
      if (s.charAt(pos) === ":") { pos++; skipNumber(); }
      if (isLeaf) out.push(name);
    }
    function label() {
      if (s.charAt(pos) === "'") {
        pos++;
        var start = pos;
        while (pos < s.length && s.charAt(pos) !== "'") pos++;
        var q = s.slice(start, pos);
        pos++;
        return q;
      }
      var begin = pos;
      while (pos < s.length && ":,()".indexOf(s.charAt(pos)) === -1) pos++;
      return s.slice(begin, pos).trim();
    }
    function skipNumber() {
      while (pos < s.length && ":,()".indexOf(s.charAt(pos)) === -1) pos++;
    }
    node();
    return out;
  }

  function consensus(block, order) {
    var cols = block.columns;
    var cons = new Array(cols);
    for (var j = 0; j < cols; j++) {
      var counts = {};
      for (var r = 0; r < order.length; r++) {
        var seq = block.rows[order[r]];
        if (!seq) continue;
        var b = seq.charAt(j);
        if (isBase(b)) counts[b] = (counts[b] || 0) + 1;
      }
      var best = null, bestN = -1;
      for (var k in counts) { if (counts[k] > bestN) { bestN = counts[k]; best = k; } }
      cons[j] = best;
    }
    return cons;
  }

  // -- tree (TheiaPhyloViewer) -----------------------------------------------

  var viewer = null;

  function clusterStyles(clusters) {
    var styles = {};
    if (!clusters) return styles;
    Object.keys(clusters).forEach(function (name) {
      styles[name] = { fillColour: hexToRgba(clusterHex(clusters[name])) };
    });
    return styles;
  }

  function renderTree(newick, clusters, heightPx) {
    var host = el("tree");
    var size = { width: host.clientWidth || 360, height: heightPx };
    var styles = clusterStyles(clusters);
    if (!TPV_CTOR) {
      host.innerHTML = "<p class='warn'>Tree viewer failed to load.</p>";
      return null;
    }
    if (viewer) {
      viewer.setProps({ source: newick, styles: styles, size: size });
    } else {
      viewer = new TPV_CTOR(host, {
        source: newick,
        size: size,
        showLabels: true,
        showLeafLabels: true,
        showShapes: true,
        scalebar: true,
        interactive: { tooltip: true, highlight: true },
        styles: styles,
        fontSize: 12,
        nodeSize: 9,
        padding: 14,
      });
      installSelectionBridge();
    }
    return viewer;
  }

  function highlightTip(name) {
    if (viewer && typeof viewer.highlightNode === "function") {
      try { viewer.highlightNode(name); } catch (e) { /* viewer not ready */ }
    }
  }

  // Mirror the tree's persistent selection (clicked tips) onto the alignment.
  // TheiaPhyloViewer exposes no click callback, but its internal `_onClick`
  // fires exactly when `props.selectedIds` changes (leaf ids === sample names),
  // so we wrap it once: run the original, then copy the selection across. This
  // is synchronous with the selection change, so no rAF/timing guesswork. A
  // leaf id maps straight to an alignment row because both share the Newick
  // leaf order. Drag-panning never reaches `_onClick` (TPV gates it on a
  // not-dragged flag), so panning the tree won't spuriously clear the row.
  function syncSelectionFromTree() {
    var ids = (viewer && viewer.props && viewer.props.selectedIds) || [];
    var next = Object.create(null);
    for (var i = 0; i < ids.length; i++) next[ids[i]] = true;
    aln.selected = next;
    alnDraw();
    drawStrip();
  }

  function installSelectionBridge() {
    if (!viewer || viewer.__skSelBridge || typeof viewer._onClick !== "function") return;
    viewer.__skSelBridge = true;
    var orig = viewer._onClick.bind(viewer);
    viewer._onClick = function () {
      orig.apply(null, arguments);
      syncSelectionFromTree();
    };
  }

  // -- alignment (virtualized, zoomable canvas) ------------------------------

  var aln = {
    block: null, order: [], cons: null, diff: false, scope: "core",
    rowHeight: 12, colW: 6, minColW: 0.5, maxColW: 26, scrollX: 0,
    cssW: 0, rowsH: 0, dpr: 1, ctx: null, raf: null, dragging: false,
    dragStartX: 0, dragStartScroll: 0, ready: false, showBases: false,
    selected: Object.create(null),
    // Vertical grid is adopted from the tree's rendered tip positions
    // (applyTreeRowGeometry): rowTop0 = top of row 0, rowHeight = tip pitch.
    rowTop0: 0,
  };

  // Cells must be at least this wide/tall (CSS px) before a base letter fits
  // legibly; below it the toggle stays on but letters are suppressed until the
  // user zooms in. The label colour contrasts with the cell underneath.
  var LETTER_MIN_COLW = 7;
  var LETTER_MIN_ROWH = 9;
  var LETTER_DARK = "#1f2430";
  var LETTER_LIGHT = "#ffffff";

  function alnMaxScroll() {
    return Math.max(0, aln.block.columns * aln.colW - aln.cssW);
  }

  function alnSizeCanvas() {
    var canvas = el("aln");
    aln.dpr = window.devicePixelRatio || 1;
    aln.cssW = el("aln-wrap").clientWidth || 600;
    canvas.style.width = aln.cssW + "px";
    canvas.style.height = aln.rowsH + "px";
    canvas.width = Math.max(1, Math.floor(aln.cssW * aln.dpr));
    canvas.height = Math.max(1, Math.floor(aln.rowsH * aln.dpr));
  }

  function alnDrawNow() {
    aln.raf = null;
    if (!aln.block) return;
    var ctx = aln.ctx, colW = aln.colW, cols = aln.block.columns;
    ctx.setTransform(aln.dpr, 0, 0, aln.dpr, 0, 0);
    ctx.clearRect(0, 0, aln.cssW, aln.rowsH);
    var first = Math.max(0, Math.floor(aln.scrollX / colW));
    var last = Math.min(cols - 1, Math.ceil((aln.scrollX + aln.cssW) / colW));
    var w = Math.max(1, Math.ceil(colW));
    var h = Math.max(1, Math.ceil(aln.rowHeight));
    // Letters only render once cells are big enough to hold them; the toggle can
    // be left on while zoomed out (cheaper than re-querying on every zoom step).
    var letters = aln.showBases && colW >= LETTER_MIN_COLW &&
      aln.rowHeight >= LETTER_MIN_ROWH;
    if (letters) {
      ctx.font = Math.max(7, Math.floor(Math.min(aln.rowHeight, colW) * 0.82)) +
        "px ui-monospace, SFMono-Regular, Menlo, Consolas, monospace";
      ctx.textAlign = "center";
      ctx.textBaseline = "middle";
    }
    for (var r = 0; r < aln.order.length; r++) {
      var seq = aln.block.rows[aln.order[r]];
      if (seq === undefined) continue;
      var y = Math.round(aln.rowTop0 + r * aln.rowHeight);
      for (var j = first; j <= last; j++) {
        var b = seq.charAt(j);
        var fill;
        if (aln.diff) {
          fill = !isBase(b) ? GAP_COLOR : (b === aln.cons[j] ? MATCH_COLOR : baseColor(b));
        } else {
          fill = baseColor(b);
        }
        var x = Math.floor(j * colW - aln.scrollX);
        ctx.fillStyle = fill;
        ctx.fillRect(x, y, w, h);
        if (letters && isBase(b)) {
          // Dark glyph on the pale match cell, white on the saturated bases.
          ctx.fillStyle = fill === MATCH_COLOR ? LETTER_DARK : LETTER_LIGHT;
          ctx.fillText(b, x + colW / 2, y + aln.rowHeight / 2);
        }
      }
    }
    drawSelectedRows(ctx, h);
  }

  // Ring every row whose sample is currently selected in the tree. Drawn last so
  // it sits above the painted cells. A 4px white stroke laid under a 2px ink
  // stroke on the same rectangle leaves a 1px white fringe outside the ink, so
  // the ring stays visible whether the row is mostly saturated bases or pale
  // match cells. Rows always fit the viewport (no vertical scroll), so a
  // selected row is on screen the instant its tip is clicked.
  function drawSelectedRows(ctx, h) {
    var sel = aln.selected;
    var any = false;
    for (var key in sel) { any = true; break; }
    if (!any) return;
    var w = Math.max(0, aln.cssW - 3);
    for (var r = 0; r < aln.order.length; r++) {
      if (!sel[aln.order[r]]) continue;
      var y = Math.round(aln.rowTop0 + r * aln.rowHeight);
      var rh = Math.max(1, h) - 3;
      if (rh <= 0) continue;
      ctx.lineWidth = 4;
      ctx.strokeStyle = SELECT_HALO;
      ctx.strokeRect(1.5, y + 1.5, w, rh);
      ctx.lineWidth = 2;
      ctx.strokeStyle = SELECT_INK;
      ctx.strokeRect(1.5, y + 1.5, w, rh);
    }
  }

  function alnDraw() {
    if (aln.raf == null) aln.raf = requestAnimationFrame(alnDrawNow);
  }

  function alnClampScroll() {
    aln.scrollX = clamp(aln.scrollX, 0, alnMaxScroll());
  }

  function alnZoomAt(clientPx, factor) {
    var worldCol = (clientPx + aln.scrollX) / aln.colW;
    aln.colW = clamp(aln.colW * factor, aln.minColW, aln.maxColW);
    aln.scrollX = worldCol * aln.colW - clientPx;
    alnClampScroll();
    alnDraw();
  }

  function drawStrip() {
    var strip = el("strip");
    strip.innerHTML = "";
    if (!aln.clusters) return;
    var canvas = document.createElement("canvas");
    var dpr = window.devicePixelRatio || 1;
    canvas.style.width = "14px";
    canvas.style.height = aln.rowsH + "px";
    canvas.width = Math.floor(14 * dpr);
    canvas.height = Math.max(1, Math.floor(aln.rowsH * dpr));
    var ctx = canvas.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    var h = Math.max(1, Math.ceil(aln.rowHeight));
    for (var r = 0; r < aln.order.length; r++) {
      var c = aln.clusters[aln.order[r]];
      if (c === undefined) continue;
      ctx.fillStyle = clusterHex(c);
      ctx.fillRect(0, Math.round(aln.rowTop0 + r * aln.rowHeight), 14, h);
    }
    // Same ring as the alignment row, so a clicked tip lights up tree → strip →
    // alignment as one connected band across the three panels.
    for (var s = 0; s < aln.order.length; s++) {
      if (!aln.selected[aln.order[s]]) continue;
      var sy = Math.round(aln.rowTop0 + s * aln.rowHeight);
      var sh = h - 3;
      if (sh <= 0) continue;
      ctx.lineWidth = 3;
      ctx.strokeStyle = SELECT_HALO;
      ctx.strokeRect(1.5, sy + 1.5, 11, sh);
      ctx.lineWidth = 1.5;
      ctx.strokeStyle = SELECT_INK;
      ctx.strokeRect(1.5, sy + 1.5, 11, sh);
    }
    strip.appendChild(canvas);
  }

  function alnSetData(data, scope, order, diff) {
    var prevBlock = aln.block;
    var prevScope = aln.scope;
    aln.block = data.alignment[scope];
    aln.clusters = data.clusters;
    aln.order = order;
    aln.diff = diff;
    aln.scope = scope;
    aln.cons = diff ? consensus(aln.block, order) : null;
    // Horizontal fit only; the vertical grid is adopted from the tree by
    // applyTreeRowGeometry(). Re-fit the column width and rewind to the start
    // only when the column set actually changes -- first load or a core<->all
    // scope switch. Toggling "differences only" recolours the same block, so the
    // user's current zoom and scroll position are kept intact.
    aln.cssW = el("aln-wrap").clientWidth || 600;
    if (!prevBlock || prevScope !== scope) {
      var fit = aln.cssW / aln.block.columns;
      aln.colW = clamp(fit, aln.minColW, 8);
      aln.scrollX = 0;
    } else {
      alnClampScroll();
    }
  }

  // Adopt the tree's rendered tip positions as the shared vertical grid for the
  // strip and alignment, so tip N, strip cell N and alignment row N line up on
  // one horizontal axis. TheiaPhyloViewer owns leaf layout — its padding, scale
  // bar and label metrics set the exact pitch — so rather than guess we read the
  // tip-marker centres straight from the SVG. #tree, #strip and #aln-wrap are
  // flex siblings pinned to the same top, so a y measured against the tree host
  // applies unchanged inside the canvases. Falls back to uniform rows from the
  // top if the markup ever stops giving exactly one circle per leaf.
  function applyTreeRowGeometry() {
    var n = aln.order.length;
    var host = el("tree");
    var svg = host ? host.querySelector("svg") : null;
    var circles = svg ? svg.querySelectorAll("circle") : [];
    if (svg && circles.length === n && n > 0) {
      var hostTop = host.getBoundingClientRect().top;
      var ys = [];
      for (var i = 0; i < circles.length; i++) {
        var rc = circles[i].getBoundingClientRect();
        ys.push(rc.top + rc.height / 2 - hostTop);
      }
      ys.sort(function (a, b) { return a - b; });
      var pitch = n > 1 ? (ys[n - 1] - ys[0]) / (n - 1) : host.clientHeight;
      aln.rowHeight = pitch;
      aln.rowTop0 = ys[0] - pitch / 2;
      aln.rowsH = host.clientHeight;
    } else {
      aln.rowHeight = rowHeightFor(n);
      aln.rowTop0 = 0;
      aln.rowsH = aln.rowHeight * n;
    }
    alnSizeCanvas();
    if (aln.block) alnClampScroll();
    alnDraw();
  }

  function alnColAt(clientX) {
    var rect = el("aln").getBoundingClientRect();
    return Math.floor((clientX - rect.left + aln.scrollX) / aln.colW);
  }
  function alnRowAt(clientY) {
    var rect = el("aln").getBoundingClientRect();
    return Math.floor((clientY - rect.top - aln.rowTop0) / aln.rowHeight);
  }

  function initAlnHandlers() {
    if (aln.ready) return;
    aln.ready = true;
    var canvas = el("aln");
    aln.ctx = canvas.getContext("2d");
    var tip = el("tooltip");

    canvas.addEventListener("wheel", function (ev) {
      ev.preventDefault();
      var rect = canvas.getBoundingClientRect();
      var factor = ev.deltaY < 0 ? 1.18 : 1 / 1.18;
      alnZoomAt(ev.clientX - rect.left, factor);
    }, { passive: false });

    canvas.addEventListener("mousedown", function (ev) {
      aln.dragging = true;
      aln.dragStartX = ev.clientX;
      aln.dragStartScroll = aln.scrollX;
      canvas.classList.add("dragging");
      ev.preventDefault();
    });
    window.addEventListener("mousemove", function (ev) {
      if (!aln.dragging) return;
      aln.scrollX = aln.dragStartScroll - (ev.clientX - aln.dragStartX);
      alnClampScroll();
      alnDraw();
      tip.hidden = true;
    });
    window.addEventListener("mouseup", function () {
      if (!aln.dragging) return;
      aln.dragging = false;
      canvas.classList.remove("dragging");
    });

    canvas.addEventListener("mousemove", function (ev) {
      if (aln.dragging || !aln.block) return;
      var col = alnColAt(ev.clientX);
      var row = alnRowAt(ev.clientY);
      var name = aln.order[row];
      if (col < 0 || col >= aln.block.columns || name === undefined) {
        tip.hidden = true;
        return;
      }
      highlightTip(name);
      tip.hidden = false;
      tip.style.left = (ev.clientX + 12) + "px";
      tip.style.top = (ev.clientY + 12) + "px";
      tip.innerHTML = "<b>" + name + "</b><br>" + aln.scope + " SNP #" +
        (col + 1) + " / " + aln.block.columns;
    });
    canvas.addEventListener("mouseleave", function () {
      tip.hidden = true;
      highlightTip(null);
    });

    var basesBtn = el("bases-toggle");
    basesBtn.onclick = function () {
      aln.showBases = !aln.showBases;
      basesBtn.classList.toggle("active", aln.showBases);
      basesBtn.setAttribute("aria-pressed", aln.showBases ? "true" : "false");
      alnDraw();
    };

    el("zoom-in").onclick = function () { alnZoomAt(aln.cssW / 2, 1.4); };
    el("zoom-out").onclick = function () { alnZoomAt(aln.cssW / 2, 1 / 1.4); };
    el("zoom-fit").onclick = function () {
      aln.colW = clamp(aln.cssW / aln.block.columns, aln.minColW, 8);
      aln.scrollX = 0;
      alnDraw();
    };

    var resizeRaf = null;
    window.addEventListener("resize", function () {
      if (resizeRaf != null) return;
      resizeRaf = requestAnimationFrame(function () {
        resizeRaf = null;
        if (!aln.block) return;
        // Re-fit the tree to the new viewport first, then re-adopt its tip grid
        // so the strip and alignment track it. Horizontal zoom (colW) is
        // independent, so the user's alignment zoom is preserved across resize.
        var n = aln.order.length;
        if (viewer) {
          viewer.setProps({ size: { width: el("tree").clientWidth || 360, height: rowHeightFor(n) * n + TREE_CHROME } });
        }
        applyTreeRowGeometry();
        drawStrip();
        syncSidebarHeight();
      });
    });
  }

  // -- panels ----------------------------------------------------------------

  function row(k, v) { return "<tr><td>" + k + "</td><td>" + v + "</td></tr>"; }

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
    legend += " <span class='swatch' style='background:" + MATCH_COLOR +
      "'></span>matches consensus";
    el("legend").innerHTML = "<h2>Bases</h2><div class='chips'>" + legend + "</div>";

    if (data.clusters) {
      var ids = {};
      Object.keys(data.clusters).forEach(function (n) { ids[data.clusters[n]] = true; });
      var swatches = Object.keys(ids).map(Number).sort(function (a, b) { return a - b; })
        .map(function (id) {
          return "<span class='swatch' style='background:" + clusterHex(id) +
            "'></span>" + id;
        }).join(" ");
      el("legend").innerHTML += "<h2>Clusters</h2><div class='chips'>" + swatches + "</div>";
    }

    renderFiles(data.outputs);
  }

  function esc(s) {
    return String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
  }

  function renderFiles(outputs) {
    var host = el("files-panel");
    if (!outputs || !outputs.length) { host.innerHTML = ""; return; }
    var items = outputs.map(function (f) {
      return "<li><code class='fname'>" + esc(f.name) + "</code>" +
        "<span class='fdesc'>" + esc(f.description) + "</span></li>";
    }).join("");
    host.innerHTML = "<h2>Output files</h2>" +
      "<p class='files-note'>Written alongside this report:</p>" +
      "<ul class='files'>" + items + "</ul>";
  }

  function truncationNote(data, scope) {
    var note = el("trunc-note");
    var trunc = data.alignment.all && data.alignment.all.truncatedTo;
    if (scope === "all" && trunc != null) {
      note.textContent = "Showing " + trunc +
        " most-informative of the full variable-SNP set (truncated).";
    } else {
      note.textContent = "";
    }
  }

  var current = { data: null, state: null };

  // Make the side panel exactly as tall as the visualization card so the two
  // columns share one bottom edge. The tree's tip pitch is width-limited (its
  // column is split with the alignment), so the visualization can't usefully
  // grow taller -- instead the sidebar is pinned to the visualization's height
  // and the Output files card (last in the column, flex: 1) absorbs the slack,
  // scrolling internally when the output list is long. Only on the wide
  // side-by-side layout; in the stacked layout the panel sits below the tree at
  // its natural height, so the forced height is cleared.
  function syncSidebarHeight() {
    var panels = el("panels");
    if (getComputedStyle(panels.parentNode).flexDirection !== "row") {
      panels.style.height = "";
      return;
    }
    panels.style.height = el("viz").offsetHeight + "px";
  }

  function fullRender(data, state) {
    current.data = data;
    current.state = state;
    var n = data.samples.length;
    var rowHeight = rowHeightFor(n);
    var order = orderedLeafNames(data.trees[state.tree]);
    renderTree(data.trees[state.tree], data.clusters, rowHeight * n + TREE_CHROME);
    alnSetData(data, state.aln, order, state.diff);
    // Adopt the freshly-rendered tree's tip grid so the strip and alignment rows
    // line up with the tips, then draw the strip on that same grid.
    applyTreeRowGeometry();
    drawStrip();
    syncSidebarHeight();
    truncationNote(data, state.aln);
    // The tree's selectedIds survive a re-render; re-derive the alignment rows
    // from it so a previously clicked tip stays ringed after switching trees or
    // toggling the diff/all-SNP views.
    syncSelectionFromTree();
  }

  window.SKTREE_RENDER = function (data) {
    var state = { tree: Object.keys(data.trees)[0], aln: "core", diff: false };
    initAlnHandlers();

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
    toggle.onclick = function () {
      var on = toggle.getAttribute("aria-pressed") !== "true";
      toggle.setAttribute("aria-pressed", on ? "true" : "false");
      state.aln = on ? "all" : "core";
      fullRender(data, state);
    };

    var diff = el("diff-toggle");
    diff.onclick = function () {
      var on = diff.getAttribute("aria-pressed") !== "true";
      diff.setAttribute("aria-pressed", on ? "true" : "false");
      state.diff = on;
      fullRender(data, state);
    };

    panels(data);
    fullRender(data, state);
    window.__SKTREE_READY = true;
  };
})();
