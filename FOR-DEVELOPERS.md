# skTree — A Developer's Guide

> The whole project, in plain language: what it does, how it's built, why it's
> built that way, and the potholes we hit so you don't have to.

---

## 1. What problem are we solving?

Imagine you run a hospital lab. Last week three patients turned up with the same
nasty *Klebsiella*. Are these three independent infections, or is one ward
quietly seeding the other two? To answer that you need to know how genetically
close the bugs are — ideally a tree showing who's cousin to whom.

The classic way to build that tree is: pick a reference genome, align every
isolate to it, read off the SNPs (single-base differences), and feed those SNPs
to a tree-builder. The trouble is the **reference**. Bacteria have wildly plastic
genomes — genes come and go — so a single reference is a straitjacket. Pick the
wrong one and you throw away exactly the regions that distinguish your isolates.

The **split-k-mer idea** sidesteps this entirely: forget references and
alignments. Chop every genome into overlapping odd-length *k*-mers. Two
*k*-mers that are identical except for their **middle base** are pointing at the
same genomic locus with a SNP sitting in the centre. Collect all such loci across
your samples and you've discovered SNPs *reference-free* and *alignment-free*.

**SKA / SKA2** took that split-*k*-mer insight and made it scream-fast in
Rust. A split *k*-mer is an odd *k*-mer cut into two equal flanks around a
variable middle base — `XXXX·XXXX`. SKA stores a hash map from *flank* → *middle
base*; comparing samples is then just comparing what middle base each one parks
in the same flank slot.

**skTree** is a thin wrapper around SKA2: it drives SKA2 as the fast engine for
the heavy *k*-mer work, then turns its raw output into the *deliverables* an
end-to-end workflow needs — optimal-*k* selection, core/majority SNP
partitioning, a SNP matrix, and a phylogenetic tree — in approachable Python.

```
        bacterial assemblies / reads
                    │
                    ▼
        ┌───────────────────────┐
        │   SKA2 (Rust binary)  │   build → align → distance
        │   split-k-mer engine  │
        └───────────┬───────────┘
                    │ SNP alignment (FASTA), distances (TSV)
                    ▼
        ┌───────────────────────┐
        │   skTree (Python)     │   classify SNPs, build trees, report
        └───────────┬───────────┘
                    ▼
     alignment · snp_matrix · trees · summary
```

---

## 2. Technical architecture

### 2.1 The big decision: wrap, don't reimplement

SKA2 is a mature, benchmarked Rust tool with no Python bindings. We had three
options:

1. **Port** the split-*k*-mer algorithm to Python — slow and a maintenance trap.
2. **Bind** to the Rust library via FFI — brittle, needs a build toolchain.
3. **Wrap** the `ska` CLI as a subprocess — simple, robust, version-tolerant.

We chose **#3**. The entire coupling to SKA lives behind one class,
`SkaRunner` (`sktree/engine/ska.py`), and *every* subprocess call funnels through
a single private method, `_run()`. That gives us exactly one place to handle
errors, capture stdout, and (in tests) mock the boundary. If SKA changes a flag
tomorrow, there's one file to touch.

> **Design rule we kept:** the engine layer knows *nothing* about phylogenetics,
> and the science layer knows *nothing* about subprocesses. Each can be tested in
> isolation. This is the [hexagonal / ports-and-adapters](https://en.wikipedia.org/wiki/Hexagonal_architecture_(software))
> idea in miniature — SKA is an "adapter" plugged into a port.

### 2.2 The pipeline as a thin conductor

`run_pipeline()` (`sktree/pipeline.py`) is deliberately boring. It's a conductor,
not a performer: it sequences the steps and owns *no* algorithms itself.

```
inputs ──▶ build .skf ──▶ align (SNP FASTA) ──▶ distance (TSV)
                                                      │
                          parse + classify ◀──────────┘
                                  │
         ┌────────────────────────┼─────────────────────────┐
         ▼                        ▼                          ▼
   SNP matrix TSV          core SNP FASTA            NJ / MP / ML trees
                                  │
                                  ▼
                             summary.txt
```

Each box is a separate, independently-tested module. The conductor just wires
their inputs and outputs and adds cross-cutting concerns (logging, the divergence
warning).

### 2.3 Module map

| Module | Responsibility | Depends on |
|--------|----------------|------------|
| `inputs.py` | resolve CLI files + manifest into named `Sample`s (auto-pair reads) | (stdlib only) |
| `logconfig.py` | console + always-on DEBUG file logging into the output dir | (stdlib only) |
| `engine/ska.py` | subprocess wrapper for the `ska` binary | `inputs` (`Sample`) |
| `engine/ml.py` | optional IQ-TREE / RAxML-NG wrapper | (stdlib only) |
| `snps.py` | parse SNP FASTA, classify core/majority/variable loci | NumPy |
| `kselect.py` | Kchooser-style optimal odd-*k* selection | (stdlib only) |
| `phylo.py` | neighbor-joining tree + SNP distance matrix | DendroPy, NumPy |
| `parsimony.py` | pure-Python maximum-parsimony (Fitch + search) | DendroPy |
| `annotate.py` | reference-coordinate SNP → gene + codon effect | Biopython, pyrodigal* |
| `engine/cluster.py` | optional `fastbaps` population-structure wrapper | (stdlib only) |
| `report.py` | SNP matrix TSV + human summary | (none) |
| `html_report.py` | assemble + render the self-contained HTML report | NumPy, `report_assets/` |
| `pipeline.py` | orchestration (`sktree run`) | all of the above |
| `cli.py` | argparse front door | `pipeline` |

The dependency arrows only ever point *down* this table — no cycles. `cli` →
`pipeline` → {everything else} → {engine, libraries}. The `*` on pyrodigal marks
an **optional** dependency: `annotate.py` imports cleanly without it and only
reaches for pyrodigal inside `predict_genes`, so the core pipeline never pays for
a feature it isn't using.

---

## 3. The interesting bits (where the science meets the code)

### 3.1 Why *odd* k, and why it matters in code

A split *k*-mer needs two **equal** flanks around one central base, so
`k = flank + 1 + flank` must be odd. `SkaRunner.build()` raises `ValueError` if
you hand it an even *k*. This isn't pedantry — an even *k* has no single centre
base, so the whole "SNP in the middle" abstraction collapses. We enforce the
invariant at the boundary so it can never leak downstream.

### 3.2 Core vs majority vs pan SNPs

Not every SNP is seen in every sample — a *k*-mer can be absent because of a real
deletion, low coverage, or divergence. skTree sorts loci into:

- **core** — present in *all* samples. The gold standard for trees; no missing
  data to confuse distance/likelihood.
- **majority** — present in at least a threshold fraction (default 0.5).
- **variable / pan** — variable in at least two samples, however sparse.

In `snps.py` this is pure NumPy vectorisation. The SNP alignment is loaded into a
2-D array of single characters (`dtype='<U1'`). A locus is "present" for a sample
if its base is one of `A/C/G/T` (`np.isin`), and "variable" if it shows ≥2
distinct real alleles. `present_fraction` is then a column-wise mean, and the
core/majority masks are simple boolean comparisons. No Python loops over loci —
the work happens in C inside NumPy.

> **Insight:** representing the alignment as a typed character matrix instead of a
> list of `SeqRecord`s is what makes the classification a three-line vectorised
> operation. The data structure *is* the algorithm.

### 3.3 Optimal-k selection

`kselect.py` answers "what *k* should I use?" with a simple uniqueness scan:
the best *k* is the **smallest odd *k* where almost all *k*-mers are unique**.
Too small and unrelated loci collide (false SNPs); too large and a single error
splinters a locus into many.

The subtle bit is **canonicalisation**. DNA is double-stranded, so a *k*-mer and
its reverse-complement are the same physical thing. Before counting uniqueness we
map each *k*-mer to the lexicographically smaller of {itself, its RC} using a
`str.translate` complement table. Skip this and you'd double-count every *k*-mer
and badly misjudge uniqueness.

### 3.4 Three trees, three philosophies

skTree can build the same tree-shaped answer three ways:

- **Neighbor-joining** (`phylo.py`) — distance-based, fast, deterministic. We
  compute a pairwise SNP-distance matrix (with *pairwise deletion* of missing
  sites) and hand it to **DendroPy**'s `nj_tree()`. Great default.
- **Maximum parsimony** (`parsimony.py`) — "the simplest history wins." Scored
  with **Fitch's algorithm**, searched with stepwise addition + NNI hill-climbing.
  Written in pure Python so it has zero external-binary dependencies (see §3.5).
- **Maximum likelihood** (`engine/ml.py`) — the statistical heavyweight. Real ML
  search is a solved problem in **IQ-TREE** / **RAxML-NG**, so we *wrap* one if
  present and **skip gracefully** if not. We don't reinvent a likelihood engine.
  We also don't second-guess it: with no `-m`, IQ-TREE runs its own
  **ModelFinder Plus** default and picks the substitution model from the data.
  On the demo's core-SNP alignment it selects `TVM+F+ASC+R3` — note the `+ASC`
  ascertainment-bias correction it adds on its own, exactly what a
  variable-sites-only alignment needs and something a hardcoded `GTR+G` would
  have silently suppressed. The pipeline collects the whole `tree_ml.*` family
  the engine writes (tree, the human-readable `.iqtree` report, the run `.log`,
  plus checkpoint/distance side files), not just the bare Newick tree.

### 3.5 Fitch parsimony, done iteratively

The parsimony scorer is the most algorithm-heavy code in the project, so it's
worth understanding. Fitch's algorithm scores one site on a tree in a single
post-order sweep:

- At each leaf, the state set is `{its base}` (or the full `{A,C,G,T}` set if the
  base is missing — that way missing data never adds cost).
- At each internal node, intersect the children's state sets. If the
  intersection is **non-empty**, keep it. If it's **empty**, the lineages
  disagree → add **1** to the score and take the **union** instead.

The total over all sites is the parsimony length. Two implementation choices are
deliberate:

1. **Iterative, not recursive.** SNP alignments can have thousands of sites and
   genomes can make deep trees; Python's recursion limit and per-call overhead
   would bite. The post-order is done with an explicit stack.
2. **`frozenset` states.** Immutable, hashable, and `frozenset.intersection`/
   `union` read exactly like the algorithm's prose.

The tree *search* uses a compact integer-adjacency representation (`_Tree`) for
speed during the thousands of candidate-scoring calls, then exports a DendroPy
tree at the very end so it shares output plumbing with NJ. Leaves are integers
`0..n-1`, internal nodes `≥ n` — cheap to create, link, and unlink.

> **Honest limitation:** this pure-Python search is roughly `O(n³ · sites)`. It's
> perfect for outbreak-scale sets (tens of isolates) but will crawl past ~30
> samples — so the pipeline *warns* above that count rather than silently hanging.
> For big sets, use `--ml` with IQ-TREE instead.

### 3.6 Stepping out of the reference-free world (SNP annotation)

Everything above is deliberately *reference-free* — that's SKA's whole premise.
But "which gene does this SNP land in, and does it change the protein?" is a
question you can only answer against coordinates, so `--reference` opts into a
reference for that final labelling step. The flow (`annotate.py`):

1. **`ska map -f vcf`** projects the split-k-mer SNPs onto the chosen reference,
   giving us positions in *reference* coordinates as a standard VCF.
2. **pyrodigal** predicts genes on that same reference. We run it in
   **metagenomic mode** (`GeneFinder(meta=True)`) so there's no self-training
   pass — that matters because a reference can be a single short contig, which is
   too little data to train a species-specific model on.
3. For each SNP we find the containing gene (if any) and classify it
   **intergenic** or **coding**; coding SNPs get a **synonymous /
   non-synonymous** call by translating the affected codon before and after the
   substitution.

The codon logic has one genuinely fiddly part: **minus-strand genes**. We always
slice the codon out of the forward reference first, apply the substitution in
forward coordinates, *then* reverse-complement the codon before translating — so
the frame and the strand are handled in exactly one place. The codon index math
differs per strand: on `+` it's `(pos - gene.start)//3`, on `-` it's
`(gene.end - pos)//3` counting inward from the gene's 3′ end.

> **Design payoff:** the pure logic (`parse_vcf`, `annotate_variants`,
> `_codon_effect`) takes plain dataclasses, so the whole synonymous/non-synonymous
> story is unit-tested on tiny hand-built ORFs (`ATGAAATTTGGGTAA`) with **no ska
> and no pyrodigal** in the loop. Only `predict_genes` touches pyrodigal, and the
> tests `importorskip` it. Push the I/O to the edges, keep the core pure.

The same `ska map` call has a second mode — `-f aln` — that emits a
**reference-anchored pseudo-alignment**: one row per sample, every column a
reference coordinate. `--map-tree` turns that into trees by feeding it straight
back through the *same* SNP layer the reference-free path uses
(`parse_alignment` → `classify` → NJ/parsimony on the variable subset, ML on the
core subset). The map-tree path adds **zero new tree engines** — it is just a
second input source for the machinery we already have, which is why it cost a
wrapper flag and a handful of pipeline lines, not a new module. One pitfall to
keep in mind: `SnpMatrix.is_variable()` scans columns with `np.unique`, which is
free on a SNP-only alignment but `O(sites)` on a full reference-length
pseudo-alignment (~millions of columns). It hasn't bitten us at outbreak scale,
but if map-trees ever feel slow, that loop is the first place to look — vectorise
it before reaching for anything cleverer.

### 3.7 The HTML report — a message in a bottle

Picture a genomic-surveillance analyst on an air-gapped HPC cluster. They run
skTree, scp a single file to a locked-down review laptop, and double-click it.
No internet, no install, no "please enable JavaScript from this CDN." The tree
draws, the SNP alignment scrolls alongside it, the cluster colours light up. Open
that same file in 2030 and it still works. That is the whole design brief for
`report.html`: a **message in a bottle** — everything it will ever need sealed
inside, thrown into the future, readable wherever it washes up.

That ambition dictates almost every decision in `html_report.py` and the vendored
`report_assets/` folder:

- **Inline everything; trust no CDN.** `render_html` reads each asset off disk
  (`theiaphyloviewer.iife.js`, our own `app.js`, `styles.css`) and substitutes
  them into `template.html` via `string.Template.safe_substitute`. The JSON
  payload of trees, alignment, stats, and clusters is inlined the same way
  (`$DATA_JSON`). We use `safe_substitute` rather than `substitute` on purpose:
  minified JavaScript is full of `$` tokens, and a strict substitution would
  choke on every one of them. The output references **zero** external URLs — a
  property we pin with a test that greps the rendered document for `http`.

  > **Why this is the right paranoia:** a CDN `<script src>` is a dependency on
  > someone else's uptime, someone else's TLS cert, and someone else's decision
  > not to delete a version. Air-gapped clusters can't reach it *today*; link rot
  > guarantees the rest of the world can't reach it *eventually*. Inlining trades
  > a few hundred KB of disk for a report that is correct forever and offline.

- **The report names its own siblings.** A self-contained HTML file is great for
  sharing, but it hides the fact that a run also drops a `.fasta` alignment, a
  `.tsv` matrix, Newick trees, and so on next to it. So the payload carries an
  `outputs` list — `{name, description}` for every result file — rendered as an
  "Output files" panel. The list is built in the **pipeline**, not the report
  code, from the actual `Path` objects that were written (`path.name`), with
  optional steps (MP/ML trees, clusters, annotation) dropped when their path is
  `None`. That keeps the displayed names honest: the report can only ever list a
  file that genuinely exists on disk, even if a filename later changes. We don't
  hyperlink them — a report opened from `file://` or moved elsewhere would have
  dangling links — so plain names plus descriptions is the durable choice.

- **The phylogram is TheiaPhyloViewer — a real viewer, not a hand-rolled pen.**
  The tree is drawn by **TheiaPhyloViewer** (TPV), a pure-ESM, zero-dependency
  SVG phylo viewer we vendor under its owner's grant. We bundle its 36-file
  source with Vite in library mode (IIFE), which tree-shakes the reachable render
  path down to ~60 KB (gzip ~18 KB) and exposes a single global,
  `TPV.TheiaPhyloViewer`; `app.js` does `new TPV.TheiaPhyloViewer(host, {source:
  newick, styles, ...})`. TPV draws branches *to scale*, colours each tip from a
  `styles[leafId] = {fillColour: [r,g,b,a]}` map (our cluster palette, converted
  hex → rgba), and renders its own scale bar — so a near-clonal clade visibly
  hugs the root while a divergent isolate stretches far to the right. This
  replaced an earlier hand-written d3 phylogram; dropping d3 removed ~280 KB of
  vendored JavaScript outright, because the only other thing that ever needed a
  library — the alignment — is now a plain `<canvas>` (see below).

  > **Lesson — the two-tier label switch.** TPV gates *all* text behind a master
  > `showLabels` boolean, with `showLeafLabels` / `showBlockLabels` /
  > `showInternalLabels` as sub-switches beneath it; **both** the master and the
  > specific switch must be true, and both default to `false`. Passing only
  > `showLeafLabels: true` rendered a tree with shapes but no tip names — silent,
  > because each switch is individually valid. The fix is one line (`showLabels:
  > true`), but the diagnosis took reading the renderer: when a library's flag
  > has no effect, suspect a gating flag above it before suspecting your data.

- **Always-on pan/zoom forced a decoupled architecture — then a *snapshot* re-coupled
  the rows.** TPV's camera pans and zooms regardless of any `interactive` prop — a
  deliberate design choice in the viewer. That ruled out *continuously* binding the
  alignment rows to the tree's live y-coordinates: the moment the user zooms the
  tree, a shared coordinate space would smear. So the panels stay independent
  cameras, linked semantically by **leaf order**: `orderedLeafNames` parses the
  (ladderized) Newick into the source-order list of tip names, and both the cluster
  strip and the alignment paint their rows in exactly that order. Hovering an
  alignment row calls the viewer's public `highlightNode(name)`, so the views talk
  through a stable name, not a fragile coordinate.

  On top of that semantic link sits a thin **vertical-registration** layer so the
  three panels read as one horizontal grid — tip *N*, its cluster cell, and its SNP
  row share a baseline. The trick is to *snapshot, don't bind*: after every tree
  render `applyTreeRowGeometry()` measures the rendered tip-marker circles
  (`getBoundingClientRect` on the SVG `<circle>`s, relative to the tree host),
  derives the exact row pitch and the top of row 0, and lays the strip + alignment
  rows on that grid. Because it samples once per render at the default camera — not
  on every frame of a live zoom — it gets pixel-precise registration without ever
  reintroducing the smear. The tree owns the layout (its padding, label metrics and
  scale bar set a pitch we could never guess); the other two panels adopt it. The
  alignment's zoom toolbar is `position: absolute` so the canvas can start at the
  column's very top and its first row lines up with the first tip; the toolbar
  floats in the thin headroom band above row 0 with a top-down fade, and is
  click-through except for its buttons so a drag begun in the band still pans.

- **Reading a selection out of a viewer that offers no callback.** The reverse
  link — click a tree tip, ring its alignment row — looked impossible at first:
  TPV exposes `highlightNode` but *no* click callback and *no* selection getter.
  The way through was to notice that a leaf's `id` **is** its sample name (the
  normalize step sets `id` from `name`), and that the internal `_onClick`
  mutates `props.selectedIds` to exactly those names on every click. So we wrap
  the instance's `_onClick` once: run the original, then copy `selectedIds` into
  the alignment's `selected` set and redraw. Wrapping the method (not listening
  on the container) means the sync is *synchronous with the selection change* —
  no `requestAnimationFrame` guesswork — and TPV already gates the click on a
  not-dragged flag, so panning the tree never clears the row. The ring itself is
  a 4px white stroke under a 2px ink stroke on the same rect: the white fringe
  keeps it legible whether the row is saturated bases or pale match cells.

- **Ladderize, or the tree rotates under you.** Neighbor-joining tie-breaking can
  swap two sibling clades run-to-run — same topology, same branch lengths, but a
  different top-to-bottom leaf order — because `ska align` emits SNP columns in a
  thread-dependent order and the resulting ULP-level wobble flips a tie. For an
  *unrooted* NJ tree that rotation is scientifically meaningless, but it makes the
  rendered report non-reproducible and harder to read. `to_newick` therefore calls
  `tree.ladderize(ascending=True)` before serialising, sorting each node's children
  by descendant count. The same data now always renders the same canonical,
  ladder-shaped tree — the convention every phylogenetics tool follows.

- **The zoomable alignment is a *virtualized* canvas.** The SNP alignment is its
  own `<canvas>`, and it zooms horizontally — scroll to zoom, drag to pan, plus
  `−`/`fit`/`+` buttons — so an analyst can go from "whole matrix at a glance" to
  "this one polymorphic column" without a separate viewer. The trick that keeps
  it smooth on a 20 000-column matrix is **virtualization**: `alnDrawNow` paints
  only the columns currently visible in the viewport (first..last under the
  current `colW` and `scrollX`), never the whole matrix. Zooming changes `colW`
  (the per-column pixel width) and repaints; the canvas *element* stays viewport-
  sized, so cost scales with what's on screen, not with the alignment length.
  `alnZoomAt` keeps the column under the cursor fixed while zooming — the small
  detail that makes zoom feel like a map instead of a slider.

- **Difference mode: let the eye find the signal.** A 32 000-column core alignment
  rendered in full ACGT colour is a beautiful, useless wall of confetti. The
  `highlight differences only` toggle computes a per-column **consensus** (most
  common unambiguous base) and mutes every cell that matches it to light grey,
  leaving only the *differences* in base colour. Suddenly a clonal clade collapses
  to near-blank rows while a divergent isolate lights up — the visualisation now
  shows you population structure instead of hiding it in noise.

- **`min-width: 0`, the one-line CSS that unbroke the page.** Even virtualized,
  the canvas lives inside a flex row alongside the tree and the stats panel, and
  flex items default to `min-width: auto` ("never shrink below my content"). That
  let the alignment column refuse to shrink and pushed its neighbours around.
  Setting `min-width: 0` on the alignment wrapper and its flex parent lets the box
  shrink to the viewport, so the canvas re-fits instead of overflowing — a classic
  flexbox trap, and a reminder that "let it shrink" must be granted explicitly.

- **Narrow viewports stack instead of squeeze.** Three columns competing in one
  flex row (tree · cluster strip · alignment) plus a side panel is fine at desktop
  width and miserable on a laptop split-screen: the alignment got crushed to a
  sliver. Two complementary fixes — the tree and alignment now *share* the row
  proportionally (`flex: 1 1 0` each, with a `min-width` floor on the tree so it
  never vanishes), and a `@media (max-width: 900px)` breakpoint drops the stats
  panel **below** the visualization so the alignment reclaims the full width. The
  reliable responsive primitive is still the oldest one: when columns won't fit
  side by side, stack them.

- **A hard cap on width, by informativeness.** A core alignment is small, but the
  *all-SNP* matrix can be enormous — and a browser will not happily paint a
  million-column canvas. So `build_report_data` caps the "all" view at
  `MAX_ALL_SNP_COLUMNS = 20_000`, and crucially it doesn't just truncate from the
  left: `select_informative_columns` scores every column by **minor-allele
  presence** (how many samples carry something other than the column's majority
  base) and keeps the most informative 20 000, re-sorted into genomic order. The
  pipeline logs that the view was capped, so a truncated report never lies about
  being complete. Monomorphic and near-monomorphic columns — the ones that tell
  you nothing — are exactly the ones dropped first.

- **fastbaps: fail loud or skip clean, never half-do.** Population-structure
  clusters come from `fastbaps`, an external tool we *wrap* rather than rebuild —
  the same call we made for IQ-TREE (§3.4) and SKA itself. `engine/cluster.py`
  resolves the binary with `shutil.which` at construction and raises
  `FastbapsNotAvailable` if it's missing, so the pipeline decides **once** whether
  clustering is possible. `--cluster` writes a standalone `clusters.csv`; `--html`
  opportunistically runs fastbaps too (clusters colour the tree) but degrades
  gracefully with a logged "skipping clusters: …" if it isn't installed. The
  report is still perfectly useful without colours — it just doesn't pretend it
  had data it didn't.

- **Clean names, fixed at the source.** SKA names every sample by the *verbatim
  path* you hand it, so `ska build /data/SAMD00052601.fa.gz …` would stamp that
  whole path onto every tip label, alignment row, matrix column, and cluster id —
  ugly and unwieldy. Rather than patch the labels at display time in three
  different renderers, we fix it once, upstream in `inputs.py`: `resolve_inputs`
  turns the raw CLI files (and any `--manifest`) into named `Sample`s, where
  `clean_sample_name` has already stripped the directory and the trailing
  compression/FASTA suffixes (`SAMD00052601.fa.gz → SAMD00052601`), de-duplicating
  any collisions. `engine/ska.py` then writes the tab-separated **file list**
  (`name<TAB>file…`) straight from those samples and passes it with `ska build -f`.
  Every downstream artifact inherits the clean stem for free, because the name was
  right before any of them ever saw it.

> **The throughline:** a report that has to survive offline machines and a decade
> of bit-rot is a *packaging* problem as much as a rendering one. Inline the
> assets, own your rendering, bound the payload, fix names at the source, and make
> every optional input fail loud or skip clean. Do that and the bottle floats.

### 3.8 Reads and assemblies, side by side (`inputs.py`)

SKA2 has always called SNPs straight from FASTQ reads — `ska build` takes a file
list where a line can be either `name<TAB>assembly` or `name<TAB>fwd<TAB>rev` for
paired reads, and it filters sequencing error with `--min-count`. The gap was
never in the engine; it was in *how a user spells the request*. A shell glob
hands you `strainB_R1.fastq.gz strainB_R2.fastq.gz` as two unrelated paths, but
they are **one sample**. Something has to know that.

We kept that knowledge out of the engine. `inputs.py` is a small, pure module
whose entire job is to turn "whatever the user typed" into a list of `Sample`s —
a frozen dataclass of `(name, files, is_reads)` — that `engine/ska.py` can render
into SKA's native file list without ever thinking about pairing again. The
separation matters: `build()` takes `Sequence[Sample]`, so it is identical
whether a sample came from a glob, a `--manifest`, or a future input source we
haven't imagined.

The auto-pairing itself is a tuple of regexes tried **most-specific first**:
`_R1/_R2` (with an optional Illumina `_001` lane suffix), then `.R1/.R2`, then the
bare `_1/_2`. Order is the whole game — a greedy `_1/_2` rule would wrongly split
`sample_R1` on the `_1`, so the `_R[12]` pattern has to win first. Two FASTQs that
share a base and carry opposite mate numbers collapse into one read `Sample`; a
lone FASTQ with no mate is still valid (SKA treats single FASTQ as assembly-like).

Two sharp edges worth knowing:

- **A paired sample is one input, not two.** The "need at least two genomes" guard
  counts *samples*, so `strainB_R1 strainB_R2` alone is correctly rejected as a
  single-sample run — a test pins exactly this (`test_single_paired_read_sample_counts_as_one_input`).
- **`--auto-k` needs an assembly.** `kselect` scans plaintext FASTA to score
  *k*-mer uniqueness; it cannot parse FASTQ. So `_resolve_k` filters to the
  *assembly* inputs before probing, and on a **reads-only** run it logs a warning
  and falls back to the fixed `-k` (default 31) rather than silently doing the
  wrong thing. Read-error flags (`--min-count`/`--min-qual`/`--qual-filter`) are
  likewise only emitted when the build actually contains reads — SKA ignores them
  for assemblies, but we don't pass noise we don't need.

> **The throughline:** when a tool you wrap already supports a feature, the work
> is rarely in the engine — it's in the *input grammar*. Give that grammar its own
> module, make it pure and exhaustively tested, and the engine layer stays as dumb
> as it should be.

### 3.9 Logging for the moment it breaks (`logconfig.py`)

A tool that shells out to `ska` and a handful of optional engines has a specific
failure signature: *the useful evidence is the exact command, its stderr, and
which step we were on* — and all three are trivially lost if they only ever
reached the terminal. By the time a run dies on a cluster at 2 a.m., the
scrollback is gone.

So skTree logs to **two sinks with different jobs**. The console handler follows
`-v`/`--debug` for live use. The file handler is the interesting one: it writes
`<outdir>/sktree.log` at **DEBUG, always**, no matter how quiet the console is.
That asymmetry is deliberate — a default, flagless run still leaves a full
post-mortem trail on disk. You never have to "re-run with `--debug`" after the
fact, which is exactly the moment you can't (the inputs may be gone, the failure
may be intermittent).

The single highest-value line in the whole feature lives in `engine/ska.py`'s
`_run`: before every subprocess it logs `running: <shlex-joined argv>`, and on a
non-zero exit it folds **both the command and stderr into the `SkaError`
itself**:

```
ska exited with code 101
  command: ska build -o out/combined -k 31 -f out/combined.filelist.tsv
  stderr : ... real.fasta has no valid sequence
```

Before this, a failed build told you only the return code — you couldn't even see
*what* was run. Now the error is self-contained and the same detail is in the log.

Two design choices worth calling out:

- **`logconfig.configure_logging` is idempotent.** It tears down existing handlers
  before adding new ones, so a second run in the same process (or a test that
  configures twice) never double-writes every line. The file opens in `"w"` mode:
  one log *per run*, not an ever-growing append no one prunes.
- **The CLI is the only caller, and it owns the catch-all.** `run_pipeline` stays
  free of logging-setup concerns; `cli.main` configures the sinks, logs a context
  header (version, argv, resolved `ska` version), and wraps the run in a
  two-tier `except`: known `SkaError`/`ValueError` get a one-line console message,
  anything unexpected gets `logger.exception(...)` so the **traceback lands in the
  file**. Either way the last thing printed is *where the log is*. A crash you
  can't see is a crash you can't fix.

> **The throughline:** logging isn't decoration you sprinkle on at the end — it's
> the contract that the tool will be debuggable by someone who wasn't watching.
> Decide what evidence a failure needs, then guarantee that evidence reaches disk
> *before* you need it.

---

## 4. War stories (bugs and lessons)

These are real things that went wrong during development. They're the most
valuable part of this document.

### 4.1 The 98% divergence that wasn't — *the* lesson of the project

A smoke test screamed **"Max pairwise divergence 98.36%"** on genomes that
differed by a handful of SNPs. The warning was reading the wrong column of
`ska distance`.

`ska distance` emits two different numbers:
- **Distance** — the SNP *count* between two samples.
- **Mismatches (proportion)** — the fraction of compared *split-k-mers* that
  differ.

We naively warned on the second. But here's the split-*k*-mer gotcha: **one SNP
invalidates roughly *k* split-k-mers** (every window that spans it). So the
*k*-mer mismatch proportion runs about *k*-fold higher than the true per-base
divergence. On a genome with very few shared *k*-mers, it can pin to nearly 100%.

The fix: estimate true divergence as `SNP_count / (match_count + mismatch_count)`
— SNPs over the number of comparable split-*k*-mer positions. That gave 0.1%,
matching reality. The lesson, framed generally:

> **Never trust a column name; verify what the number *means*.** A "proportion"
> that's computed in *k*-mer space is not the "proportion" your biology intuition
> expects. We confirmed semantics straight from `ska distance --help` before
> trusting either column.

This is now pinned by `test_max_divergence_uses_snp_rate_not_kmer_proportion`, a
*unit* test that needs no `ska` binary — it feeds the real TSV layout and asserts
we don't re-trip the old false alarm.

### 4.2 The tandem-repeat fixture that found zero SNPs

An early smoke test used a genome built as `"ACGT" * 75` — a pure tandem repeat.
The run found **0 SNPs** and produced IUPAC ambiguity codes (`H`, `Y`) in the
alignment. Panic? No — *garbage in, garbage out*. A tandem repeat has almost no
*unique* split-*k*-mers (they all collide), so SKA can't anchor anything.

> **Lesson:** when validating a *k*-mer tool, your synthetic data must be
> *non-repetitive*, or you're testing the pathology, not the pipeline. We switched
> to a deterministic SHA-256-seeded pseudo-random genome and immediately got
> sensible SNPs and distances. Determinism (seeded, not `random`) keeps the test
> reproducible.

### 4.3 Relative imports in test files

`from .conftest import requires_ska` blew up with *"attempted relative import with
no known parent package."* pytest collects test files as top-level modules, so
they have no package context for relative imports. The fix was to stop sharing
helpers via import and inline the skip guard:
`pytestmark = pytest.mark.skipif(shutil.which("ska") is None, ...)`. Shared
*fixtures* still live in `conftest.py` (pytest injects those by name, no import
needed) — but shared *code* should be imported absolutely or inlined.

### 4.4 The empty-string path that became a directory

`parse_alignment("")` crashed with `IsADirectoryError: Is a directory: '.'`. The
source-type detection fell through to `Path("").read_text()`, and an empty path
resolves to the current directory. The guard now explicitly checks the string is
non-empty *and* names a real file before treating it as a path:
`elif "\n" not in source and source and Path(source).is_file()`. Tiny bug, classic
shape: an empty/edge value silently taking a code path meant for something else.

### 4.5 Fail loud, never silent

A recurring theme (and a house rule): the ML wrapper raises `MlNotAvailable` at
*construction* if no engine is on PATH, so the pipeline decides **once** whether
ML is possible and logs a clear "skipping ML tree: …". It never half-runs and
hides the reason. Likewise `SkaRunner._run` raises `SkaError` with captured
stderr on any non-zero exit. Errors should be impossible to miss.

---

## 5. Technologies, and why these and not others

| Tool | Role | Why it, over the alternative |
|------|------|------------------------------|
| **SKA2** (`ska` Rust) | split-*k*-mer engine | The whole premise. Orders of magnitude faster than a Python port; battle-tested. |
| **DendroPy** | trees & distance matrices | Mature phylogenetics library with a clean `nj_tree()` and Newick I/O. Biopython's `Phylo` is thinner for distance trees. |
| **NumPy** | the SNP matrix | Vectorised column stats turn per-locus loops into C-speed array ops. |
| **Biopython** | FASTA parsing | The de-facto standard; no reason to hand-roll a parser. |
| **IQ-TREE / RAxML-NG** | ML search (optional) | Reinventing a likelihood engine is a multi-year project; wrap the best-in-class instead. |
| **pytest** | tests | Fixtures + markers (`@needs_ska`) make the integration/unit split clean. |
| **uv** | dependency / venv management | Fast, reproducible, lockfile-friendly. |
| **ruff** | lint + format | One fast tool replaces flake8 + isort + more. |
| **hatchling** | build backend | Minimal, standards-based `pyproject` packaging. |
| **TheiaPhyloViewer** | report tree renderer | Pure-ESM, zero-dependency SVG phylo viewer with branch-scaled layout, per-tip styling, a scale bar and built-in pan/zoom. Vendored under its owner's grant; replaced a hand-rolled d3 phylogram and removed d3 entirely. |
| **Vite** (build-time only) | bundle TPV → one IIFE | Library mode tree-shakes TPV's 36-file source to a ~60 KB self-contained global, so the report stays single-file and offline. Not a runtime dependency. |

---

## 6. Testing philosophy

We wrote tests **before** implementation (TDD) and split them two ways:

- **Unit tests** mock the subprocess boundary (`unittest.mock.patch`) or feed
  in-memory data. They run in milliseconds and need *no* external binaries — so
  CI stays green even where `ska`, `iqtree`, and `raxml-ng` aren't installed.
- **Integration tests** drive the *real* `ska` binary and **auto-skip** when it's
  absent (`@pytest.mark.skipif(shutil.which("ska") is None, …)`). On a dev box
  with SKA they exercise the genuine end-to-end path.

This is why the suite has tests like "argv construction for IQ-TREE" (pure unit,
no binary) right next to "full run produces all outputs" (integration). Each
class of bug has a test that can actually catch it.

> **Why mock the boundary and not the science?** Mocking `subprocess.run` lets us
> assert the *exact* command line we build without needing the tool installed —
> the contract we own. We never mock our own SNP/parsimony logic; that we test
> against known-answer inputs (e.g. a hand-computed Fitch score of 1).

---

## 7. How to find your way around

- Start at `cli.py` → `pipeline.run_pipeline()`. Read it top to bottom; it's the
  table of contents for the whole system.
- Each step in the conductor maps to one module in §2.3. Open the module to see
  the algorithm; open its `tests/test_*.py` to see what "correct" means.
- The two design documents — `docs/RESEARCH.md` (why SKA2, how split k-mers work,
  the design rationale) and `docs/PLAN.md` (milestones, output contract) — capture
  the *reasoning* behind the code.

## 8. Known limitations & roadmap

- **Reference-based SNP annotation** is now built (§3.6, `--reference`): `ska map`
  + pyrodigal label each SNP intergenic/coding with a synonymous/non-synonymous
  call. It stays an *optional* path behind the `annotate` extra so the
  reference-free core has zero new required dependencies. Possible next steps:
  multi-codon (MNP) effects, frameshift/indel handling, and gene *names* via a
  supplied annotation (GFF) rather than de-novo prediction.
- **Divergent samples**: SKA's exact-flank matching loses recall beyond ~1%
  divergence. skTree *warns*; it can't fix the underlying limit.
- **Parsimony scale**: pure-Python search is for tens of isolates, not hundreds
  (§3.5). Use `--ml` for large sets.
