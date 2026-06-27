// Bundles TheiaPhyloViewer into a single self-contained IIFE exposing `TPV`.
//
// Set TPV_SRC to the absolute path of the TheiaPhyloViewer source entry
// (its `src/index.js`), then:  npx vite build --config vite.config.mjs
//   -> dist/tpv.iife.js   (copy to ../theiaphyloviewer.iife.js)
import { defineConfig } from "vite";
import { fileURLToPath } from "node:url";

const TPV_SRC =
  process.env.TPV_SRC ||
  fileURLToPath(new URL("./tpv-src/src/index.js", import.meta.url));

export default defineConfig({
  resolve: { alias: { TPV_SRC } },
  build: {
    lib: {
      entry: fileURLToPath(new URL("./entry.js", import.meta.url)),
      name: "TPV",
      formats: ["iife"],
      fileName: () => "tpv.iife.js",
    },
    outDir: fileURLToPath(new URL("./dist", import.meta.url)),
    minify: true,
    emptyOutDir: true,
  },
});
