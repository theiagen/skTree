// IIFE bundle entry: re-export TheiaPhyloViewer so Vite library mode exposes it
// as the `TPV` global. TPV_SRC is resolved at build time (see vite.config.mjs)
// to the TheiaPhyloViewer source's `src/index.js`.
import { TheiaPhyloViewer } from "TPV_SRC";
export { TheiaPhyloViewer };
