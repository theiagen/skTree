# syntax=docker/dockerfile:1

# skTree bundles two toolchains: the Rust `ska` binary (the split-k-mer engine)
# and the Python `sktree` package that conducts it. We compile ska in a throwaway
# Rust stage and copy only the resulting binary into a slim Python runtime, so the
# final image ships no Rust toolchain. Both stages share the same Debian release
# (bookworm) so the binary's glibc linkage stays valid after the copy.

ARG SKA_VERSION=0.5.1
ARG FASTBAPS_VERSION=0.1.1
ARG PYTHON_VERSION=3.12

# --- Stage 1: build the SKA2 engine -----------------------------------------
FROM rust:1-slim-bookworm AS ska-builder
ARG SKA_VERSION
RUN cargo install ska --version "${SKA_VERSION}" --locked --root /opt/ska

# --- Stage 2: Python runtime ------------------------------------------------
FROM python:${PYTHON_VERSION}-slim-bookworm AS runtime
ARG FASTBAPS_VERSION

# iqtree provides the optional ML engine so `sktree run --ml` works out of the
# box; without it the pipeline skips ML gracefully rather than failing.
RUN apt-get update \
    && apt-get install -y --no-install-recommends iqtree \
    && rm -rf /var/lib/apt/lists/*

COPY --from=ska-builder /opt/ska/bin/ska /usr/local/bin/ska

WORKDIR /opt/sktree
COPY pyproject.toml README.md ./
COPY sktree ./sktree

# Install the package with the annotate extra, plus fastbaps for the optional
# `--cluster` / `--html` population-structure step. Both pyrodigal and fastbaps
# ship manylinux wheels (fastbaps is an abi3 PyO3 wheel exposing a `fastbaps`
# console script on PATH), so no compiler is needed in the runtime stage.
RUN pip install --no-cache-dir ".[annotate]" "fastbaps==${FASTBAPS_VERSION}"

# Run as an unprivileged user; /data is the conventional bind-mount for inputs.
RUN useradd --create-home --uid 1000 sktree
USER sktree
WORKDIR /data

ENTRYPOINT ["sktree"]
CMD ["--help"]
