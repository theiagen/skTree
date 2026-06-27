#!/usr/bin/env bash
# Head-to-head benchmark: skTree vs kSNP4 on one cohort of genomes.
#
# Both tools run in their own container under IDENTICAL caps (same CPUs, same
# memory ceiling with swap disabled) so wall-time and peak-RSS are comparable and
# a real out-of-memory is recorded as a failure rather than hidden by swap. Each
# container reads its own cgroup-v2 memory.peak high-water mark, robust even for
# sub-second runs.
#
# Usage: bench.sh <cohort_dir> <k> <cpus> <mem_gb> <out_root>
set -uo pipefail

COHORT_DIR=${1:?cohort dir}
K=${2:?k}
CPUS=${3:?cpus}
MEM_GB=${4:?mem_gb}
OUT_ROOT=${5:?out_root}

KSNP_IMAGE=${KSNP_IMAGE:-staphb/ksnp4:4.1}
SKTREE_IMAGE=${SKTREE_IMAGE:-sktree:latest}

cohort=$(basename "$COHORT_DIR")
work="$OUT_ROOT/$cohort"
inputs="$work/inputs"
mkdir -p "$inputs"
rm -rf "$work/ksnp_out" "$work/sktree_out"
mkdir -p "$work/ksnp_out" "$work/sktree_out"

# ── stage genomes with clean .fasta names (deref symlinks) ──
shopt -s nullglob
for f in "$COHORT_DIR"/*.fasta "$COHORT_DIR"/*.fa "$COHORT_DIR"/*.fna; do
    name=$(basename "$f"); name=${name%.*}
    cp -L "$f" "$inputs/${name}.fasta"
done
n_genomes=$(ls "$inputs"/*.fasta | wc -l)
echo "[bench] cohort=$cohort genomes=$n_genomes k=$K cpus=$CPUS mem=${MEM_GB}g"

# Containers run as root (kSNP4) / uid-1000 (sktree); make outputs writable by both.
chmod -R 777 "$work"

run_in_container() {
    # run_in_container <label> <image> <extra_docker_args...> -- <inner shell cmd>
    # Force --entrypoint bash so the image's own ENTRYPOINT (e.g. sktree's) does
    # not swallow our wrapper script.
    local label=$1 image=$2; shift 2
    local docker_args=(); while [ "$1" != "--" ]; do docker_args+=("$1"); shift; done; shift
    local inner=$1
    local metrics="$work/${label}.metrics.tsv"
    local wrapped="$inner; rc=\$?; cat /sys/fs/cgroup/memory.peak > /work/${label}.peak 2>/dev/null || echo NA > /work/${label}.peak; exit \$rc"
    local start end wall ec peak
    start=$(date +%s.%N)
    docker run --rm \
        --cpus="$CPUS" --memory="${MEM_GB}g" --memory-swap="${MEM_GB}g" \
        -v "$work:/work" -w /work \
        --entrypoint bash \
        "${docker_args[@]}" "$image" \
        -c "$wrapped" \
        > "$work/${label}.run.log" 2>&1
    ec=$?
    end=$(date +%s.%N)
    wall=$(awk -v e="$end" -v s="$start" 'BEGIN{printf "%.3f", e-s}')
    peak=$(cat "$work/${label}.peak" 2>/dev/null || echo NA)
    printf 'tool\tcohort\tk\tn_genomes\tcpus\tmem_gb\texit_code\twall_seconds\tpeak_rss_bytes\n' > "$metrics"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$label" "$cohort" "$K" "$n_genomes" "$CPUS" "$MEM_GB" "$ec" "$wall" "$peak" >> "$metrics"
    echo "[bench] $label: exit=$ec wall=${wall}s peak=${peak}B"
}

# ── kSNP4: build TSV (abs_path<TAB>name), then kSNP4 -core -vcf ──
: > "$work/ksnp_input.tsv"
for f in "$inputs"/*.fasta; do
    printf '/work/inputs/%s\t%s\n' "$(basename "$f")" "$(basename "$f" .fasta)" >> "$work/ksnp_input.tsv"
done
# Route HOME and TMPDIR onto /work as well: kSNP4 keeps a FilteredKmersCache and
# PyInstaller temp under HOME and writes Jellyfish hashes into the outdir, so
# pointing everything at /work keeps all of kSNP4's scratch on the mounted
# benchmark volume (which may be a dedicated disk) instead of the host root.
mkdir -p "$work/ksnp_home"
run_in_container ksnp4 "$KSNP_IMAGE" \
    -e HOME=/work/ksnp_home -e TMPDIR=/work/ksnp_home -- \
    "kSNP4 -in /work/ksnp_input.tsv -outdir /work/ksnp_out -k $K -core -vcf -CPU $CPUS"

# kSNP4's Jellyfish intermediates dominate disk and accumulate across cohorts;
# scrub them once metrics are read so sequential scaling cells don't fill the disk.
find "$work/ksnp_out" -maxdepth 1 -type f \
    \( -name 'kmers.fsplit*' -o -name 'fsplit*' -o -name '*.Jelly' -o -name 'Jelly.fsplit*' \
       -o -name 'SNPs_all*' -o -name 'core_SNPs_matrix*' \) -delete 2>/dev/null

# ── skTree: same cohort, same k. Tree flags default to parsimony+NJ, but the
# pure-Python parsimony is O(n^3*sites) — set SKTREE_TREE_FLAGS='' for large n to
# run NJ only (NJ is also the tree used for the kSNP4 topology comparison). ──
sktree_inputs=$(for f in "$inputs"/*.fasta; do printf '/work/inputs/%s ' "$(basename "$f")"; done)
run_in_container sktree "$SKTREE_IMAGE" \
    -- \
    "sktree run $sktree_inputs -o /work/sktree_out -k $K ${SKTREE_TREE_FLAGS---parsimony} --threads $CPUS"

echo "[bench] cohort $cohort complete -> $work"
