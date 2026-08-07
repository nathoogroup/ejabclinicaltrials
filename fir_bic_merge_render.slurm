#!/bin/bash
#SBATCH --job-name=ejab-bic-merge
#SBATCH --account=def-nathoo_cpu
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -euo pipefail

module load r/4.4.0 r-bundle-bioconductor/3.20
export R_LIBS_USER="${SLURM_TMPDIR}/empty-r-library"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
mkdir -p "$R_LIBS_USER"

: "${EJAB_WORK_DIR:?EJAB_WORK_DIR must be set}"
: "${EJAB_SHARD_DIR:?EJAB_SHARD_DIR must be set}"

Rscript "$EJAB_WORK_DIR/fir_bic_merge_render.R" \
  "$EJAB_SHARD_DIR" "$EJAB_WORK_DIR"
