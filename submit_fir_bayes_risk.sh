#!/bin/bash
set -euo pipefail

work_dir="$(pwd -P)"
if [[ ! -f "$work_dir/bayes_risk_analysis.Rmd" ]]; then
  echo "Run this script from the ejabclinicaltrials repository root." >&2
  exit 1
fi

run_root="$work_dir/slurm_bic_runs"
mkdir -p "$run_root"
run_dir="$(mktemp -d "$run_root/run-XXXXXXXX")"
shard_dir="$run_dir/shards"
log_dir="$run_dir/logs"
mkdir -p "$shard_dir" "$log_dir"

export EJAB_WORK_DIR="$work_dir"
export EJAB_SHARD_DIR="$shard_dir"

array_submission="$(
  sbatch --parsable \
    --array=1-420 \
    --chdir="$work_dir" \
    --output="$log_dir/%x-%A_%a.out" \
    --error="$log_dir/%x-%A_%a.err" \
    --export=ALL,EJAB_WORK_DIR,EJAB_SHARD_DIR \
    "$work_dir/fir_bic_array.slurm"
)"
array_job_id="${array_submission%%;*}"

merge_submission="$(
  sbatch --parsable \
    --dependency="afterok:$array_job_id" \
    --chdir="$work_dir" \
    --output="$log_dir/%x-%j.out" \
    --error="$log_dir/%x-%j.err" \
    --export=ALL,EJAB_WORK_DIR,EJAB_SHARD_DIR \
    "$work_dir/fir_bic_merge_render.slurm"
)"
merge_job_id="${merge_submission%%;*}"

printf 'run_dir=%s\narray_job_id=%s\nmerge_job_id=%s\n' \
  "$run_dir" "$array_job_id" "$merge_job_id" | tee "$run_dir/submission.txt"
