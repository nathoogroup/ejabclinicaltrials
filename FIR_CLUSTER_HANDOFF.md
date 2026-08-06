# Fir handoff: dense Bayes-risk analysis at BF01 = 1

## Goal

Run the Bayes-risk analysis on Fir with the common decision threshold
`BF01 <= 1` for eJAB, BIC, and TSBF. Use Slurm task parallelism for the dense
empirical BIC calculation, merge the task outputs, and render the final HTML,
tables, and plots.

The source of truth is `bayes_risk_analysis.Rmd`. It already contains the
BF01 = 1 source changes, including the threshold-dependent analytic TSBF
critical values and plot labels. The expensive BIC cache has not been fully
regenerated locally.

## Required implementation on Fir

1. Parameterize the empirical BIC grid in `bayes_risk_analysis.Rmd`. A useful
   dense target is 48 logarithmically spaced sample sizes from 100 to 10,000
   and 64 effect sizes from 0.002 to 0.9 (20 logarithmic points from 0.002 to
   0.02 plus 44 linear points from 0.025 to 0.9).
2. Extract or expose the BIC task functions so one Slurm array task can run a
   deterministic subset of the Cartesian product of the seven BIC test classes
   and sample-size grid. Each task should evaluate every effect size for its
   assigned `(test_key, n)` cell.
3. Use a Slurm job array, not a single-node `mclapply`, for the expensive BIC
   grid. One array element per `(test_key, n)` cell is straightforward; chunked
   ranges are also fine if Fir limits array size. Default each array element to
   one CPU and set BLAS/OpenMP thread counts to one.
4. Preserve the existing deterministic seed rule:
   `params$seed + 100000L + task_row_number`. A retry must reproduce the same
   shard.
5. Write one RDS shard per array element to a dedicated temporary/shard
   directory. Write atomically (temporary file followed by rename) so a
   preempted task cannot leave a valid-looking partial shard.
6. Submit a merge/render job with an `afterok` dependency on the array. Before
   merging, verify all expected task IDs are present, the `(test_key, n,
   theta)` keys are unique and complete, row counts match the requested grid,
   and required BIC risks are finite.
7. Save the merged cache at
   `bayes_risk_outputs/revised_three_comparisons/large_sample_empirical_bic_simulation.rds`,
   then render `bayes_risk_analysis.Rmd` with `bic_refresh = false`. Refresh the
   non-BIC simulation and dense analytic grid as needed so every result uses
   BF01 = 1.

## Current statistical settings

- BIC test classes: two-sample t, linear regression, logistic regression, Cox
  proportional hazards, one-way ANOVA, repeated-measures ANOVA, and chi-square
  independence.
- BIC null replicates per `(test_key, n)` task: 2,000.
- BIC alternative replicates per effect-size cell: 400.
- Equal prior model probabilities and equal error costs.
- Rejection rule for eJAB, BIC, and TSBF: `BF01 <= 1`.
- Cox uses the observed number of events as `n_eff`.

## Suggested Slurm layout

- `fir_bic_array.slurm`: maps `SLURM_ARRAY_TASK_ID` to one or more rows of the
  BIC task grid and creates RDS shards.
- `fir_bic_merge_render.slurm`: validates and merges shards, writes the cache,
  and renders the report.
- `submit_fir_bayes_risk.sh`: submits the array, captures its job ID, and
  submits the merge/render job with `--dependency=afterok:<array_job_id>`.

Add Fir-specific account, partition, quality-of-service, module, scratch, and
time/memory directives based on the current site documentation. Keep all
runtime logs and shards under a job-specific directory so concurrent or retried
runs cannot mix results.

## Outputs to return to the repository

- `bayes_risk_analysis.html`
- `bayes_risk_outputs/revised_three_comparisons/large_sample_empirical_bic_simulation.rds`
- `bayes_risk_outputs/revised_three_comparisons/large_sample_empirical_bic_risk_grid.csv`
- `bayes_risk_outputs/revised_three_comparisons/dense_bayes_risk_grid.rds`
- `bayes_risk_outputs/revised_three_comparisons/dense_bayes_risk_grid.csv`
- All four comparison PNGs, the focused chi-square PNG, and their scaling and
  validation CSVs in the same output directory.

## Acceptance checks

- No `1/3`, `log(3)`, or `log(12)` decision thresholds remain in the active Rmd
  source or rendered plot labels.
- The rank-test thresholds satisfy BF01 = 1 at their critical values:
  `sqrt(2) * exp(-c^2 / 4) = 1` for Wilcoxon/Mann-Whitney and
  `4 * exp(-c / 4) = 1` for Kruskal-Wallis.
- The merged BIC cache has exactly one row per requested
  `(test_key, n, theta)` cell and no missing/non-finite in-scope BIC risks.
- Every figure title and comparison label displays `BF01 <= 1`.
- The R Markdown render completes without errors, and `git diff --check`
  passes before committing generated results.
