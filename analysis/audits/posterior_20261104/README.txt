Posterior discrepancy audit, 2026-11-04

This directory reproduces the initial submission and revision at eJAB > 1/3
and eJAB > 1, separates the original p < .05 boundary bug from the repeated
study-level prior effect, and traces the change in empirical proportions.
Production data, figures, and manuscript files are not modified.

Run from the repository root:
  Rscript analysis/audits/posterior_20261104/audit.R
  python analysis/audits/posterior_20261104/verify.py
  Rscript analysis/audits/posterior_20261104/plot_audit.R
  python analysis/audits/posterior_20261104/build_report.py

Dependencies: R packages dplyr, tidyr, readr, stringr, ggplot2, scales, and
gtools (needed by the unmodified original posterior function); Python pandas,
NumPy, SciPy and ReportLab. In this session, missing dependencies were installed
in isolated tmp/posterior_audit_Rlib and tmp/posterior_audit_venv directories.

The original function is loaded from Figure 3/setup.RData and run unchanged.
Its 5,000-draw result and the raw preprocessing reconstruction are cached to
avoid repeating expensive work. To audit changed source data or source code,
use a new dated audit directory or invalidate those caches explicitly.
The original and revised record counts are asserted, and current output columns
and the complete eligibility flow are compared with a fresh reconstruction.
The raw CSV reader uses base R's connection reader to avoid readr's temporary
275 MB extraction on a disk with very little free space. Both produced the
same final eligible data and flow. The original ZIP passes its CRC check.

Main outputs:
  output/pdf/posterior_audit_20261104.pdf (four-page report)
  results/posterior_comparison.csv (empirical, posterior, prior decomposition)
  results/cohort_bridge.csv (isolate cohort and calculation changes)
  results/independent_prior_decomposition.csv (independent Python arithmetic)
  results/legacy_exact_summary.csv (unmodified original function)
  results/source_manifest.json (repository commits and SHA-256 source hashes)
  audit.log (replication progress and checks)

Definitions:
  p_rule is applied identically to the empirical count and its posterior.
  The original as-executed code uses p < .05; other comparisons use <= .05.
  Revision p means p_nhst, equal to p for two-sided tests and 2*p1 for one-sided.
  There are no BF values exactly equal to 1/3 or 1 in either cohort.
  prior_per_cell is the a in independent Dirichlet(a,a,a,a) study priors.
  posterior_sd describes posterior uncertainty. mcse_mean is simulation error.
  ratio_of_expectations is E(numerator)/E(denominator), not E(their ratio).
  pooled_beta_* columns are diagnostic, ignoring within-study clustering.
  Prior sensitivity values are diagnostic and are not chosen replacement priors.

The manuscript needs an estimand/model decision before using these high
posterior values to make phase-comparison claims. This audit makes no such
replacement and performs no git commit or push.
