# Revised early-phase ClinicalTrials.gov analysis

This directory contains the executable analysis used for the revised manuscript
and response to reviewers. It uses the archived August 7, 2025 snapshot in
`jtc/study3.zip` and retains 16,563 primary- and secondary-outcome analyses from
1,790 exact Early Phase I, Phase I, and Phase II trials. The analysis is a
retrospective reanalysis; the screening rules are documented in Appendix C of
`Supplementary Information/SI.pdf`.

For a one-sided scalar test, reconstruct `p = 2 * min(p1, 1-p1)` and compute

```text
eJAB_0+ = sqrt(n) * exp(-qchisq(1-p, df=1)/2) / (2*(1-p1))
```

The exponent omits the finite-sample factor `(n-1)/n`. The Wagenmakers
directional factor `1/[2*(1-p1)]` is retained. Two-sided/unrestricted tests
retain the usual `1 - n^(-1/q)` factor. Significance uses `p1 < .025` for
one-sided tests and `p < .05` otherwise, retaining reported inequalities.

Run from the repository root, with R packages `dplyr`, `ggplot2`, `readr`,
`stringr`, `tidyr`, `lme4`, `cowplot`, `ggpubr`, and `httr` installed:

```sh
Rscript analysis/early_phase/early_phase_analysis.R
Rscript analysis/early_phase/application_figure_revised.R
Rscript analysis/early_phase/figure4_revised.R
Rscript analysis/early_phase/supplementary_application_figures.R
Rscript analysis/early_phase/verify_one_sided.R
```

The scripts write analysis-level values, screening counts, phase/test summaries,
posterior draws, and the candidate catalog to `outputs/`, and manuscript,
reviewer-response, and supplementary figures to `figures/`. The cached MeSH
mapping in `outputs/mesh_tree_numbers.csv` avoids network requests during figure
generation. The scripts use fixed random seeds for the posterior summaries.

Removing the one-sided finite-sample factor changes the one-sided Bayes factors
but leaves the observed totals at 461 Bayes/NHST contradictions and 19
Jeffreys–Lindley cases. The two-sided values and significance classifications
are unchanged. The alpha-sweep figures retain the existing reconstructed
two-sided-p-value convention; this update changes only the one-sided Bayes
factor calculation.
