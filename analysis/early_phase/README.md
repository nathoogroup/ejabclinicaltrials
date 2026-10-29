# Revised early-phase ClinicalTrials.gov analysis

This directory contains the executable analysis used for the revised manuscript
and response to reviewers. It uses the archived August 7, 2025 snapshot in
`jtc/study3.zip` and retains 16,563 primary- and secondary-outcome analyses from
1,790 exact Early Phase I, Phase I, and Phase II trials. The analysis is a
retrospective reanalysis; the screening rules are documented in Appendix C of
`Supplementary Information/SI.pdf`.

For a one-sided scalar test, reconstruct `p = 2 * min(p1, 1-p1)` and compute

```text
eJAB_01 = sqrt(n) * exp(-(1-1/n) * qchisq(1-p, df=1)/2)
eJAB_0+ = eJAB_01 / (2*(1-p1))
```

The one-sided definition applies the Wagenmakers directional factor
`1/[2*(1-p1)]` to the two-sided eJAB and retains its finite-sample factor
`(n-1)/n`. Two-sided/unrestricted tests retain the usual `1 - n^(-1/q)` factor. Significance uses `p1 <= .025` for
one-sided tests and `p <= .05` otherwise, retaining reported inequalities.

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

Numerical p-values are classified as displayed, including equality at the
direction-specific threshold. Literal reports such as `p < .05` retain their
reported inequality and use the numerical upper bound for eJAB. Rounding at
the threshold can obscure the underlying full-precision p-value.

Including equality adds 65 significant results and 32 Bayes/NHST contradictions.
The updated totals are 6,474 significant results, 493 contradictions, and 20
Jeffreys-Lindley cases. The eligible set is unchanged. One-sided eJAB values
now follow the direct multiplicative correction above; two-sided values are
unchanged.
The alpha-sweep figures retain the reconstructed two-sided-p-value convention
and include equality at each alpha threshold.
