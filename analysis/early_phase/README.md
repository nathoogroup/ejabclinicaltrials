# Revised early-phase ClinicalTrials.gov analysis

This directory contains the executable analysis used for the revised manuscript
and response to reviewers. It uses the archived August 7, 2025 snapshot in
`jtc/study3.zip` and retains 16,563 primary- and secondary-outcome analyses from
1,790 exact Early Phase I, Phase I, and Phase II trials. The analysis is a
retrospective reanalysis; the screening rules are documented in Appendix C of
the companion manuscript repository's `REVISION/SI.tex`.

For a one-sided scalar test, reconstruct `p = 2 * min(p1, 1-p1)` and compute

```text
eJAB_0+ = sqrt(n) * exp(-qchisq(1-p, df=1)/2) / (2*(1-p1))
```

The one-sided definition combines the null-centered approximation, omitting
`(n-1)/n` from its exponent, with the Wagenmakers directional factor
`1/[2*(1-p1)]`. Two-sided/unrestricted tests retain the usual `1 - n^(-1/q)` factor. Significance uses `p1 <= .025` for
one-sided tests and `p <= .05` otherwise, retaining reported inequalities.

Run from the repository root, with R packages `dplyr`, `ggplot2`, `readr`,
`stringr`, `tidyr`, `lme4`, `cowplot`, `ggpubr`, and `httr` installed:

```sh
Rscript analysis/early_phase/early_phase_analysis.R
Rscript analysis/early_phase/jlp_catalog.R
Rscript analysis/early_phase/application_figure_v3.R
Rscript analysis/early_phase/supplementary_application_figures.R
Rscript analysis/early_phase/verify_one_sided.R
```

The active V3 renderer writes `Figure3_V3.png` for `REVISION/REVISED_V3.tex`;
`application_figure_revised.R` and `Figure3.png` retain the V2 version.
The V3 application renderer currently builds the combined Figure 3: A and C share the top
row, scatterplot B spans the middle row, and outcome-distribution panels D-G
occupy the bottom row. It sources `figure4_revised.R` for the outcome panels;
that helper can still be run directly to reproduce the legacy Figure 4 layout.

The scripts write analysis-level values, screening counts, phase/test summaries,
trial-cluster bootstrap results, prior-audit draws, and screening catalogs to `outputs/`, and manuscript,
reviewer-response, and supplementary figures to `figures/`. The cached MeSH
mapping in `outputs/mesh_tree_numbers.csv` avoids network requests during figure
generation. The scripts use fixed random seeds for the bootstrap and prior-audit summaries.

Numerical p-values are classified as displayed, including equality at the
applicable NHST threshold. Literal reports such as `p < .05` retain their
reported inequality and use the numerical upper bound for eJAB. Rounding at
the threshold can obscure the underlying full-precision p-value.

Including equality adds 65 significant results and 32 Bayes/NHST contradictions.
The updated totals are 6,474 significant results and 493 contradiction flags.
The stronger screen (`p_nhst <= .05` and applicable `eJAB > 3`) returns 20
records from nine trials, comprising 11 numerical reports and nine left-censored
reports. The latter yield upper bounds on eJAB and do not establish `eJAB > 3`;
four also describe equivalence or bioequivalence hypotheses. Thus these are
screening records, not 20 confirmed Jeffreys--Lindley examples.
The eligible set is unchanged. One-sided eJAB values
follow the null-centered approximation above; two-sided values are
unchanged.
All clinical NHST comparisons and plots use `p_nhst = p` for two-sided tests
and `p_nhst = 2*p1` for one-sided tests. The alpha sweep counts a contradiction
when `p_nhst <= alpha` and the applicable eJAB exceeds one. Thus a one-sided
test must satisfy `p1 <= alpha/2`, preserving its specified direction. The
comparison value can exceed one; it is not a reconstructed two-sided p-value.
The folded `p_two_sided = 2*min(p1, 1-p1)` is used only to reconstruct the
unsigned chi-squared statistic inside eJAB. Using it for NHST would incorrectly
admit opposite-tail findings. At alpha .05 this formerly added two Phase II
findings to the sweep, although the reported total of 493 contradictions
already used the directional rule and is unchanged. Both alpha-sweep figures
now read the same `early_phase_alpha_profile.csv` output.

Both one-sided and two-sided/unrestricted tests contribute to the phase summaries.
At alpha .05, significant denominators are 94, 927,
and 5,453 in Early Phase I, Phase I, and Phase II, respectively, including 4, 27,
and 104 one-sided results. The observed contradiction counts are 6, 28, and 459;
none of the 135 significant one-sided results has eJAB greater than one.

The response's Figure 2 (manuscript Figure 3A) plots `C_g(alpha) / N_g`,
where `N_g` is all eligible analyses in phase `g` and `C_g(alpha)` counts
NHST-significant results with the applicable eJAB above one. The response's
Figure 3 (manuscript Figure 3C) reports the observed conditional proportion
`C_g(.05) / S_g(.05)`, where `S_g` counts significant results. For example,
Early Phase I contributes 6/270 overall and 6/94 among significant results.

The shared `clinical_phase_summary.R` implementation supplies both primary
figures. It resamples whole eligible trials uniformly with replacement within
each phase, retains all their analyses together, and recomputes the ratio of
total contradictions to total significant records. It uses 20,000 replicates
and seed 20261105; the 2.5th and 97.5th percentiles give the displayed intervals.
Trials with no significant results remain in the resampling frame. A zero
significant denominator is undefined and triggers a review rather than being
counted as zero; none occurred in this run. The intervals are 1.23--13.75%,
1.14--5.86%, and 6.62--10.12% for Early Phase I, Phase I, and Phase II.
They describe resampling variability, not posterior probabilities or uncertainty
from selective reporting, censored values, or the eJAB approximation.

Outputs `clinical_phase_bootstrap_summary.csv`, `clinical_phase_trial_counts.csv`,
`clinical_phase_bootstrap_draws.csv`, and `clinical_phase_leave_one_trial_out.csv`
retain the summary, cluster counts, resamples, and influence checks. The primary
standalone figure is `early_phase_observed.png`. The observed proportions give
each significant analysis equal weight; they do not average trial-specific
proportions. Leave-one-trial-out ranges are diagnostics, not confidence intervals.

The previous independent per-study Dirichlet(1,1,1,1) model is retained for
audit and prior sensitivity, with its existing output filenames preserved.
Its prior-sensitive means (approximately 24.1%, 16.3%, and 18.1%) are not the
primary phase estimates. The audited densities appear in supplementary
`FigureS10_posterior_audit.png`, with dashed lines marking observed proportions.
Appendix C.8 and `analysis/audits/posterior_20261104/` document the matched
original/revised threshold comparisons and the repeated prior contribution.
`posterior_density_plot.R` remains a renderer for these legacy audit densities.

The JLP screening catalog is reproduced by `jlp_catalog.R`, which reads the
eligible set and archived group metadata without changing either. It writes
`outputs/jlp_cases.csv` and `outputs/jlp_summary.csv`, and generates
`tmp/jlp_20261105/jlp_cases_table.tex` for supplementary Table S3
(`REVISION/JLP_cases_table.tex`, label `tab:jlp_cases`). Row codes J01--J20 map
to complete archived analysis identifiers. The catalog preserves reported
inequalities, group denominators, test-tail information, hypothesis-context
flags, and independent BF checks. The table rounds upper bounds upward for
display. Numerical reports remain subject to rounding and do not independently
verify every hypothesis or analysis-specific sample size. Historical
`csv/JLP05.csv` is not used: its old export mistakenly wrote `CTE05` instead of
`JLP05`; the archived files are preserved.

The next figure revision will retain both contradictions and JLPs, with phase-specific
threshold curves and kernel-density displays. That redesign is still in progress;
the current V3 checkpoint retains the verified observed-proportion panel.
