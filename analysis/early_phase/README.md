# Revised early-phase ClinicalTrials.gov analysis

This directory contains the executable analysis used for the revised manuscript
and response to reviewers. It uses the archived August 7, 2025 snapshot in
`jtc/study3.zip` and retains 16,563 primary- and secondary-outcome analyses from
1,790 exact Early Phase I, Phase I, and Phase II trials. The analysis is a
retrospective reanalysis; the screening rules are documented in Appendix C of
the companion manuscript repository's `SI/SI.tex`.

For a one-sided scalar test, reconstruct `p = 2 * min(p1, 1-p1)` and compute

```text
eJAB_0+ = sqrt(n) * exp(-qchisq(1-p, df=1)/2) / (2*(1-p1))
```

The one-sided definition combines the null-centered approximation, omitting
`(n-1)/n` from its exponent, with the Wagenmakers directional factor
`1/[2*(1-p1)]`. Two-sided/unrestricted tests retain the usual `1 - n^(-1/q)` factor. Significance uses `p1 <= .025` for
one-sided tests and `p <= .05` otherwise, retaining reported inequalities.

Run from the repository root, with R packages `dplyr`, `ggplot2`, `readr`,
`stringr`, `tidyr`, `lme4`, `cowplot`, `ggpubr`, and `httr` installed. To regenerate
the active clinical figures and posterior summaries from the retained eligible
set:

```sh
Rscript analysis/early_phase/jlp_catalog.R
python3 analysis/early_phase/render_jlp_tables.py
Rscript analysis/early_phase/application_figure_v3.R
Rscript analysis/early_phase/supplementary_application_figures.R
Rscript analysis/early_phase/posterior_prior_sensitivity.R
python3 analysis/early_phase/render_prior_results.py
Rscript analysis/early_phase/verify_one_sided.R
```

`early_phase_analysis.R` rebuilds the eligible set, counts, and significance-level
sweep from the archived data; it also regenerates legacy summaries that are not
included in the active documents. The scripts write numerical outputs to
`outputs/` and figures to `figures/`. The cached MeSH mapping in
`outputs/mesh_tree_numbers.csv` avoids network requests during figure generation.

The active manuscript is `REVISION/REVISED_V3.tex`, with response
`response_to_reviewers/Response_to_Reviewers_V2.tex`. Their result tables and
figures are consolidated into these main TeX files. `render_jlp_tables.py` and
`render_prior_results.py` write regenerated manuscript/response fragments to
the ignored `tmp/early_phase_latex/` directory at the repository root for review
and manual transfer; they do not overwrite either consolidated document or
recreate its archived inputs. The V3 application renderer
writes `figures/Figure3_V3.png`, copied to `REVISION/Figure3_V3.png` for
Supplementary Figure S7. Its top row contains the BF-threshold profile (A) and
scatterplot (B); its bottom row contains primary-outcome panels C--D and
secondary-outcome panels E--F, with unadjusted and adjusted results adjacent.
It sources `figure4_revised.R` for these outcome panels. The same renderer writes
`figures/Clinical_scatter_large.png`, copied to `SI/Clinical_scatter_large.png`
for the standalone enlarged Supplementary Figure S8, including its inset and
legends.

The supplement source is `SI/SI.tex`. Appendix C retains the full inclusion and
exclusion criteria and the JLP screening catalogue. Appendix F contains only the
small-a posterior sensitivity figure and table; Appendix G contains the NHST
significance-level sweep (Figure S6); Appendix H contains Figures S7--S8.
`SI/ONLINE_SUPPLEMENT.tex` contains detailed preprocessing/code and auxiliary
MeSH and diagnostic plots (Online Figures O1--O5), using the current generated
`FigureS4.png`--`FigureS8.png` assets copied into `SI/`. These asset filenames do
not determine their displayed figure numbers.

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
Clinical plots use reported p-values for two-sided/unrestricted tests and
`p_two_sided = 2*min(p1, 1-p1)` for one-sided tests. The scatterplots span
`0 <= p <= 1`, and the adjusted outcome panels use logarithms of these same
plotted values. NHST classifications retain the specified test direction:
`p <= alpha` for two-sided/unrestricted tests and `p1 <= alpha/2` for one-sided
tests. The stored `p_nhst` column encodes this classification rule; it is not the
plotted p-value. Neither folding the plotted values nor enlarging the scatterplot
changes the retained cohort, eJAB values, or significance flags.

Both one-sided and two-sided/unrestricted tests contribute to the phase summaries.
At alpha .05, significant denominators are 94, 927,
and 5,453 in Early Phase I, Phase I, and Phase II, respectively, including 4, 27,
and 104 one-sided results. The observed contradiction counts are 6, 28, and 459;
none of the 135 significant one-sided results has eJAB greater than one.

Supplementary Figure S6 (Appendix G), also included in the response, plots
`C_g(alpha) / N_g`, where `N_g` counts all eligible analyses in phase `g` and
`C_g(alpha)` counts significant results with applicable eJAB above one. It uses
`early_phase_alpha_profile.csv`; its image retains the filename
`FigureS11_alpha_profile.png`. Supplementary Figure S7A instead fixes the NHST
rules at `.05` (two-sided/unrestricted) and `.025` (one-sided), and plots
`R_g(K) / N_g` against the BF threshold `K`. Here `R_g(K)` counts significant
records with strict `eJAB > K`. Markers at `K = 1` and `K = 3` identify
contradictions and the stronger JLP screen; the inset enlarges thresholds 3--5.

The posterior displays in Supplementary Appendix F and the active response use
independent per-trial Dirichlet(a,a,a,a) priors with `a = 1e-2, 1e-4, 1e-6`.
Each phase panel overlays kernel densities from 50,000 posterior draws per prior,
with dashed lines marking the empirical proportion of significant records with
eJAB above one. Trial contributions retain eligible-analysis-count weights.
`posterior_prior_sensitivity.R` writes `prior_sensitivity_*.csv` outputs and
`figures/FigureS13_prior_sensitivity.png`; the latter is copied as
`SI/Prior_sensitivity.png` and `response_to_reviewers/Figures/Prior_sensitivity.png`.
`render_prior_results.py` regenerates the live `SI/Prior_sensitivity_results.tex`
input and writes `tmp/early_phase_latex/Prior_sensitivity_results_response.tex`
as a review copy for the consolidated response. These are proportions of analysis records,
not proportions of trials with a flag.

The JLP screening catalog is reproduced by `jlp_catalog.R`, which reads the
eligible set and archived group metadata without changing either. It writes
`outputs/jlp_cases.csv` and `outputs/jlp_summary.csv`, and generates
`tmp/jlp_20261105/jlp_cases_table.tex` for supplementary Table S3
(`SI/JLP_cases_table.tex`, label `tab:jlp_cases`). Row codes J01--J20 map
to complete archived analysis identifiers. The catalog preserves reported
inequalities, group denominators, test-tail information, hypothesis-context
flags, and independent BF checks. The table rounds upper bounds upward for
display. Numerical reports remain subject to rounding and do not independently
verify every hypothesis or analysis-specific sample size. Historical
`csv/JLP05.csv` is not used: its old export mistakenly wrote `CTE05` instead of
`JLP05`; the archived files are preserved.

The stronger-screen phase counts are 0/94, 2/927, and 18/5,453 among significant
records, or 0/270, 2/2,098, and 18/14,195 overall. Both Phase I records are
bound-only reports; Phase II has 11 numerical and seven bound-only reports.

Historical bootstrap and a=1 posterior code, draws, and audit outputs are retained
for reproducibility, including `clinical_phase_summary.R`, legacy functions in
`jlp_phase_panels.R`, `posterior_density_plot.R`, and
`analysis/audits/posterior_20261104/`. They are not displayed in the active
manuscript, response, SI, or online supplement; the former SI discrepancy audit
(Appendix C.8) was removed. The active V3 renderer uses only the threshold-plot
helper from `jlp_phase_panels.R` and does not recompute bootstrap results.
Farouk's preserved V2 manuscript and original response now reside in the
companion repository's `REVISION/cache/` and `response_to_reviewers/cache/`.
