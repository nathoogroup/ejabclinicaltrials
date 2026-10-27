# R2: realized values of D_n in four replication projects

Appendix E of the Supplementary Material, and the table and figure in the response to Reviewer 1.

## Files
- `Scripts/build_strict_data.R`: builds `Data/` from the projects' raw files.
- `Scripts/plot_strict.py`: draws the figure from `Data/`.
- `Data/{EERP,RPP,RPCB,SCORE}_orig.csv`, `_repl.csv`: one row per original result or replication attempt in the H1 subset (`id`, `p`, `logp`, `n`, `censored`).
- `Data/counts.csv`: counts at each step of the selection.
- `Output/Dn_vs_n_points_maxn_uncensored_logy_strict005_FINALFINAL.png`: the figure.

## Rule
- Original p ≤ .005, and every replication attempt p ≤ .005.
- Only p and n are required; the degrees of freedom q are not.
- p-values reported as a bound (`censored = TRUE`) are left out of the figure.
- Details, including how bounds and p = 0 are handled, are in the header of `build_strict_data.R`.

## Raw data (not included)
Place under one directory, `RAW_DIR`:
- `EERP/camerer2016_studydetails.csv`: Camerer et al. (2016), transcribed from `create_studydetails.do`, OSF bzm54.
- `RPP/rpp_data_updates.csv`: Open Science Collaboration (2015), study-level data.
- `RPCB/rpcb_effect_level.csv`: Errington et al. (2021), final analysis, effect level.
- `SCORE/raw/analyst data.RData`: Tyner et al. (2025), OSF g5sny.

## Run
```
Rscript Scripts/build_strict_data.R RAW_DIR
python Scripts/plot_strict.py      # pandas, matplotlib, statsmodels
```
`plot_strict.py` also writes four variants of the figure (linear axis, quantile band, polynomial fit, bounds kept) that are not used.
