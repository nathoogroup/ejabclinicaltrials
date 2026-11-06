"""Independent arithmetic and source provenance for the posterior audit."""
from pathlib import Path
import hashlib
import json
import subprocess

import numpy as np
import pandas as pd
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).parent / "results"
old = pd.read_csv(OUT / "original_analysis_records.csv")
new = pd.read_csv(ROOT / "analysis/early_phase/outputs/eligible_early_phase_analyses.csv")

# Independent reconstruction of both archived and current Bayes factors.
q_old = np.where(old.model.isin(["anova", "kruskal_wallis", "repeated_measures"]), old.I - 1, 1)
q_old = np.where(old.model.eq("chi_squared"), (old.R - 1) * (old.C - 1), q_old)
old_bf = np.sqrt(old.N) * np.exp(-.5 * (1 - old.N ** (-1 / q_old)) * chi2.ppf(1 - old.pValue, q_old))
assert np.allclose(old_bf, old.JAB, rtol=1e-9, atol=1e-12)
p_fold = np.where(new.one_sided, 2 * np.minimum(new.p_bound, 1 - new.p_bound), new.p_bound)
p_nhst = np.where(new.one_sided, 2 * new.p_bound, new.p_bound)
new_bf = np.sqrt(new.sample_size) * np.exp(-.5 * (1 - new.sample_size ** (-1 / new.dimension)) * chi2.ppf(1 - p_fold, new.dimension))
one_bf = np.sqrt(new.sample_size) * np.exp(-.5 * chi2.ppf(1 - p_fold, 1)) / (2 * (1 - new.p_bound))
new_bf = np.where(new.one_sided, one_bf, new_bf)
assert np.allclose(new_bf, new.ejab_01, rtol=1e-9, atol=1e-12)
assert np.allclose(p_nhst, new.p_nhst, rtol=0, atol=1e-15)
assert np.array_equal(p_nhst <= .05, new.statistically_significant)
assert np.array_equal((p_nhst <= .05) & (new_bf > 1), new.contradiction)

details = []
for phase, d in new.groupby("phase", sort=False):
    ct = d.assign(s=d.p_nhst.le(.05), c=d.p_nhst.le(.05) & d.ejab_01.gt(1)).groupby("nctId").agg(m=("analysisId", "size"), s=("s", "sum"), c=("c", "sum"))
    w = ct.m / (ct.m + 4)
    A, B, H = (w * ct.c).sum(), (w * ct.s).sum(), w.sum()
    details.append(dict(phase=phase, studies=len(ct), analyses=int(ct.m.sum()),
        significant=int(ct.s.sum()), contradictions=int(ct.c.sum()), raw=ct.c.sum()/ct.s.sum(),
        median_all_per_study=float(ct.m.median()), median_sig_per_study=float(ct.s.median()),
        studies_no_sig=int((ct.s==0).sum()), studies_one_result=int((ct.m==1).sum()),
        studies_le4_results=int((ct.m<=4).sum()),
        nominal_pseudocounts_all=4*len(ct), nominal_pseudocounts_sig=2*len(ct),
        weighted_data_num=A, weighted_data_den=B, weighted_prior_num=H, weighted_prior_den=2*H,
        data_only_ratio=A/B, prior_share_in_expected_sig_mass=2*H/(B+2*H),
        ratio_of_expectations=(A+H)/(B+2*H),
        expected_flagged_mass_no_observed_sig=float(w[ct.s==0].sum()),
        p05_count=int((d.p_nhst==.05).sum()),
        censored_p05_count=int(((d.p_nhst==.05) & d.p_left_censored).sum())))
pd.DataFrame(details).to_csv(OUT / "independent_prior_decomposition.csv", index=False)
print(pd.DataFrame(details).to_string(index=False))

phase_columns = ["EARLY_PHASE1", "PHASE1", "PHASE2", "PHASE3", "PHASE4"]
flags = old[phase_columns].fillna(0).astype(int)
exact_early = flags.iloc[:, :3].sum(axis=1).eq(1) & flags.sum(axis=1).eq(1)
restricted = old.loc[exact_early & old.studyType.eq("INTERVENTIONAL") & old.outcomeType.isin(["PRIMARY", "SECONDARY"])]
bridge = []
for label, d, pcol, bfcol in [
    ("Original full", old, "pValue", "JAB"),
    ("Original restricted eligibility", restricted, "pValue", "JAB"),
    ("Original on shared IDs", old.loc[old.analysisId.isin(new.analysisId)], "pValue", "JAB"),
    ("Revised on shared IDs", new.loc[new.analysisId.isin(old.analysisId)], "p_nhst", "ejab_01"),
    ("Revision full", new, "p_nhst", "ejab_01"),
]:
    for t in (1/3, 1):
        s = d[pcol].le(.05)
        c = s & d[bfcol].gt(t)
        bridge.append(dict(step=label, threshold=t, analyses=len(d), significant=int(s.sum()), flagged=int(c.sum()), raw=float(c.sum()/s.sum())))
pd.DataFrame(bridge).to_csv(OUT / "independent_cohort_bridge.csv", index=False)
print(pd.DataFrame(bridge).to_string(index=False))

files = [ROOT / "Figure 3/setup.RData", ROOT / "jtc/study3.zip",
         ROOT / "rfuncs/DirichletMultinomial.R", ROOT / "Figure 3/Figure 3.R",
         ROOT / "analysis/early_phase/early_phase_analysis.R",
         ROOT / "analysis/early_phase/application_figure_revised.R"]
files.extend(sorted((ROOT / "analysis/early_phase/outputs").glob("*.csv")))
man = ROOT / "Reassessing_the_Statistical_Evidence_in_Cancer_Research_Using_Approximate_Objective_Bayes_Factors"
files.extend([man / "manuscript/submitted1022.tex", man / "REVISION/REVISED.tex",
              man / "REVISION/SI.tex", man / "response_to_reviewers/Response_to_Reviewers.tex"])
manifest = dict(root_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                manuscript_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=man, text=True).strip(),
                sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
                checks="Independent Python reconstruction of all old/current eJABs, directional NHST flags, contradiction flags, and exact prior decomposition passed.")
(OUT / "source_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
print(manifest["checks"])
