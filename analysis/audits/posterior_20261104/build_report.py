"""Build a reviewable PDF from the completed replication, not manuscript edits."""
from pathlib import Path
import json
import pandas as pd
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak, Image

ROOT = Path(__file__).resolve().parents[3]
RESULTS = Path(__file__).parent / "results"
DEST = ROOT / "output/pdf/posterior_audit_20261104.pdf"
d = pd.read_csv(RESULTS / "posterior_comparison.csv")
decomp = pd.read_csv(RESULTS / "independent_prior_decomposition.csv")
bridge = pd.read_csv(RESULTS / "cohort_bridge.csv")
manifest = json.loads((RESULTS / "source_manifest.json").read_text())
legacy_result = pd.read_csv(RESULTS / "legacy_exact_summary.csv").iloc[0]
W = 7.25 * inch
styles = getSampleStyleSheet()
styles.add(ParagraphStyle(name="AuditTitle", fontName="Helvetica-Bold", fontSize=18, leading=22, spaceAfter=10, textColor=colors.HexColor("#17324D")))
styles.add(ParagraphStyle(name="AuditHead", fontName="Helvetica-Bold", fontSize=12, leading=15, spaceBefore=10, spaceAfter=7))
styles.add(ParagraphStyle(name="AuditBody", fontName="Helvetica", fontSize=9.8, leading=13.5, spaceAfter=7))
styles.add(ParagraphStyle(name="AuditSmall", fontName="Helvetica", fontSize=8.3, leading=11, spaceAfter=5))
styles.add(ParagraphStyle(name="Cell", fontName="Helvetica", fontSize=8.3, leading=10.8))
styles.add(ParagraphStyle(name="CellHead", fontName="Helvetica-Bold", fontSize=8.1, leading=10.5, textColor=colors.white))
story = []
def p(text, style="AuditBody"):
    story.append(Paragraph(text, styles[style]))
def h(text): p(text, "AuditHead")
def table(headers, rows, widths):
    values = [[Paragraph(str(x), styles["CellHead"]) for x in headers]]
    values += [[Paragraph(str(x), styles["Cell"]) for x in row] for row in rows]
    t = Table(values, colWidths=widths, repeatRows=1, hAlign="LEFT")
    t.setStyle(TableStyle([
        ("BACKGROUND", (0,0), (-1,0), colors.HexColor("#17324D")),
        ("VALIGN", (0,0), (-1,-1), "TOP"),
        ("LEFTPADDING", (0,0), (-1,-1), 6), ("RIGHTPADDING", (0,0), (-1,-1), 6),
        ("TOPPADDING", (0,0), (-1,-1), 6), ("BOTTOMPADDING", (0,0), (-1,-1), 6),
        ("LINEBELOW", (0,0), (-1,0), .6, colors.HexColor("#17324D")),
        ("ROWBACKGROUNDS", (0,1), (-1,-1), [colors.white, colors.HexColor("#F1F4F6")]),
        ("LINEBELOW", (0,-1), (-1,-1), .4, colors.HexColor("#A6B3BD")),
    ]))
    story.extend([t, Spacer(1,8)])
def case(dataset, group, threshold, rule="<= .05", prior=1):
    z = d[(d.dataset==dataset)&(d.group==group)&(abs(d.threshold-threshold)<1e-12)&(d.p_rule==rule)&(d.prior_per_cell==prior)]
    assert len(z)==1
    return z.iloc[0]
def pct(x, digits=2): return f"{100*x:.{digits}f}%"
def count_rate(z): return f"{int(z.flagged):,}/{int(z.significant):,}<br/>{pct(z.empirical)}"
def posterior(z): return f"{pct(z.posterior_mean)}<br/>[{pct(z.lower)}, {pct(z.upper)}]"

p("Audit of empirical and posterior proportions", "AuditTitle")
p("ClinicalTrials.gov analysis | 4 November 2026", "AuditSmall")
p("<b>Finding:</b> the original counts and posterior are reproducible, and the current counts and posterior are also reproducible. The current densities use eJAB &gt; 1. Their displacement from the empirical proportions is caused by the repeated study-level priors and the aggregation used by the model. A separate p-value boundary error is present in the original posterior code.")
h("1. Compare like with like")
p("Every row below uses p<sub>NHST</sub> &le; .05 for both the empirical count and posterior. In the revision, p<sub>NHST</sub> = p for two-sided/unrestricted tests and 2p<sub>1</sub> for one-sided tests. Posterior entries are means [95% credible intervals] from 20,000 independent draws, using the existing per-study Dirichlet(1,1,1,1) model.")
rows=[]
for dataset,group,label in [("Original","All","Original: all trials"),("Revision","All","Revision: all early phases"),("Revision","Early Phase I","Revision: Early Phase I"),("Revision","Phase I","Revision: Phase I"),("Revision","Phase II","Revision: Phase II")]:
    a=case(dataset,group,1/3); b=case(dataset,group,1)
    rows.append([label,count_rate(a),posterior(a),count_rate(b),posterior(b)])
table(["Cohort","eJAB &gt; 1/3<br/>Observed","eJAB &gt; 1/3<br/>Posterior","eJAB &gt; 1<br/>Observed","eJAB &gt; 1<br/>Posterior"],rows,[112,103,103,101,103])
p("eJAB &gt; 1 favors the null. The broader event eJAB &gt; 1/3 also includes weak evidence favoring the alternative; it is not the same event as a contradiction in direction of evidence.","AuditSmall")
h("2. What the initial submission actually computed")
p(f"The unmodified function saved with the original Figure 3 data reproduces a posterior mean of <b>{pct(legacy_result.posterior_mean,4)}</b> and posterior SD <b>{100*legacy_result.posterior_sd:.4f} percentage points</b> (5,000 draws), matching the reported 35.5% (0.22%). However, it uses <b>p &lt; .05</b> and <b>eJAB &ge; 1/3</b>. The figure and paper describe p &le; .05.")
a=case("Original","All",1/3,"< .05"); b=case("Original","All",1/3)
table(["Original analysis at eJAB &gt; 1/3","Observed","Posterior mean"],[["As executed: p &lt; .05",count_rate(a),pct(legacy_result.posterior_mean)],["Matched to stated rule: p &le; .05",count_rate(b),pct(b.posterior_mean)]],[235,147,140])
p("There are 3,285 original results at the stored p-value boundary .05, including 1,690 reported as bounds. Of the 3,285, 3,281 have eJAB &gt; 1/3 and 1,468 have eJAB &gt; 1. No archived or revised eJAB equals exactly 1/3 or 1, so &gt; versus &ge; at the BF threshold makes no numerical difference here.","AuditSmall")

story.append(PageBreak())
p("Why the Phase II posterior moves to 18.1%", "AuditTitle")
p("The relevant sample size for each study's categorical model is its number of eligible reported analyses, not its clinical participant count. Phase II has 14,195 eligible analyses across <b>1,463 separately modeled studies</b>. The median study has <b>3 eligible analyses and 1 significant analysis</b>; 580 studies have no significant analyses. Four prior pseudocounts are added separately to each study.")
h("3. Exact accounting of the prior contribution")
p("Let m<sub>i</sub> be all eligible analyses, s<sub>i</sub> the significant analyses, and c<sub>i</sub> the significant analyses with eJAB &gt; 1 in study i. With prior parameter a in each category, the model draws each study's four probabilities and computes:")
p("<b>R = sum[m<sub>i</sub> P<sub>i4</sub>] / sum[m<sub>i</sub>(P<sub>i3</sub> + P<sub>i4</sub>)]</b>")
p("Writing w<sub>i</sub> = m<sub>i</sub> / (m<sub>i</sub> + 4a), the exact expected weighted numerator is sum[w<sub>i</sub> c<sub>i</sub>] + a sum[w<sub>i</sub>], and the expected weighted denominator is sum[w<sub>i</sub> s<sub>i</sub>] + 2a sum[w<sub>i</sub>]. Their ratio approximates the posterior mean of R; it is not an identity for the expectation of a ratio.")
ph2=decomp[decomp.phase=="Phase II"].iloc[0]
table(["Phase II quantity (a = 1)","Observed-count contribution","Prior contribution","Total"],[
    ["Expected weighted numerator",f"{ph2.weighted_data_num:,.2f}",f"{ph2.weighted_prior_num:,.2f}",f"{ph2.weighted_data_num+ph2.weighted_prior_num:,.2f}"],
    ["Expected weighted denominator",f"{ph2.weighted_data_den:,.2f}",f"{ph2.weighted_prior_den:,.2f}",f"{ph2.weighted_data_den+ph2.weighted_prior_den:,.2f}"],
],[162,132,114,114])
p(f"The ratio is <b>{pct(ph2.ratio_of_expectations)}</b>, agreeing with the simulated posterior mean of <b>{pct(case('Revision','Phase II',1).posterior_mean)}</b>. The empirical value is 459/5,453 = <b>8.42%</b>. The prior supplies {pct(ph2.prior_share_in_expected_sig_mass,1)} of the expected weighted significant mass and pulls its conditional contradiction proportion toward 50%.")
p("This is not a single small prior applied to 5,453 significant results. There are 5,852 nominal pseudocounts across the Phase II study models, including 2,926 in the two significant categories. Weighting reduces their effective contribution to the quantities above. Studies with no observed significant results still contribute positive posterior mass to both significant categories.")
h("4. Prior sensitivity with identical data and thresholds")
rows=[]
for ph in ["Early Phase I","Phase I","Phase II"]:
    b=case("Revision",ph,1)
    rows.append([ph,pct(b.empirical),pct(b.posterior_mean),pct(case("Revision",ph,1,prior=.1).posterior_mean),pct(case("Revision",ph,1,prior=.01).posterior_mean)])
table(["Phase","Observed","a = 1","a = 0.1","a = 0.01"],rows,[130,98,98,98,98])
p("All posterior columns use Dirichlet(a,a,a,a) independently for each study. These are sensitivity diagnostics, not selected replacements. A prior should be justified by its implications for the target quantity; choosing a solely to reproduce the empirical percentage would not justify the analysis.","AuditSmall")
p("As another diagnostic, a single pooled Beta(1,1) prior on the significant-result proportion gives 460/5,455 = 8.43% in Phase II. This ignores clustering and is not proposed as a final replacement model.","AuditSmall")

story.append(PageBreak())
p("Where the old-to-new decrease comes from", "AuditTitle")
p("The revision reprocesses and restricts the same August 7, 2025 archived extraction. It is not a newly collected independent CTG sample. The original cohort covers 71,126 analyses from 9,002 studies; the revision contains 16,563 analyses from 1,790 studies.")
story.append(Image(str(RESULTS / "empirical_posterior_comparison.png"), width=W, height=W*4.2/9))
p("Matched comparisons at p<sub>NHST</sub> &le; .05. Diamonds show empirical proportions; blue points and intervals show posterior means and 95% credible intervals under the existing model. Connecting lines emphasize the difference. The original as-executed p &lt; .05 analysis is shown separately on page 1.","AuditSmall")
h("5. Hold the eJAB threshold fixed while changing the cohort")
labels=["Original full cohort","Original restricted to interventional exact early phases and primary/secondary outcomes","Original values on shared analysis IDs","Revised values on shared analysis IDs","Full revised cohort"]
rows=[]
for label,step in zip(labels,bridge.step.drop_duplicates()):
    a=bridge[(bridge.step==step)&(bridge.threshold<1)].iloc[0]
    b=bridge[(bridge.step==step)&(bridge.threshold==1)].iloc[0]
    rows.append([label,f"{int(a.analyses):,}",f"{int(a.significant):,}",f"{int(b.flagged):,} ({pct(b.empirical)})",f"{int(a.flagged):,} ({pct(a.empirical)})"])
table(["Cohort / calculation","Eligible","Significant","eJAB &gt; 1","eJAB &gt; 1/3"],rows,[218,60,65,89,90])
p("At eJAB &gt; 1, the broad cohort restriction alone reduces 13.28% to 7.82%; full revised preprocessing yields 7.62%. There are 16,463 shared analysis IDs. Their NHST comparison values are unchanged. Of these, 348 eJAB values change; six results leave the eJAB &gt; 1 category and 29 leave the eJAB &gt; 1/3 category. The revision adds 100 eligible IDs, four significant, with one above 1 and three above 1/3.","AuditSmall")

story.append(PageBreak())
p("Verification and implications for the revision", "AuditTitle")
h("6. Checks completed")
for text in [
    "Reproduced the original 71,126 analyses, 30,790 significant results, and 4,088 contradictions exactly from the archived Figure 3 dataset. Ran its unmodified saved posterior function at the original threshold and significance rule.",
    "Rebuilt all current eligibility and eJAB calculations from the raw zipped CSV using the production preprocessing expressions. Every column of the resulting 16,563 eligible records agrees with the cached current dataset. The raw ZIP integrity check passed.",
    "Used an independent Python/SciPy implementation to reconstruct every original and revised eJAB and every revised NHST/contradiction indicator. All checks passed. Independent Python and R counts also agree for the cohort comparison.",
    "Replicated the posterior using a separately implemented gamma/Dirichlet sampler with 20,000 draws per case. Comparisons with the saved current posterior and the unmodified original posterior pass Monte Carlo error checks. Posterior standard deviations and Monte Carlo standard errors of the means are stored separately.",
    "Checked BF threshold ties, p-value boundary counts, per-study category counts, weighting, and the prior contribution analytically. The revised significant set combines both test types: 135 significant results are one-sided, and none of those has eJAB &gt; 1.",
]: p("- " + text)
h("7. What needs correction or reconsideration")
p("<b>Original code/caption inconsistency:</b> the original posterior uses p &lt; .05 while the caption and empirical count use p &le; .05. Applying the stated inclusive rule changes the old eJAB &gt; 1/3 posterior from about 35.5% to about 39.8%.")
p("<b>Current threshold:</b> the revised phase densities really use eJAB &gt; 1. The separate roughly 35% statement uses eJAB &gt; 1/3. Confusing these two events explains part of the historical comparison, but does not explain the current same-threshold discrepancy.")
p("<b>Model sensitivity:</b> the implementation assigns independent fixed Dirichlet priors to studies; it does not learn common hyperparameters or partially pool study proportions. Calling the study-level prior flat does not make its aggregate influence small. The current phase comparisons are heavily affected by these repeated priors.")
p("<b>Interpretation:</b> the empirical percentage is a descriptive proportion in the observed eligible records. The current posterior concerns a weighted mixture of latent study category probabilities conditional on significance. A valid posterior under that specified model is not, by itself, evidence that the model is appropriate for the manuscript's inferential purpose.")
p("The defensible next step is to retain the verified empirical counts, resolve the target estimand, and justify a model and prior with an explicit sensitivity analysis before relying on the posterior phase comparisons. The current high posterior percentages should not be defended by the one-sentence pseudocount explanation alone. No manuscript, production figure, or production data was changed by this audit.")
h("8. Reproducibility")
p("Analysis repository: " + manifest["root_commit"][:12] + ". Manuscript repository: " + manifest["manuscript_commit"][:12] + ".", "AuditSmall")
p("Sources: original Figure 3/setup.RData and Figure 3/Figure 3.R; rfuncs/DirichletMultinomial.R; manuscript/submitted1022.tex; jtc/study3.zip; analysis/early_phase/early_phase_analysis.R and application_figure_revised.R. The source manifest records SHA-256 hashes.", "AuditSmall")
p("Audit directory: analysis/audits/posterior_20261104. audit.R contains the replication and sensitivity calculations; verify.py independently checks the arithmetic. Results include posterior_comparison.csv, cohort_bridge.csv, independent_prior_decomposition.csv, draw files, and source_manifest.json. Seeds: 20261104 for the original function, 20261105 for comparison simulations. The pooled Beta calculation is a diagnostic only.","AuditSmall")

def footer(canvas, doc):
    canvas.setStrokeColor(colors.HexColor("#C6D0D8")); canvas.setLineWidth(.4)
    canvas.line(45,36,567,36)
    canvas.setFont("Helvetica",8); canvas.setFillColor(colors.HexColor("#4C5E6D"))
    canvas.drawString(45,24,"eJAB posterior audit | 4 November 2026")
    canvas.drawRightString(567,24,str(doc.page))
DEST.parent.mkdir(parents=True,exist_ok=True)
doc=SimpleDocTemplate(str(DEST),pagesize=(612,792),rightMargin=45,leftMargin=45,topMargin=40,bottomMargin=48,
    title="Audit of empirical and posterior contradiction proportions",author="eJAB analysis")
doc.build(story,onFirstPage=footer,onLaterPages=footer)
print(DEST)
