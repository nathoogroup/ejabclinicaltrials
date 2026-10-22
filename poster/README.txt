eJAB methods poster and manuscript draft

Editable documents
  ejab_methods_draft.tex     Methods manuscript with derivations and references
  ejab_methods_poster.tex    One-page poster, 48 x 36 inches, landscape
  shared_math.tex            Equations shared between the two documents
  abstract.tex               Outcome-focused abstract shared by both drafts
  references.bib             Manuscript references and TSBF optimality reference
  figures/                  Decision error profiles, two QR codes and header logos
  citations.txt             Seven citations and DOI links encoded in the QR
  clinical_examples.csv      Three manuscript examples underlying the clinical table
  THRESHOLD_NOTES.txt        TSBF optimality argument, assumptions and sources
  make_citations_qr.py       Rebuilds the QR PNG and SVG from citations.txt
  make_arxiv_qr.py           Rebuilds the separate link to arXiv:2510.10358
  make_header_logos.py       Rebuilds the vector header wordmarks
  render_profiles.R         Regenerates plot titles from saved risk-grid CSVs
  fonts/                    Crimson Pro Bold and its open-font license
  LOGO_NOTES.txt             Logo sources and saved assets

Build from this directory
  latexmk -pdf -interaction=nonstopmode -halt-on-error ejab_methods_draft.tex
  latexmk -pdf -interaction=nonstopmode -halt-on-error ejab_methods_poster.tex

Requires a LaTeX installation with amsmath, lmodern, natbib, graphicx,
geometry, caption and TikZ. The active revised journal manuscript and
reviewer response now also use the updated profile name, Normalized EWDE
legends and BIC comparison. Archived submissions are unchanged.

Scope
The draft contains the one-sided null-centered construction, its explicit
finite-sample convention, regularity conditions and a consistency proof
sketch, Equally Weighted Decision Error Profiles, and a compact early-phase clinical-trial
contradiction table containing the three examples from the revised manuscript.
The poster has four profile panels: A, estimated-effect TSBF; B, true-effect
TSBF; C, p < .05 across ten tests; D, BIC across six parametric tests.
The companion methods draft retains the three-panel main figure and the
separate p < .005 and BIC comparison figures. The revised journal
manuscript shares the poster's four-panel A-D arrangement.

SOURCE_MAP.txt records the source sections and the limited editorial changes.

Authorship and citations update
Zhengxiao Wei is first author, followed by Puneet Velidi. There is no joint
first-authorship or equal-contribution designation in either draft.

The citations QR encodes the seven references in citations.txt as plain text,
including their DOI links. It does not require a hosted references page.
Standalone PNG and SVG copies are in figures/. Regenerate them with
  python3 make_citations_qr.py
using the Python qrcode package.
The separate "Paper on arXiv" QR links directly to
https://arxiv.org/abs/2510.10358 and is also clickable in the PDF.
Rebuild its PNG and SVG with python3 make_arxiv_qr.py.

Header logos
The UVic Mathematics and Statistics mark and the Nathoo Group logo appear
at the top right of the poster. The UVic mark uses the university's official
white vector crest/wordmark with a white department label. The Nathoo Group
mark follows the website's horizontal icon-and-name layout and uses its
Crimson Pro Bold wordmark. Per the user's update, the butterfly retains
the website's original rectangular background and full ASCII artwork.
The UVic mark and Nathoo wordmark sit directly on the navy header.
SOURCE_MAP.txt and LOGO_NOTES.txt record the sources.

The vector logo PDFs are included and require no extra tools for LaTeX
compilation. Rebuild them, if desired, with
  python3 make_header_logos.py
using the Python reportlab, svglib and Pillow packages.

Title and abstract update
Both drafts use the title "Further developments of generalized Jeffreys's
approximate Bayes factor". The abstract follows the user's proposed wording,
with outcomes checked against the existing risk grids. The left column
contains the abstract, one-sided construction and regularity panels; the
results occupy the right column. The profile image has been moved left within
that column to make room for the clinical table, without changing its
physical size, fonts or plotted values. Abstract and body font sizes remain
unchanged.

Updated regularity conditions
The reviewer-response revision work plan contains R1' (left-tail calibration)
and the weakened R2: sqrt(n)*p -> 0 for integer q >= 2; q = 1 retains
-sqrt(n)*p*log(p) -> 0, with p -> 0 explicit. Both drafts now state these
conditions and update the proof. SOURCE_MAP.txt records the exact source.

Clinical-trial table
The table contains the three illustrative cases from the revised manuscript:
NCT00494806 (Early Phase I, n = 66, p < .05, eJAB = 1.23),
NCT03907683 (Phase I, n = 80, p = .037, eJAB = 1.04), and
NCT00539539 (Phase II, n = 1,586, p = .016, eJAB = 2.19).
It shows the clinical comparison, test, sample size, p-value and eJAB.
For the censored first case, eJAB uses the reported bound .05.
These are candidate contradictions for follow-up, not established Type I errors.

Panel A effect-size estimate
The poster includes the plug-in estimator omega_hat = sqrt(X^2/n) and
the resulting scale tau_hat^2 = X^2/q above the decision error profiles, with the
symbols defined. This is transferred from bayes_risk_analysis.Rmd;
SOURCE_MAP.txt records its location. The plotted results are unchanged.

Decision threshold and caption
The results panel combines the manuscript's threshold-1 rationale with
the Neyman-Pearson connection for TSBF. Under the chi-squared model, the
prespecified-effect TSBF gives a UMP upper-tail rule among X^2-based tests
at its attained size. Equal model probabilities and 0-1 loss select
cutoff 1 for minimizing prior-predictive error. THRESHOLD_NOTES.txt gives
the proof and qualifications; Fowlie (2021) is added to the citations QR.
Panel A remains a plug-in sensitivity analysis. The caption is a single
paragraph. All existing font sizes and the physical normalized EWDE figure size are
preserved.

Corner footer
The bottom-right corner is blank. The draft label, date and replacement
sentence have been removed.

Plot naming
"Equally Weighted Decision Error Profiles" is the common name in the
poster and manuscript headings, subtitles, figure titles and captions,
and in bayes_risk_analysis.Rmd. The color scale is labeled "Normalized EWDE".
EWDE is (Type I error + Type II error)/2; the plotted comparison is the
difference in the two rules' EWDE divided by its panelwise maximum
absolute value. The numerical formula is unchanged.
The existing image dimensions, font sizes and normalization are preserved.
The poster's A and B use one row each, followed by the C and D grids, to
fit the additional BIC panel in the same physical figure area.
Plot images were regenerated from the saved result CSVs,
without rerunning simulations or changing risk data. To rebuild them
from a full repository checkout, run Rscript poster/render_profiles.R.
The LaTeX files continue to compile using the bundled PNG assets.

Detailed comparison wording
The abstract now cites Velidi and Wei (2025), Johnson (2005), and
Johnson et al. (2023), and calls the p-value comparisons NHST significance
rules. It distinguishes 2x2/3x3 from 5x5 Pearson independence tests,
estimated from true omega, and the BIC and p < .05 regimes. The abstract
panel is taller; the hypothesis statement is inline and display spacing
is tightened to preserve all poster font sizes and the figure area.
