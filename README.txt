

Current early-phase clinical-trial reanalysis

The revised manuscript application is reproduced by the scripts in
analysis/early_phase/. See analysis/early_phase/README.md for the run order,
input data, generated results, and one-sided eJAB definition. The one-sided
calculation omits (n-1)/n from the exponent and retains the Wagenmakers
factor 1/[2(1-p1)]. The two-sided calculation retains its finite-sample factor.
The earlier top-level CTGJAB.R and figure scripts describe the original
analysis; use analysis/early_phase/ for the revised early-phase application.

Guidelines for Using the Provided HTML Files

Please download the HTML files and view them using the following instructions to ensure proper functionality:


(1) Rendering LaTeX Expressions: These will not render correctly without an active internet connection. Please ensure that you are connected to the internet when viewing the files.


(2) Equation Display Issues on macOS Ventura: Users accessing these files via the Chrome browser on macOS Ventura might experience issues with equation display. To resolve this, right-click on any equation, navigate to "Math Settings" > "Math Renderer", and select an alternative renderer such as "Common HTML".


(3) Hyperlinks in RPubs Versions: If you are viewing the RPubs version of HTML-knitted R Markdown files, some hyperlinks may not be directly clickable. To open these links, right-click on the hyperlink and select "Open link in new tab".


(4) Reproducibility: Stan results are contingent upon the consistency of several configurations, including the Stan version, computer hardware, and the C++ compiler. https://mc-stan.org/docs/reference-manual/reproducibility.html

See also `.RNG.name` in JAGS. https://www.rdocumentation.org/packages/rjags/versions/4-15/topics/jags.model
