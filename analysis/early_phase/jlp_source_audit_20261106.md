# Source audit of the 11 numeric JLP screening records

Checked 6 November 2026 against the archived `jtc/study3.zip`, the existing
`outputs/jlp_cases.csv`, and primary registry records. The companion CSV records
the exact source fields and findings. No eligible records, screening counts,
p-values, sample sizes or Bayes factors were changed by this audit.

The archive has `groupDescription`, `pValue`, `statisticalMethod` and
outcome-level denominator rows. It does **not** contain `pValueComment`,
`nonInferiorityComment`, effect estimates, confidence intervals, or outcome
category-specific denominators. Several consequential details are therefore
absent from the extraction rather than contradicted by its stored values.

## Findings requiring explicit flags

| Records | Flag | Primary-source finding |
|---|---|---|
| J03 | C | Reported two-sided p=.05 accompanies difference −.463 and 95% CI [−.731,−.194]; this does not match an ordinary corresponding t/Wald calculation. |
| J04 | S | Month-24 category gives 344+354=698 evaluable participants for the selected arms; extraction uses overall 485+486=971. |
| J08 | C | OR .334, CI [.199,.938], p=.0373 require reconciliation; the interval is not log-symmetric and ordinary log-Wald reconstruction gives p≈.0056. The CI method is unspecified. |
| J18 | D | Official p-value comment identifies p=.0456 as one-sided, so 2p1=.0912 and the manuscript's directional significance rule is not met. |
| J19–J20 | M | Both p=.044 values are Benjamini–Hochberg adjusted; the SAP specifies moderated t-tests. They are not unadjusted individual-test p-values. |
| J19–J20 | S | Biomarker-specific sample counts are 461+433=894 for CCN4 and 463+433=896 for TGF beta-1, versus overall mBFAS n=921 used by the extraction. |

C denotes a source-method consistency question, not a proven reporting error.
Different confidence-interval constructions can matter. S denotes a documented
denominator mismatch; D denotes direction; M denotes multiplicity adjustment.
These flags supplement the existing bound/equivalence flags and do not imply
that every unflagged record has been fully validated.

The direction and multiplicity comments are particularly explicit:

- [NCT04376684 API](https://clinicaltrials.gov/api/v2/studies/NCT04376684), first Part-1 alive/free-of-respiratory-failure analysis: “p-value is from a one-sided test.” The published [OSCAR report](https://publications.ersnet.org/content/erj/61/2/2101870) reports two-sided p=.09 for the overall Part-1 comparison.
- [NCT05013008 API](https://clinicaltrials.gov/api/v2/studies/NCT05013008), CCN4 and TGF beta-1 analyses: “Benjamini-Hochberg adjusted P-value was calculated using the log2-scaled NPX value.” The individual no-difference null is explicit; the problem is the adjusted input, not an equivalence hypothesis. [SAP section 6.2.1](https://cdn.clinicaltrials.gov/large-docs/08/NCT05013008/SAP_001.pdf) specifies moderated tests and multiplicity correction. [EU registry primary-endpoint rows](https://www.clinicaltrialsregister.eu/ctr-search/trial/2021-003053-37/results) provide the biomarker-specific counts.
- [NCT00405821 API](https://clinicaltrials.gov/api/v2/studies/NCT00405821) provides J03's effect and interval. Its p=.05 is retained as reported; it should not be silently replaced by a value inferred from the interval.
- [NCT00565812 API](https://clinicaltrials.gov/api/v2/studies/NCT00565812) gives J04's Month-24 denominators in the outcome `classes` title, separately from the overall `denoms` counts.

Thus the 11 numeric reports comprise six with specific source flags
(J03, J04, J08, J18, J19 and J20), and five without
these particular flags (J05, J06, J07, J09, J10). The unflagged designation
means only that these checks did not reveal the identified input problems.
The nine bound reports are outside this numeric-input audit.

## Recommended illustrations

**J05: joint-space narrowing.** At Month 24 in NCT00565812, progression was
reported for 136/419 participants assigned SD-6010 50 mg daily and 165/429
assigned placebo. The adjusted logistic OR was .745 (95% CI .557–.997), p=.047.
Using the extracted n=848 gives eJAB=4.06. Progression meant a decrease in
joint-space width exceeding .199 mm. This is a secondary superiority/other
comparison, not an equivalence test. The clinical trial ran during 2007–2011;
these results were posted in 2016. [Registry source](https://clinicaltrials.gov/study/NCT00565812?tab=results)

**J06: major bleeding.** In NCT01225822, the safety-set comparison was
dabigatran 50 mg twice daily versus enoxaparin 40 mg daily, with 1/389 versus
8/392 major bleeding events. The reported OR .121 (95% CI .015–.974) and
p=.0472 are approximately consistent on a log-Wald scale. The extraction gives
n=781 and eJAB=3.91. The outcome window extended from approximately 14 days
before surgery to 4–6 weeks afterward. Only nine events underpin this contrast,
so asymptotic Gaussian approximations warrant particular caution.
[Registry source](https://clinicaltrials.gov/study/NCT01225822?tab=results)

**J09: thromboembolism/death composite.** A second comparison in that same
trial used dabigatran 150 mg twice daily versus enoxaparin 40 mg daily. The
composite counted 49/375 versus 72/382 events; the treatment-and-centre logistic
model gave OR .647 (95% CI .427–.980), p=.0401. Extracted n=757 yields
eJAB=3.36. J10's related VTE endpoint has identical observed counts and reported
contrast; these are not independent replications. Both are superiority/other
analyses. [Registry source](https://clinicaltrials.gov/study/NCT01225822?tab=results)

These are useful **as-reported** illustrations of the difference between a
thresholded significance decision and the evidence supplied by an approximate
Bayes factor. Under the stated eJAB construction, values around 3–4 favor the
point null relative to its specified alternative. They do not establish that a
treatment effect is absent, that a significant finding is a Type I error, or
that a different clinical decision should have been made. For adjusted models,
the outcome-level counts still do not independently verify the exact
complete-case analysis population or all likelihood assumptions.

The original proposed CCN4 example should be replaced in the narrative. Keep
its archived screening row with the multiplicity and sample-size flags. J18
and J03 are also unsuitable uncomplicated substitutes. If only two clinical
illustrations are desired, J05 and J06 provide distinct trials without adding
the dependence between J06 and J09.
