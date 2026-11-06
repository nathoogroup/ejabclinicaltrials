#!/usr/bin/env Rscript

# Audit/catalog of the current early-phase screening results with p_nhst <= .05
# and applicable eJAB > 3. This does not rerun or alter the eligibility or
# posterior analyses. Literal p-value bounds and hypothesis-context concerns
# are retained, rather than silently treating every screened record as a
# verified Jeffreys--Lindley example.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else
  "analysis/early_phase/jlp_catalog.R"
script_dir <- dirname(normalizePath(script_file, mustWork = TRUE))
project_dir <- normalizePath(file.path(script_dir, "../.."), mustWork = TRUE)
output_dir <- file.path(script_dir, "outputs")
draft_dir <- file.path(project_dir, "tmp", "jlp_20261105")
dir.create(draft_dir, recursive = TRUE, showWarnings = FALSE)
eligible_file <- file.path(output_dir, "eligible_early_phase_analyses.csv")
archive_file <- file.path(project_dir, "jtc", "study3.zip")
input_hashes <- tools::md5sum(c(eligible_file, archive_file))

eligible <- read_csv(eligible_file, show_col_types = FALSE, progress = FALSE)
stopifnot(!anyDuplicated(eligible$analysisId))
selected <- eligible %>%
  filter(p_nhst <= 0.05, ejab_01 > 3) %>%
  mutate(phase_order = match(phase, c("Early Phase I", "Phase I", "Phase II"))) %>%
  arrange(phase_order, nctId, outcomeTitle, analysisId) %>%
  mutate(catalog_row = row_number(), catalog_id = sprintf("J%02d", catalog_row), .before = 1)

raw <- read_csv(unz(archive_file, "study3.csv"),
  col_types = cols(.default = col_character()), progress = FALSE,
  show_col_types = FALSE)
raw_selected <- raw %>% filter(analysisId %in% selected$analysisId) %>%
  mutate(value_num = parse_number(value, locale = locale(grouping_mark = ",")))
metadata <- raw_selected %>% group_by(analysisId) %>% summarise(
  briefTitle = first(briefTitle),
  group_description = paste(unique(na.omit(groupDescription)), collapse = " | "),
  denominator_groups = paste(paste0(groupId, "=", value), collapse = "; "),
  denominator_units = paste(unique(units), collapse = "; "),
  group_rows = n(), distinct_group_ids = n_distinct(groupId),
  repeated_group_id = anyDuplicated(groupId) > 0,
  value_sum = sum(value_num), value_first = first(value_num),
  value_max = max(value_num),
  denominator_values_valid = all(is.finite(value_num) & value_num > 0),
  raw_study_type = first(studyType),
  raw_outcome_type = first(outcomeType),
  raw_method = first(statisticalMethod),
  raw_p_reported = first(pValue),
  raw_nctId = first(nctId),
  exact_early_phase_i = first(EARLY_PHASE1) == "1" &
    all(first(PHASE1) == "0", first(PHASE2) == "0", first(PHASE3) == "0", first(PHASE4) == "0"),
  exact_phase_i = first(PHASE1) == "1" &
    all(first(EARLY_PHASE1) == "0", first(PHASE2) == "0", first(PHASE3) == "0", first(PHASE4) == "0"),
  exact_phase_ii = first(PHASE2) == "1" &
    all(first(EARLY_PHASE1) == "0", first(PHASE1) == "0", first(PHASE3) == "0", first(PHASE4) == "0"),
  consistent_analysis_metadata = all(
    n_distinct(nctId) == 1, n_distinct(outcomeTitle) == 1,
    n_distinct(outcomeType) == 1, n_distinct(statisticalMethod) == 1,
    n_distinct(pValue) == 1), .groups = "drop")

catalog <- selected %>% left_join(metadata, by = "analysisId") %>% mutate(
  raw_phase = case_when(exact_early_phase_i ~ "Early Phase I",
    exact_phase_i ~ "Phase I", exact_phase_ii ~ "Phase II", TRUE ~ NA_character_),
  raw_sample_size = case_when(
    test_class %in% c("One-sample/paired t-test", "Wilcoxon signed-rank test") ~ value_first,
    test_class == "Conditional logistic regression" ~ value_max,
    TRUE ~ value_sum),
  raw_dimension = if_else(test_class %in% c("One-way ANOVA", "Kruskal-Wallis test"),
    group_rows - 1L, 1L),
  p_bound_recomputed = as.numeric(str_extract(raw_p_reported,
    "[0-9]+[.]?[0-9]*(?:[eE][-+]?[0-9]+)?")),
  one_sided_recomputed = raw_dimension == 1 & str_detect(str_to_lower(raw_method),
    "one[ -]?sided|1[ -]?sided|one[ -]?tailed|1[ -]?tailed"),
  p_nhst_recomputed = if_else(one_sided_recomputed, 2 * p_bound_recomputed, p_bound_recomputed),
  # Use the survival quantile directly, independently of the pipeline's
  # Q(1-p) implementation, avoiding cancellation in 1-p.
  statistic_p = if_else(one_sided_recomputed,
    2 * pmin(p_bound_recomputed, 1 - p_bound_recomputed), p_bound_recomputed),
  chisq_statistic_recomputed = qchisq(statistic_p, raw_dimension, lower.tail = FALSE),
  ejab_recomputed = if_else(one_sided_recomputed,
    sqrt(raw_sample_size) * exp(-chisq_statistic_recomputed / 2) /
      (2 * (1 - p_bound_recomputed)),
    exp(log(raw_sample_size) / 2 -
      (1 - raw_sample_size^(-1 / raw_dimension)) * chisq_statistic_recomputed / 2)),
  bf_absolute_difference = abs(ejab_recomputed - ejab_01),
  p_reporting = if_else(p_left_censored, "left-censored bound", "numeric report"),
  bf_is_upper_bound = p_left_censored & !one_sided,
  bf_interpretation = if_else(bf_is_upper_bound,
    "Upper bound evaluated at reported p-value bound; BF > 3 is not established",
    "Estimate from numeric reported p-value; subject to reporting precision and model assumptions"),
  hypothesis_context = case_when(
    str_detect(str_to_lower(group_description), "equivalen") ~ "explicit equivalence context",
    analysisId == "0c7899487f9b7c316126c2c964c697ee" ~ "bioequivalence context by cross-reference",
    TRUE ~ "no conflicting hypothesis identified in archived description"),
  hypothesis_context_flag = hypothesis_context !=
    "no conflicting hypothesis identified in archived description",
  tail_information = case_when(
    str_detect(str_to_lower(raw_method), "two[ -]?sided|2[ -]?sided|two[ -]?tailed|2[ -]?tailed") ~
      "Explicit two-sided method label",
    one_sided_recomputed ~ "Explicit one-sided method label",
    TRUE ~ "Tail unspecified in method; pipeline uses two-sided/unrestricted default"),
  hypothesis_context_evidence = case_when(
    hypothesis_context == "explicit equivalence context" ~ group_description,
    hypothesis_context == "bioequivalence context by cross-reference" ~ paste(
      "Description says 'Same as the Day 1 analyses'. Sibling analysis",
      "7a0702f36080dca4f0194692283e8157 explicitly describes two one-sided",
      "bioequivalence hypotheses; ordinary point-null interpretation needs review."),
    TRUE ~ ""),
  denominator_check = if_else(repeated_group_id,
    "Repeated group identifier: denominator needs review",
    "Distinct participant-group rows; sum reproduces pipeline n"),
  denominator_scope_caveat = paste("Counts are outcome-level participant denominators selected by analysis group IDs.",
    "The archived extraction does not establish the complete-case analysis-specific n."),
  registry_url = paste0("https://clinicaltrials.gov/study/", nctId),
  table_flag = case_when(hypothesis_context_flag & p_left_censored ~ "B,E",
    p_left_censored ~ "B", TRUE ~ "N"),
  # Illustrations use numeric reports, distinguish the biomarker identities,
  # and make no claim about unrecorded treatment effects or error status.
  illustrative_rank = case_when(
    analysisId == "9e0c4bca23cae4a61441fa6c37ca8353" ~ 1L,
    analysisId == "0d9665b5c025261d0b134ee96d1d6e11" ~ 2L,
    analysisId == "387f7c1eadf70b5084aea02301b7ade4" ~ 3L,
    analysisId == "d45c0322e68167a7aa6ea2034e5cecec" ~ 4L,
    analysisId == "2828aa15694934e977a66ae338aafd1e" ~ 5L,
    TRUE ~ NA_integer_),
  table_outcome = case_when(
    analysisId == "387f7c1eadf70b5084aea02301b7ade4" ~
      "Plasma biomarker change: CCN family member 4 (O95388)",
    analysisId == "dc6573fbe37eade71e369bcc5b991ac6" ~
      "Plasma biomarker change: TGF beta-1 proprotein (P01137)",
    analysisId == "0d9665b5c025261d0b134ee96d1d6e11" ~
      "Major bleeding events: BIBR 1048 50 mg twice daily versus enoxaparin 40 mg daily",
    analysisId == "392859fe7a68badbe8be57d387feb046" ~
      "Proximal DVT: BIBR 1048 225 mg versus 50 mg, both twice daily",
    analysisId == "c41ae4e13b5aac3224a1efea55c829f9" ~
      "Proximal DVT, pulmonary embolism and VTE-related mortality: BIBR 1048 225 mg versus 50 mg, both twice daily",
    analysisId == "ff42bcee42d42c5b1b488443bffb185b" ~
      "VTE events and all-cause mortality: BIBR 1048 150 mg twice daily versus enoxaparin 40 mg daily",
    analysisId == "05870cbc4e509cba6181955bc11e9514" ~
      "VTE events: BIBR 1048 150 mg twice daily versus enoxaparin 40 mg daily",
    analysisId == "9f2dc6f681f7757f63b93ba61869b07d" ~
      "Meningococcal C hSBA geometric mean titers at day 28: MenC-CRM ROS versus EMV",
    analysisId == "e135961a30dcdd7977f9f150889b4036" ~
      "Meningococcal C hSBA geometric mean titers at day 28: MenC-CRM LIQ versus EMV",
    analysisId == "ef32df42461f4233ed24564185d37aef" ~
      "Meningococcal C hSBA geometric mean titers at day 28: MenC-CRM liquid versus ROS",
    TRUE ~ outcomeTitle))

bad_ids <- c("9f08f829a4d81c7b436b818dd972e1b4", "01bf66ca01ef01f21b9e3bfd18d4c81f",
  "25638018d5a7c1c3cd22fcddd0cf44d8", "b39984376c029bac02dc6de78af004cd")
stopifnot(nrow(catalog) == nrow(selected),
  all(catalog$consistent_analysis_metadata), all(catalog$denominator_values_valid),
  !any(catalog$analysisId %in% bad_ids),
  all(catalog$raw_study_type == "INTERVENTIONAL"),
  all(catalog$raw_outcome_type %in% c("PRIMARY", "SECONDARY")),
  all(catalog$phase == catalog$raw_phase), all(catalog$nctId == catalog$raw_nctId),
  all(catalog$sample_size == catalog$raw_sample_size),
  all(catalog$dimension == catalog$raw_dimension),
  all(catalog$one_sided == catalog$one_sided_recomputed),
  all(catalog$p_nhst == catalog$p_nhst_recomputed),
  all(catalog$p_nhst_recomputed <= .05), all(catalog$ejab_recomputed > 3),
  all(catalog$bf_absolute_difference < 1e-10),
  !any(str_detect(catalog$raw_p_reported, ">")))

# Ensure that the bioequivalence cross-reference is supported by the archive.
if (any(catalog$hypothesis_context == "bioequivalence context by cross-reference")) {
  linked <- raw %>% filter(analysisId == "7a0702f36080dca4f0194692283e8157")
  stopifnot(nrow(linked) > 0,
    all(str_detect(str_to_lower(linked$groupDescription), "bioequivalen")))
}

summary_for <- function(rows, scope, label) tibble(
  scope = scope, label = label, analyses = nrow(rows), trials = n_distinct(rows$nctId),
  numeric_reports = sum(!rows$p_left_censored), left_censored_reports = sum(rows$p_left_censored),
  one_sided_analyses = sum(rows$one_sided),
  hypothesis_context_flags = sum(rows$hypothesis_context_flag),
  repeated_group_id_analyses = sum(rows$repeated_group_id),
  minimum_n = if (nrow(rows)) min(rows$sample_size) else NA_real_,
  maximum_n = if (nrow(rows)) max(rows$sample_size) else NA_real_,
  minimum_ejab = if (nrow(rows)) min(rows$ejab_01) else NA_real_,
  maximum_ejab = if (nrow(rows)) max(rows$ejab_01) else NA_real_)
summary <- bind_rows(summary_for(catalog, "overall", "All screened records"),
  lapply(c("Early Phase I", "Phase I", "Phase II"), function(label)
    summary_for(filter(catalog, phase == label), "phase", label)),
  lapply(sort(unique(catalog$test_class)), function(label)
    summary_for(filter(catalog, test_class == label), "method", label)),
  lapply(sort(unique(catalog$nctId)), function(label)
    summary_for(filter(catalog, nctId == label), "trial", label)))
write_csv(catalog, file.path(output_dir, "jlp_cases.csv"), na = "")
write_csv(summary, file.path(output_dir, "jlp_summary.csv"), na = "")

latex_escape <- function(text) {
  replacements <- c("\\" = "\\textbackslash{}", "&" = "\\&", "%" = "\\%",
    "$" = "\\$", "#" = "\\#", "_" = "\\_", "{" = "\\{", "}" = "\\}",
    "~" = "\\textasciitilde{}", "^" = "\\textasciicircum{}")
  vapply(strsplit(as.character(text), "", fixed = TRUE), function(chars)
    paste0(ifelse(chars %in% names(replacements), replacements[chars], chars), collapse = ""),
    character(1), USE.NAMES = FALSE)
}
latex_p <- function(value) paste0("$", gsub(" ", "", value, fixed = TRUE), "$")
table_rows <- vapply(seq_len(nrow(catalog)), function(i) {
  row <- catalog[i, ]
  method_label <- switch(row$test_class,
    "Two-sample t-test" = "Two-sample $t$-test",
    "Logistic regression" = "Logistic regression", "One-way ANOVA" = "ANOVA",
    latex_escape(row$test_class))
  bf_display <- if (row$bf_is_upper_bound) {
    # Upward rounding preserves the displayed inequality as an upper bound.
    paste0("$\\le ", sprintf("%.2f", ceiling(row$ejab_01 * 100) / 100), "$")
  } else sprintf("%.2f", row$ejab_01)
  paste0("\\textbf{", row$catalog_id, "}\\newline ",
    "\\href{", row$registry_url, "}{", row$nctId, "}\\newline ",
    latex_escape(row$phase), " & ", latex_escape(row$table_outcome), " & ",
    method_label, "\\newline $n=", format(row$sample_size, scientific = FALSE, trim = TRUE),
    "$ & ", latex_p(row$p_reported), " & ",
    bf_display, " & ", row$table_flag, " \\\\")
}, character(1))
table <- c(
  "% Generated by analysis/early_phase/jlp_catalog.R; requires longtable and hyperref.",
  "% All 20 current screening records are retained, including bounds and hypothesis flags.",
  "\\begingroup", "\\small", "\\setlength{\\parskip}{0pt}", "\\setlength{\\tabcolsep}{3pt}",
  "\\begin{longtable}{@{}>{\\raggedright\\arraybackslash}p{.16\\linewidth}>{\\raggedright\\arraybackslash}p{.365\\linewidth}>{\\raggedright\\arraybackslash}p{.13\\linewidth}>{\\raggedright\\arraybackslash}p{.085\\linewidth}>{\\raggedright\\arraybackslash}p{.10\\linewidth}>{\\raggedright\\arraybackslash}p{.06\\linewidth}@{}}",
  paste0("\\caption{Early-phase records with reported NHST significance and computed $\\mathit{eJAB}_{01}>3$.",
    " The screening set includes literal $p$-value bounds; flagged rows require the qualifications below.}",
    "\\label{tab:jlp_cases}\\\\"),
  "\\hline", "Trial / phase & Outcome or comparison & Method / $n$ & Reported\\newline $p$ & $\\mathit{eJAB}_{01}$ & Flag \\\\",
  "\\hline", "\\endfirsthead", "\\hline",
  "Trial / phase & Outcome or comparison & Method / $n$ & Reported\\newline $p$ & $\\mathit{eJAB}_{01}$ & Flag \\\\",
  "\\hline", "\\endhead", "\\hline", "\\endfoot", table_rows,
  "\\end{longtable}",
  paste0("\\noindent N: numeric $p$ report; B: left-censored $p$ report. For B rows, the displayed eJAB is evaluated",
    " at the reported upper bound and is itself an upper bound (rounded upward for display),",
    " so $\\mathit{eJAB}_{01}>3$ is not established.",
    " E: archived descriptions indicate equivalence or bioequivalence hypotheses; an ordinary point-null",
    " interpretation is not justified without further review. NCT02260492 refers to a Day~1 analysis",
    " that explicitly describes two one-sided bioequivalence tests. All listed analyses were classified",
    " as two-sided/unrestricted in the screening pipeline and have $q=1$; significance includes $p=.05$.",
    " Participant counts are sums of the distinct outcome-level groups associated with each analysis;",
    " they do not independently verify the analysis-specific complete-case sample size.",
    " Multiple outcomes or biomarkers within one trial are separate analyses, not independent trials.",
    " Numeric reports remain subject to rounding and the approximation's assumptions.",
    " Row identifiers J01--J20 map to the complete archived analysis identifiers in ",
    "\\texttt{jlp\\_cases.csv}."),
  "\\endgroup")
writeLines(table, file.path(draft_dir, "jlp_cases_table.tex"))

# Document, without changing, the historical export error.
legacy_paths <- file.path(project_dir, "csv", c("JLP05.csv", "CTE05.csv"))
legacy_identical <- all(file.exists(legacy_paths)) &&
  identical(readLines(legacy_paths[1]), readLines(legacy_paths[2]))
audit <- c(
  "Current JLP screening catalog audit",
  paste("Analyses:", nrow(catalog), "Trials:", n_distinct(catalog$nctId)),
  paste("Numeric reports:", sum(!catalog$p_left_censored),
    "Left-censored reports:", sum(catalog$p_left_censored)),
  paste("Hypothesis-context flags:", sum(catalog$hypothesis_context_flag)),
  paste("Maximum independently recomputed BF difference:", format(max(catalog$bf_absolute_difference), digits = 16)),
  paste("Repeated group IDs within an analysis:", sum(catalog$repeated_group_id)),
  paste("Historical csv/JLP05.csv and csv/CTE05.csv identical:", legacy_identical),
  "rmd/ManuscriptFigures.Rmd creates JLP05 but writes CTE05 to csv/JLP05.csv (line 602).",
  "The historical exports are not used by this builder and were left unchanged.",
  "The full archived group descriptions, IDs, denominators, and review flags are in outputs/jlp_cases.csv.",
  "Neither eligible data nor the posterior model was changed.",
  "Input MD5:", paste(names(input_hashes), input_hashes),
  "The case table is an audit draft: upper-bound and equivalence rows are not confirmed JLP cases.")
writeLines(audit, file.path(draft_dir, "jlp_catalog_audit.txt"))
stopifnot(identical(input_hashes, tools::md5sum(c(eligible_file, archive_file))))
print(summary)
cat(paste(audit[1:7], collapse = "\n"), "\n")
