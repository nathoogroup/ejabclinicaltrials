#!/usr/bin/env Rscript
# Reproduce the submitted and revised posterior calculations without modifying
# production data, figures, or manuscript files. Run from the repository root.
audit_lib <- file.path(getwd(), "tmp/posterior_audit_Rlib")
if (dir.exists(audit_lib)) .libPaths(c(audit_lib, .libPaths()))
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})
root <- normalizePath(getwd())
out <- file.path(root, "analysis/audits/posterior_20261104/results")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
options(width = 160)
say <- function(...) { cat(..., "\n"); flush.console() }
legacy <- new.env()
load(file.path(root, "Figure 3/setup.RData"), envir = legacy)
old <- as.data.frame(legacy$study3_summary)
stopifnot(nrow(old) == 71126L, !anyDuplicated(old$analysisId),
          sum(old$pValue <= .05) == 30790L,
          sum(old$pValue <= .05 & old$JAB > 1) == 4088L)
write_csv(old, file.path(out, "original_analysis_records.csv"))

legacy_file <- file.path(out, "legacy_exact_replication.rds")
if (!file.exists(legacy_file)) {
  say("Running the unmodified posterior function stored in the original Figure 3 archive (5,000 draws).")
  set.seed(20261104)
  original_result <- legacy$posterior_analysis(old, n_sim = 5000, alpha = .05)
  saveRDS(original_result, legacy_file)
} else original_result <- readRDS(legacy_file)
say("Original exact function: posterior mean", original_result$mean,
    "posterior SD", original_result$se)
write_csv(tibble(posterior_mean = original_result$mean,
                 posterior_sd = original_result$se, draws = 5000L),
          file.path(out, "legacy_exact_summary.csv"))

# Evaluate the actual production preprocessing expressions in an isolated
# environment, stopping before output writing, mixed models, and plotting.
fresh_file <- file.path(out, "fresh_preprocessing_streamed.rds")
if (!file.exists(fresh_file)) {
  say("Rebuilding revision eligibility and eJAB values from jtc/study3.zip.")
  fresh <- new.env(parent = globalenv())
  fresh$input_zip <- file.path(root, "jtc/study3.zip")
  # readr's connection readers spool the 275 MB CSV to disk.
  # Use base R's connection reader because local free disk is limited.
  # All preprocessing expressions and supplied parsing options are unchanged;
  # every resulting cached column is checked below.
  fresh$read_csv <- function(file, ...) tibble::as_tibble(utils::read.csv(
    file, colClasses = "character", check.names = FALSE, na.strings = c("", "NA"),
    stringsAsFactors = FALSE, strip.white = TRUE, comment.char = ""
  ))
  needed <- c("bad_analysis_ids", "classify_method", "raw", "analysis_level",
              "flow_steps", "keep", "flow", "eligible")
  exprs <- parse(file.path(root, "analysis/early_phase/early_phase_analysis.R"))
  for (ex in exprs) {
    if (is.call(ex) && identical(ex[[1]], as.name("<-")) &&
        is.symbol(ex[[2]]) && as.character(ex[[2]]) %in% needed) eval(ex, fresh)
  }
  saveRDS(list(eligible = fresh$eligible, flow = fresh$flow), fresh_file)
  rm(fresh); gc()
}
fresh_result <- readRDS(fresh_file)
revised <- read_csv(file.path(root, "analysis/early_phase/outputs/eligible_early_phase_analyses.csv"), show_col_types = FALSE)
rebuilt <- fresh_result$eligible %>% select(all_of(names(revised)))
stopifnot(identical(rebuilt$analysisId, revised$analysisId))
for (nm in names(revised)) {
  stopifnot(isTRUE(all.equal(rebuilt[[nm]], revised[[nm]], check.attributes = FALSE, tolerance = 1e-12)))
}
stopifnot(nrow(revised) == 16563L, !anyDuplicated(revised$analysisId),
          sum(revised$statistically_significant) == 6474L,
          sum(revised$contradiction) == 493L)
cached_flow <- read_csv(file.path(root, "analysis/early_phase/outputs/inclusion_exclusion_flow.csv"), show_col_types = FALSE)
stopifnot(isTRUE(all.equal(as.data.frame(fresh_result$flow), as.data.frame(cached_flow), check.attributes = FALSE)))
say("Revision rebuilt from raw archive: all cached columns agree; 16,563 analyses, 6,474 significant, 493 contradictions.")
write_csv(fresh_result$flow, file.path(out, "recomputed_eligibility_flow.csv"))

old_standard <- old %>% transmute(analysisId, nctId, p = pValue, bf = JAB,
  phase = case_when(
    coalesce(EARLY_PHASE1, 0) == 1 & coalesce(PHASE1, 0) == 0 & coalesce(PHASE2, 0) == 0 & coalesce(PHASE3, 0) == 0 & coalesce(PHASE4, 0) == 0 ~ "Early Phase I",
    coalesce(PHASE1, 0) == 1 & coalesce(EARLY_PHASE1, 0) == 0 & coalesce(PHASE2, 0) == 0 & coalesce(PHASE3, 0) == 0 & coalesce(PHASE4, 0) == 0 ~ "Phase I",
    coalesce(PHASE2, 0) == 1 & coalesce(EARLY_PHASE1, 0) == 0 & coalesce(PHASE1, 0) == 0 & coalesce(PHASE3, 0) == 0 & coalesce(PHASE4, 0) == 0 ~ "Phase II",
    TRUE ~ "Other"), studyType, outcomeType, has_less_than)
new_standard <- revised %>% transmute(analysisId, nctId, p = p_nhst, bf = ejab_01, phase)

count_studies <- function(d, threshold, inclusive = TRUE) {
  sig <- if (inclusive) d$p <= .05 else d$p < .05
  fav <- d$bf > threshold
  tibble(nctId = d$nctId, sig, fav) %>% group_by(nctId) %>%
    summarise(m = n(), y1 = sum(!sig & !fav), y2 = sum(!sig & fav),
              y3 = sum(sig & !fav), y4 = sum(sig & fav), .groups = "drop")
}
draw_model <- function(ct, prior = 1, draws = 20000L) {
  num <- den <- numeric(draws)
  for (i in seq_len(nrow(ct))) {
    a <- prior + as.numeric(ct[i, c("y1", "y2", "y3", "y4")])
    z <- matrix(rgamma(draws * 4L, shape = rep(a, each = draws)), ncol = 4L)
    z <- z / rowSums(z)
    num <- num + ct$m[i] * z[, 4]
    den <- den + ct$m[i] * (z[, 3] + z[, 4])
  }
  num / den
}
summarize_case <- function(d, dataset, group, threshold, inclusive = TRUE,
                           prior = 1, draws = 20000L, seed = 20261105) {
  ct <- count_studies(d, threshold, inclusive)
  S <- sum(ct$y3 + ct$y4); C <- sum(ct$y4)
  shrink <- ct$m / (ct$m + 4 * prior)
  # Exact expectations of the weighted numerator/denominator. Their ratio
  # approximates E[num/den] and is explicitly labelled as such, not equated.
  num_data <- sum(shrink * ct$y4)
  num_prior <- sum(shrink * prior)
  den_data <- sum(shrink * (ct$y3 + ct$y4))
  den_prior <- 2 * num_prior
  set.seed(seed)
  dr <- draw_model(ct, prior, draws)
  key <- paste(dataset, gsub(" ", "_", group), if (threshold == 1) "gt1" else "gtThird",
               if (inclusive) "le" else "lt", prior, sep = "_")
  saveRDS(dr, file.path(out, paste0(key, "_draws.rds")))
  write_csv(ct, file.path(out, paste0(key, "_study_counts.csv")))
  z <- tibble(dataset, group, threshold, p_rule = if (inclusive) "<= .05" else "< .05",
    prior_per_cell = prior, analyses = nrow(d), studies = nrow(ct),
    significant = S, flagged = C, empirical = C / S,
    posterior_mean = mean(dr), posterior_sd = sd(dr), mcse_mean = sd(dr) / sqrt(draws),
    lower = unname(quantile(dr, .025)), upper = unname(quantile(dr, .975)),
    weighted_data_numerator = num_data, weighted_prior_numerator = num_prior,
    weighted_data_denominator = den_data, weighted_prior_denominator = den_prior,
    ratio_of_expectations = (num_data + num_prior) / (den_data + den_prior),
    pooled_beta_mean = (C + 1) / (S + 2),
    pooled_beta_lower = qbeta(.025, C + 1, S - C + 1),
    pooled_beta_upper = qbeta(.975, C + 1, S - C + 1),
    median_results_per_study = median(ct$m),
    median_significant_per_study = median(ct$y3 + ct$y4),
    studies_no_significant = sum(ct$y3 + ct$y4 == 0),
    studies_le4_results = sum(ct$m <= 4),
    prior_only_num_from_no_significance = sum(shrink[ct$y3 + ct$y4 == 0] * prior),
    key)
  say(key, sprintf("raw %d/%d = %.5f; posterior %.5f (SD %.5f)", C, S, C/S, mean(dr), sd(dr)))
  z
}

summary_file <- file.path(out, "posterior_comparison.csv")
rows <- list()
push_case <- function(...) {
  rows[[length(rows) + 1]] <<- summarize_case(...)
  write_csv(bind_rows(rows), summary_file)
}
for (threshold in c(1/3, 1)) {
  push_case(old_standard, "Original", "All", threshold, inclusive = FALSE)
  push_case(old_standard, "Original", "All", threshold, inclusive = TRUE)
  push_case(new_standard, "Revision", "All", threshold)
  for (phase_name in c("Early Phase I", "Phase I", "Phase II")) {
    push_case(filter(new_standard, phase == phase_name), "Revision", phase_name, threshold)
  }
}

# Verify independent implementation against the actual archived function and
# both saved versions of the current figure posterior (Monte Carlo tolerance).
res <- bind_rows(rows)
legacy_row <- filter(res, dataset == "Original", threshold == 1/3, p_rule == "< .05")
stopifnot(abs(legacy_row$posterior_mean - original_result$mean) <
  6 * sqrt(legacy_row$mcse_mean^2 + original_result$se^2 / 5000))
cache <- read_csv(file.path(root, "analysis/early_phase/outputs/posterior_summary.csv"), show_col_types = FALSE)
for (i in seq_len(nrow(cache))) {
  ph <- sub(" proportion$", "", cache$parameter[i])
  rr <- filter(res, dataset == "Revision", group == ph, threshold == 1)
  stopifnot(abs(rr$posterior_mean - cache$posterior_mean[i]) <
    6 * sqrt(rr$mcse_mean^2 + cache$posterior_sd[i]^2 / 20000))
}
say("Archived and current posterior replications pass independent Monte Carlo checks.")

# Prior sensitivity at the same threshold and fixed empirical counts.
for (prior in c(.1, .01)) {
  for (phase_name in c("Early Phase I", "Phase I", "Phase II")) {
    push_case(filter(new_standard, phase == phase_name), "Revision", phase_name, 1,
              prior = prior, draws = 20000L)
  }
}

# Separate cohort restriction from changes in p/eJAB and preprocessing.
bridge_sets <- list(
  "Original full cohort" = old_standard,
  "Original: interventional exact early phases, primary/secondary" = old_standard %>%
    filter(studyType == "INTERVENTIONAL", phase != "Other", outcomeType %in% c("PRIMARY", "SECONDARY")),
  "Original values on IDs shared with revision" = old_standard %>% filter(analysisId %in% new_standard$analysisId),
  "Revised values on IDs shared with original" = new_standard %>% filter(analysisId %in% old_standard$analysisId),
  "Revised full cohort" = new_standard
)
bridge <- bind_rows(lapply(names(bridge_sets), function(nm) {
  d <- bridge_sets[[nm]]
  bind_rows(lapply(c(1/3, 1), function(t) tibble(step = nm, threshold = t,
    analyses = nrow(d), studies = n_distinct(d$nctId), significant = sum(d$p <= .05),
    flagged = sum(d$p <= .05 & d$bf > t), empirical = flagged / significant)))
}))
write_csv(bridge, file.path(out, "cohort_bridge.csv"))

boundary <- bind_rows(lapply(list(Original = old_standard, Revision = new_standard), function(d) {
  tibble(at_p_05 = sum(d$p == .05), at_p_05_gt1 = sum(d$p == .05 & d$bf > 1),
    at_p_05_gtThird = sum(d$p == .05 & d$bf > 1/3),
    exactly_bf_1 = sum(d$bf == 1), exactly_bf_Third = sum(d$bf == 1/3))
}), .id = "dataset")
write_csv(boundary, file.path(out, "threshold_boundary_counts.csv"))
changed <- inner_join(old_standard %>% select(analysisId, nctId, p_old = p, bf_old = bf),
  new_standard %>% select(analysisId, p_new = p, bf_new = bf), by = "analysisId") %>%
  mutate(p_changed = abs(p_old - p_new) > 1e-12,
         bf_changed = abs(bf_old - bf_new) > 1e-10 * pmax(1, abs(bf_old)),
         old_sig = p_old <= .05, new_sig = p_new <= .05,
         old_gt1 = old_sig & bf_old > 1, new_gt1 = new_sig & bf_new > 1)
write_csv(filter(changed, p_changed | bf_changed), file.path(out, "changed_shared_records.csv"))
write_csv(tibble(only_original = sum(!old_standard$analysisId %in% new_standard$analysisId),
                 only_revision = sum(!new_standard$analysisId %in% old_standard$analysisId),
                 shared = nrow(changed), p_changed = sum(changed$p_changed),
                 bf_changed = sum(changed$bf_changed)), file.path(out, "record_overlap.csv"))
say("Cohort bridge:"); print(bridge, n = Inf)
say("Boundary checks:"); print(boundary)
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
say("Audit finished. No production files were changed.")
