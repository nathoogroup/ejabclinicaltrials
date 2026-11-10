# Conditional proportion of NHST-significant findings with inconclusive eJAB.
# Uses the existing study-specific Dirichlet/multinomial model, retaining
# eligible-analysis-count weights. Run from the repository root:
# Rscript analysis/early_phase/inconclusive_posterior.R

inconclusive_phases <- c("Early Phase I", "Phase I", "Phase II")
inconclusive_priors <- c(1e-2, 1e-4, 1e-6)

inconclusive_counts <- function(eligible, rule) {
  eligible |>
    dplyr::mutate(sig = .data[[rule]] <= .05,
      inconclusive = ejab_01 >= 1 / 3 & ejab_01 <= 3) |>
    dplyr::group_by(phase, nctId) |>
    dplyr::summarise(m = dplyr::n(),
      y1 = sum(!sig & !inconclusive), y2 = sum(!sig & inconclusive),
      y3 = sum(sig & !inconclusive), y4 = sum(sig & inconclusive),
      .groups = "drop") |>
    dplyr::arrange(phase, nctId)
}

inconclusive_moments <- function(counts, prior) {
  total <- counts$m + 4 * prior
  a4 <- counts$y4 + prior
  a34 <- counts$y3 + counts$y4 + 2 * prior
  data.frame(
    expected_numerator = sum(counts$m * a4 / total),
    expected_denominator = sum(counts$m * a34 / total),
    variance_numerator = sum(counts$m^2 * a4 * (total - a4) / (total^2 * (total + 1))),
    variance_denominator = sum(counts$m^2 * a34 * (total - a34) / (total^2 * (total + 1))),
    covariance = sum(counts$m^2 * a4 * (total - a34) / (total^2 * (total + 1)))
  )
}

inconclusive_draws <- function(counts, prior, rule, draws, seed, batches = 25L) {
  stopifnot(draws >= 50000L, draws %% batches == 0L,
    all(rowSums(as.matrix(counts[c("y1", "y2", "y3", "y4")])) == counts$m))
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  numerator <- denominator <- matrix(0, nrow = draws, ncol = length(inconclusive_phases),
    dimnames = list(NULL, inconclusive_phases))
  gamma_zeros <- nonzero_count_gamma_zeros <- 0
  min_row_sum <- Inf
  max_normalization_error <- 0
  for (i in seq_len(nrow(counts))) {
    observed <- unlist(counts[i, c("y1", "y2", "y3", "y4")], use.names = FALSE)
    gamma <- matrix(stats::rgamma(draws * 4L,
      shape = rep(observed + prior, each = draws)), ncol = 4L)
    zero_by_category <- colSums(gamma == 0)
    gamma_zeros <- gamma_zeros + sum(zero_by_category)
    nonzero_count_gamma_zeros <- nonzero_count_gamma_zeros + sum(zero_by_category[observed > 0])
    total <- rowSums(gamma)
    stopifnot(all(is.finite(total)), all(total > 0))
    min_row_sum <- min(min_row_sum, min(total))
    probabilities <- gamma / total
    max_normalization_error <- max(max_normalization_error,
      max(abs(rowSums(probabilities) - 1)))
    phase <- counts$phase[[i]]
    numerator[, phase] <- numerator[, phase] + counts$m[[i]] * probabilities[, 4]
    denominator[, phase] <- denominator[, phase] + counts$m[[i]] *
      (probabilities[, 3] + probabilities[, 4])
  }
  stopifnot(nonzero_count_gamma_zeros == 0, max_normalization_error < 1e-12)
  numerator <- cbind(All = rowSums(numerator), numerator)
  denominator <- cbind(All = rowSums(denominator), denominator)
  rows <- list()
  for (group in colnames(numerator)) {
    group_counts <- if (group == "All") counts else counts[counts$phase == group, ]
    values <- numerator[, group] / denominator[, group]
    stopifnot(all(is.finite(values)), all(values > 0 & values < 1))
    moments <- inconclusive_moments(group_counts, prior)
    numerator_mean_z <- (mean(numerator[, group]) - moments$expected_numerator) /
      sqrt(moments$variance_numerator / draws)
    denominator_mean_z <- (mean(denominator[, group]) - moments$expected_denominator) /
      sqrt(moments$variance_denominator / draws)
    stopifnot(abs(numerator_mean_z) < 6, abs(denominator_mean_z) < 6)
    batch <- split(values, rep(seq_len(batches), each = draws / batches))
    ci <- unname(stats::quantile(values, c(.025, .975)))
    significant <- sum(group_counts$y3 + group_counts$y4)
    inconclusive <- sum(group_counts$y4)
    raw <- inconclusive / significant
    ratio_of_expectations <- moments$expected_numerator / moments$expected_denominator
    rows[[group]] <- cbind(data.frame(
      significance_rule = rule, group, prior_per_cell = prior, draws, seed,
      eligible_trials = nrow(group_counts), eligible_analyses = sum(group_counts$m),
      significant_analyses = significant, significant_inconclusive = inconclusive,
      empirical_proportion = raw, posterior_mean = mean(values),
      posterior_sd = stats::sd(values), lower = ci[[1]], upper = ci[[2]],
      mcse_mean = stats::sd(values) / sqrt(draws),
      batch_mcse_mean = stats::sd(vapply(batch, mean, numeric(1))) / sqrt(batches),
      batch_mcse_sd = stats::sd(vapply(batch, stats::sd, numeric(1))) / sqrt(batches),
      mean_minus_empirical_percentage_points = 100 * (mean(values) - raw),
      numerator_mean_z, denominator_mean_z, ratio_of_expectations,
      second_order_ratio_mean = ratio_of_expectations +
        moments$expected_numerator * moments$variance_denominator / moments$expected_denominator^3 -
        moments$covariance / moments$expected_denominator^2
    ), moments)
  }
  list(summary = dplyr::bind_rows(rows),
    numerical_checks = data.frame(significance_rule = rule, prior_per_cell = prior,
      gamma_values = nrow(counts) * draws * 4, gamma_zeros,
      nonzero_count_gamma_zeros, minimum_dirichlet_row_sum = min_row_sum,
      maximum_normalization_error = max_normalization_error,
      minimum_pooled_denominator = min(denominator[, "All"]),
      nonfinite_proportion_draws = 0L))
}

run_inconclusive_posterior <- function(root = getwd(), draws = 50000L,
                                     base_seed = 202611100L, selected_draws = 200000L) {
  output_dir <- file.path(root, "analysis/early_phase/outputs")
  input <- file.path(output_dir, "eligible_early_phase_analyses.csv")
  initial_hash <- tools::md5sum(input)
  eligible <- readr::read_csv(input, show_col_types = FALSE)
  required <- c("analysisId", "nctId", "phase", "p_nhst", "p_two_sided", "ejab_01")
  stopifnot(all(required %in% names(eligible)), !anyNA(eligible[required]),
    !anyDuplicated(eligible$analysisId), all(eligible$phase %in% inconclusive_phases),
    nrow(eligible) == 16563L,
    !anyDuplicated(dplyr::distinct(eligible, nctId, phase)$nctId))
  target_counts <- inconclusive_counts(eligible, "p_two_sided")
  raw_target <- sum(target_counts$y4) / sum(target_counts$y3 + target_counts$y4)
  calibration <- dplyr::bind_rows(lapply(inconclusive_priors, function(prior) {
    moments <- inconclusive_moments(target_counts, prior)
    approx_mean <- moments$expected_numerator / moments$expected_denominator +
      moments$expected_numerator * moments$variance_denominator / moments$expected_denominator^3 -
      moments$covariance / moments$expected_denominator^2
    data.frame(prior_per_cell = prior, empirical_proportion = raw_target,
      second_order_ratio_mean = approx_mean,
      approximate_difference_percentage_points = 100 * (approx_mean - raw_target))
  }))
  # Choose concentration independently of Monte Carlo noise, using the
  # analytic second-order ratio-mean approximation; validate it by simulation.
  selected_prior <- max(calibration$prior_per_cell[
    abs(calibration$second_order_ratio_mean - raw_target) < .0001])
  stopifnot(is.finite(selected_prior))
  calibration$selected <- calibration$prior_per_cell == selected_prior
  all_counts <- results <- list()
  index <- 0L
  for (rule in c("p_two_sided", "p_nhst")) {
    counts <- inconclusive_counts(eligible, rule)
    stopifnot(nrow(counts) == 1790L, sum(counts$y4) == 1974L,
      sum(counts$y3 + counts$y4) == if (rule == "p_two_sided") 6476L else 6474L)
    all_counts[[rule]] <- dplyr::mutate(counts, significance_rule = rule, .before = 1)
    for (prior in inconclusive_priors) {
      index <- index + 1L
      cell_draws <- if (rule == "p_two_sided" && prior == selected_prior)
        max(draws, selected_draws) else draws
      cat(sprintf("Simulating %s, a=%g, %d draws ...\n", rule, prior, cell_draws))
      flush.console()
      results[[index]] <- inconclusive_draws(counts, prior, rule, cell_draws,
        seed = base_seed + index)
      pooled <- results[[index]]$summary[results[[index]]$summary$group == "All", ]
      cat(sprintf("  Pooled mean %.8f; SD %.8f; MCSE %.8f; raw %.8f.\n",
        pooled$posterior_mean, pooled$posterior_sd, pooled$mcse_mean,
        pooled$empirical_proportion))
      flush.console()
    }
  }
  summary <- dplyr::bind_rows(lapply(results, `[[`, "summary"))
  checks <- dplyr::bind_rows(lapply(results, `[[`, "numerical_checks"))
  selected <- summary[summary$significance_rule == "p_two_sided" &
    summary$group == "All" & summary$prior_per_cell == selected_prior, ]
  stopifnot(nrow(selected) == 1L,
    abs(selected$posterior_mean - selected$empirical_proportion) +
      3 * selected$mcse_mean < .0001)
  summary$selected_for_main_result <- summary$significance_rule == "p_two_sided" &
    summary$group == "All" & summary$prior_per_cell == selected_prior
  readr::write_csv(summary, file.path(output_dir, "inconclusive_posterior_summary.csv"))
  readr::write_csv(dplyr::bind_rows(all_counts),
    file.path(output_dir, "inconclusive_posterior_study_counts.csv"))
  readr::write_csv(checks, file.path(output_dir, "inconclusive_posterior_numerical_checks.csv"))
  readr::write_csv(calibration, file.path(output_dir, "inconclusive_posterior_prior_calibration.csv"))
  writeLines(c(
    "Target: proportion with 1/3 <= applicable eJAB <= 3 AMONG results with plotted p <= .05.",
    "Main rule p_two_sided: p for unrestricted tests; 2*min(p1,1-p1) for one-sided tests, matching the scatterplot inset.",
    "Directional sensitivity p_nhst: p for unrestricted tests; 2*p1 for one-sided tests.",
    "The full inset cohort is retained: 16,563 primary/secondary findings from 1,790 trials.",
    "Categories: nonsignificant/non-inconclusive, nonsignificant/inconclusive, significant/non-inconclusive, significant/inconclusive.",
    "For trial i: y_i | pi_i ~ Multinomial(m_i, pi_i), independently pi_i ~ Dirichlet(a,a,a,a).",
    "Thus pi_i | y_i ~ Dirichlet(y_i1+a,...,y_i4+a), independently across trials.",
    "This is the existing study-specific model, with no estimated shared concentration or hyperprior.",
    "The conditional target is rho = sum_i m_i*pi_i4 / sum_i m_i*(pi_i3+pi_i4).",
    "Weights m_i are each trial's eligible finding count, so rho is a finding-weighted conditional proportion, not an equal-trial mean.",
    "Trial-level vectors represent within-trial clustering; uncertainty is conditional on the observed trials and their fixed weights, not predictive uncertainty for new trials.",
    sprintf("Each prior/rule uses %s independent posterior draws, increased to %s for the selected main-rule prior. Phase summaries and the pooled ratio use the same trial draws.", format(draws, big.mark = ",", trim = TRUE), format(max(draws, selected_draws), big.mark = ",", trim = TRUE)),
    "Posterior SD measures precision; MCSE(mean)=posterior SD/sqrt(draws) measures numerical simulation error.",
    "Analytic numerator/denominator expectations, variances and covariance validate simulation; ratio_of_expectations is not the posterior mean of rho.",
    sprintf("Select largest a in {1e-2,1e-4,1e-6} with analytic second-order ratio mean within 0.0001 (0.01 percentage points) of raw: a=%g; validate |posterior mean - raw| + 3*MCSE < 0.0001.", selected_prior),
    "This closeness-based choice follows the requested empirical calibration; all three prior results are reported.",
    "Small-a symmetric Dirichlet priors concentrate near boundaries; their choice is not claimed to establish general weak informativeness.",
    "Gamma zeros occur only in empty categories from floating-point underflow; positive row sums and probability normalization are checked.",
    "No bootstrap is used; no raw input, JLP catalogue, existing clinical figure or posterior result is overwritten."
  ), file.path(output_dir, "inconclusive_posterior_method.txt"))
  writeLines(trimws(capture.output(sessionInfo()), which = "right"),
    file.path(output_dir, "inconclusive_posterior_session_info.txt"))
  stopifnot(identical(initial_hash, tools::md5sum(input)))
  invisible(summary)
}

if (sys.nframe() == 0L) run_inconclusive_posterior()
