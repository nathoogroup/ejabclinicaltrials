# Supplementary sensitivity analysis for the ORIGINAL weighted, independent
# per-study Dirichlet model, using only the current eligible contradiction data.
# Small symmetric Dirichlet concentration parameters favor boundary/vertex
# configurations; this is robustness within a specified prior family, not a
# claim that a small concentration is generically weakly informative.
# Run from the repository root:
# Rscript analysis/early_phase/posterior_prior_sensitivity.R

prior_sensitivity_phases <- c("Early Phase I", "Phase I", "Phase II")
prior_sensitivity_priors <- c(1e-2, 1e-4, 1e-6)

prior_sensitivity_counts <- function(eligible) {
  required <- c("analysisId", "nctId", "phase", "p_nhst", "ejab_01")
  stopifnot(all(required %in% names(eligible)), !anyNA(eligible[required]),
            !anyDuplicated(eligible$analysisId),
            all(eligible$phase %in% prior_sensitivity_phases))
  membership <- dplyr::distinct(eligible, nctId, phase)
  stopifnot(!anyDuplicated(membership$nctId))
  eligible |>
    dplyr::mutate(sig = p_nhst <= .05, fav = ejab_01 > 1) |>
    dplyr::group_by(phase, nctId) |>
    dplyr::summarise(m = dplyr::n(),
      y1 = sum(!sig & !fav), y2 = sum(!sig & fav),
      y3 = sum(sig & !fav), y4 = sum(sig & fav), .groups = "drop") |>
    dplyr::arrange(phase, nctId)
}

prior_sensitivity_analytic <- function(counts, prior) {
  total <- counts$m + 4 * prior
  a4 <- counts$y4 + prior
  a34 <- counts$y3 + counts$y4 + 2 * prior
  numerator <- sum(counts$m * a4 / total)
  denominator <- sum(counts$m * a34 / total)
  variance_numerator <- sum(counts$m^2 * a4 * (total - a4) / (total^2 * (total + 1)))
  variance_denominator <- sum(counts$m^2 * a34 * (total - a34) / (total^2 * (total + 1)))
  covariance <- sum(counts$m^2 * a4 * (total - a34) / (total^2 * (total + 1)))
  data.frame(
    expected_numerator = numerator, expected_denominator = denominator,
    weighted_data_numerator = sum(counts$m * counts$y4 / total),
    weighted_prior_numerator = sum(counts$m * prior / total),
    weighted_data_denominator = sum(counts$m * (counts$y3 + counts$y4) / total),
    weighted_prior_denominator = sum(counts$m * 2 * prior / total),
    variance_numerator, variance_denominator, covariance_numerator_denominator = covariance,
    ratio_of_expectations = numerator / denominator,
    second_order_ratio_mean = numerator / denominator +
      numerator * variance_denominator / denominator^3 - covariance / denominator^2
  )
}

prior_sensitivity_draw_cell <- function(counts, prior, draws, seed, batches = 25L) {
  stopifnot(draws >= 50000L, draws %% batches == 0L)
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  numerator <- denominator <- numeric(draws)
  gamma_zeros <- zero_count_gamma_zeros <- nonzero_count_gamma_zeros <- 0
  zero_count_gamma_values <- 0
  min_row_sum <- Inf
  max_normalization_error <- 0
  # Retain all four categories and the same count weights as the original
  # implementation, rather than replacing this model with a pooled Beta law.
  for (i in seq_len(nrow(counts))) {
    observed <- unlist(counts[i, c("y1", "y2", "y3", "y4")], use.names = FALSE)
    shape <- observed + prior
    gamma <- matrix(stats::rgamma(draws * 4L, shape = rep(shape, each = draws)), ncol = 4L)
    zero_by_category <- colSums(gamma == 0)
    gamma_zeros <- gamma_zeros + sum(zero_by_category)
    zero_count_gamma_zeros <- zero_count_gamma_zeros + sum(zero_by_category[observed == 0])
    nonzero_count_gamma_zeros <- nonzero_count_gamma_zeros + sum(zero_by_category[observed > 0])
    zero_count_gamma_values <- zero_count_gamma_values + sum(observed == 0) * draws
    total <- rowSums(gamma)
    stopifnot(all(is.finite(total)), all(total > 0))
    min_row_sum <- min(min_row_sum, min(total))
    probabilities <- gamma / total
    max_normalization_error <- max(max_normalization_error, max(abs(rowSums(probabilities) - 1)))
    numerator <- numerator + counts$m[[i]] * probabilities[, 4]
    denominator <- denominator + counts$m[[i]] * (probabilities[, 3] + probabilities[, 4])
  }
  proportion <- numerator / denominator
  stopifnot(all(is.finite(proportion)), all(denominator > 0),
            all(proportion > 0 & proportion < 1), max_normalization_error < 1e-12,
            nonzero_count_gamma_zeros == 0)
  analytic <- prior_sensitivity_analytic(counts, prior)
  numerator_mean_mcse <- sqrt(analytic$variance_numerator / draws)
  denominator_mean_mcse <- sqrt(analytic$variance_denominator / draws)
  numerator_mean_z <- (mean(numerator) - analytic$expected_numerator) / numerator_mean_mcse
  denominator_mean_z <- (mean(denominator) - analytic$expected_denominator) / denominator_mean_mcse
  stopifnot(abs(numerator_mean_z) < 6, abs(denominator_mean_z) < 6)
  batch_ids <- rep(seq_len(batches), each = draws / batches)
  by_batch <- lapply(split(proportion, batch_ids), function(x) {
    c(mean = mean(x), lower = unname(stats::quantile(x, .025)),
      upper = unname(stats::quantile(x, .975)))
  })
  batch_statistics <- data.frame(batch = seq_len(batches), do.call(rbind, by_batch), row.names = NULL)
  ci <- unname(stats::quantile(proportion, c(.025, .975)))
  summary <- cbind(data.frame(
    phase = as.character(counts$phase[[1]]), prior_per_cell = prior, draws, seed,
    eligible_trials = nrow(counts), eligible_analyses = sum(counts$m),
    significant_analyses = sum(counts$y3 + counts$y4), contradictions = sum(counts$y4),
    empirical_proportion = sum(counts$y4) / sum(counts$y3 + counts$y4),
    posterior_mean = mean(proportion), posterior_median = stats::median(proportion),
    posterior_sd = stats::sd(proportion), lower = ci[[1]], upper = ci[[2]],
    mcse_mean = stats::sd(proportion) / sqrt(draws),
    batch_mcse_mean = stats::sd(batch_statistics$mean) / sqrt(batches),
    batch_mcse_lower = stats::sd(batch_statistics$lower) / sqrt(batches),
    batch_mcse_upper = stats::sd(batch_statistics$upper) / sqrt(batches),
    batch_count = batches, draws_per_batch = draws / batches,
    empirical_numerator_mean = mean(numerator), empirical_denominator_mean = mean(denominator),
    numerator_mean_mcse, denominator_mean_mcse, numerator_mean_z, denominator_mean_z
  ), analytic)
  numeric_checks <- data.frame(
    phase = as.character(counts$phase[[1]]), prior_per_cell = prior,
    gamma_values = nrow(counts) * draws * 4,
    gamma_zeros, zero_count_gamma_values, zero_count_gamma_zeros,
    zero_count_gamma_zero_fraction = zero_count_gamma_zeros / zero_count_gamma_values,
    nonzero_count_gamma_zeros, nonpositive_dirichlet_row_sums = 0L,
    minimum_dirichlet_row_sum = min_row_sum,
    maximum_normalization_error = max_normalization_error,
    minimum_phase_denominator = min(denominator), nonpositive_phase_denominators = 0L,
    nonfinite_proportion_draws = 0L, zero_proportion_draws = sum(proportion == 0)
  )
  batch_statistics$phase <- summary$phase
  batch_statistics$prior_per_cell <- prior
  list(summary = summary, numerical_checks = numeric_checks, batch_statistics = batch_statistics,
    draws = data.frame(phase = summary$phase, prior_per_cell = prior,
      draw = seq_len(draws), proportion))
}

prior_sensitivity_pairwise <- function(summary) {
  rows <- list()
  for (phase in prior_sensitivity_phases) {
    x <- summary[summary$phase == phase, ]
    pairs <- utils::combn(seq_len(nrow(x)), 2)
    for (k in seq_len(ncol(pairs))) {
      left <- x[pairs[1, k], ]; right <- x[pairs[2, k], ]
      rows[[length(rows) + 1L]] <- data.frame(
        phase, prior_from = left$prior_per_cell, prior_to = right$prior_per_cell,
        mean_change = right$posterior_mean - left$posterior_mean,
        mcse_mean_change = sqrt(left$mcse_mean^2 + right$mcse_mean^2),
        lower_change = right$lower - left$lower,
        mcse_lower_change = sqrt(left$batch_mcse_lower^2 + right$batch_mcse_lower^2),
        upper_change = right$upper - left$upper,
        mcse_upper_change = sqrt(left$batch_mcse_upper^2 + right$batch_mcse_upper^2),
        ratio_of_expectations_change = right$ratio_of_expectations - left$ratio_of_expectations,
        absolute_mean_change_percentage_points = 100 * abs(right$posterior_mean - left$posterior_mean)
      )
    }
  }
  dplyr::bind_rows(rows)
}

prior_sensitivity_audit_check <- function(summary, audit_file) {
  audit <- readr::read_csv(audit_file, show_col_types = FALSE)
  previous <- audit[audit$dataset == "Revision" & audit$threshold == 1 &
      audit$p_rule == "<= .05" & audit$prior_per_cell == .01 &
      audit$group %in% prior_sensitivity_phases, ]
  current <- summary[summary$prior_per_cell == .01, ]
  rows <- list()
  for (phase in prior_sensitivity_phases) {
    old <- previous[previous$group == phase, ]; new <- current[current$phase == phase, ]
    stopifnot(nrow(old) == 1L, nrow(new) == 1L, old$flagged == new$contradictions,
      old$significant == new$significant_analyses, old$analyses == new$eligible_analyses,
      old$studies == new$eligible_trials,
      abs(old$ratio_of_expectations - new$ratio_of_expectations) < 1e-12)
    combined <- sqrt(old$mcse_mean^2 + new$mcse_mean^2)
    z <- (new$posterior_mean - old$posterior_mean) / combined
    stopifnot(abs(z) < 6)
    rows[[phase]] <- data.frame(phase, prior_per_cell = .01,
      old_posterior_mean = old$posterior_mean, new_posterior_mean = new$posterior_mean,
      combined_mcse = combined, standardized_mean_difference = z,
      analytic_ratio_absolute_difference = abs(old$ratio_of_expectations - new$ratio_of_expectations),
      old_lower = old$lower, new_lower = new$lower, old_upper = old$upper, new_upper = new$upper,
      counts_match = TRUE, mean_within_six_mcse = TRUE)
  }
  dplyr::bind_rows(rows)
}

make_prior_sensitivity_plot <- function(draws, summary) {
  prior_labels <- c("a = 0.01", "a = 0.0001", "a = 0.000001")
  curves <- list()
  for (phase in prior_sensitivity_phases) {
    for (prior in prior_sensitivity_priors) {
      values <- draws$proportion[draws$phase == phase & draws$prior_per_cell == prior]
      bw <- stats::bw.nrd0(values)
      density <- stats::density(c(values, -values), from = max(0, min(values) - 5 * bw), to = max(values) + 5 * bw,
                                bw = bw, n = 2048)
      curves[[length(curves) + 1L]] <- data.frame(phase,
        prior = prior_labels[match(prior, prior_sensitivity_priors)],
        proportion = density$x, density = 2 * density$y)
    }
  }
  curves <- dplyr::bind_rows(curves)
  curves$prior <- factor(curves$prior, levels = prior_labels)
  curves$phase <- factor(curves$phase, levels = prior_sensitivity_phases)
  observed <- unique(summary[c("phase", "empirical_proportion", "contradictions", "significant_analyses")])
  observed$phase <- factor(observed$phase, levels = prior_sensitivity_phases)
  labels <- stats::setNames(sprintf("%s: %s/%s observed (%.2f%%)",
    observed$phase, observed$contradictions, format(observed$significant_analyses, big.mark = ",", trim = TRUE),
    100 * observed$empirical_proportion), as.character(observed$phase))
  ggplot2::ggplot(curves, ggplot2::aes(proportion, density, colour = prior, linetype = prior)) +
    ggplot2::geom_vline(data = observed, ggplot2::aes(xintercept = empirical_proportion),
                       inherit.aes = FALSE, colour = "grey35", linetype = "longdash", linewidth = .55) +
    ggplot2::geom_line(linewidth = .8) +
    ggplot2::facet_wrap(~phase, ncol = 1, scales = "free", labeller = ggplot2::as_labeller(labels)) +
    ggplot2::scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73"), name = NULL) +
    ggplot2::scale_linetype_manual(values = c("solid", "dashed", "dotted"), name = NULL) +
    ggplot2::scale_x_continuous(labels = scales::label_percent(accuracy = .1),
                              expand = ggplot2::expansion(mult = c(0, .01))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .05))) +
    ggplot2::labs(x = "Proportion with eJAB > 1 among NHST-significant results",
                  y = "Posterior density",
                  caption = "Independent per-study Dirichlet(a, a, a, a) priors; 50,000 posterior draws per curve.\nGrey dashed lines: empirical proportions. Panel axes differ; small a favors boundary configurations.") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "top", legend.key.width = grid::unit(.85, "cm"),
                   panel.grid.minor = ggplot2::element_blank(),
                   strip.background = ggplot2::element_rect(fill = "grey96"),
                   strip.text = ggplot2::element_text(hjust = 0, size = 10.5),
                   plot.caption = ggplot2::element_text(hjust = 0, size = 9),
                   panel.spacing = grid::unit(.7, "lines"))
}

run_prior_sensitivity <- function(root = getwd(), draws = 50000L, base_seed = 202611060L) {
  output_dir <- file.path(root, "analysis/early_phase/outputs")
  figure_dir <- file.path(root, "analysis/early_phase/figures")
  eligible_file <- file.path(output_dir, "eligible_early_phase_analyses.csv")
  eligible <- readr::read_csv(eligible_file, show_col_types = FALSE)
  counts <- prior_sensitivity_counts(eligible)
  stopifnot(nrow(eligible) == 16563L, nrow(counts) == 1790L,
            sum(counts$y4) == 493L, sum(counts$y3 + counts$y4) == 6474L)
  readr::write_csv(counts, file.path(output_dir, "prior_sensitivity_study_counts.csv"))
  results <- list()
  index <- 0L
  for (phase in prior_sensitivity_phases) {
    for (prior in prior_sensitivity_priors) {
      index <- index + 1L
      cat(sprintf("Simulating %s, a=%g, %s draws ...\n", phase, prior, draws))
      flush.console()
      result <- prior_sensitivity_draw_cell(counts[counts$phase == phase, ], prior,
                                            draws, seed = base_seed + index)
      results[[index]] <- result
      cat(sprintf("  Mean %.8f, 95%% CI [%.8f, %.8f], MCSE %.8f.\n",
        result$summary$posterior_mean, result$summary$lower, result$summary$upper,
        result$summary$mcse_mean))
      flush.console()
    }
  }
  components <- c("summary", "numerical_checks", "batch_statistics", "draws")
  combined <- stats::setNames(lapply(components, function(name) {
    dplyr::bind_rows(lapply(results, `[[`, name))
  }), components)
  combined$pairwise_changes <- prior_sensitivity_pairwise(combined$summary)
  combined$audit_crosscheck <- prior_sensitivity_audit_check(combined$summary,
    file.path(root, "analysis/audits/posterior_20261104/results/posterior_comparison.csv"))
  for (name in names(combined)) readr::write_csv(combined[[name]],
    file.path(output_dir, paste0("prior_sensitivity_", name, ".csv")))
  figure <- make_prior_sensitivity_plot(combined$draws, combined$summary)
  ggplot2::ggsave(file.path(figure_dir, "FigureS13_prior_sensitivity.png"), figure,
                  width = 7.2, height = 7.4, units = "in", dpi = 400, bg = "white")
  writeLines(c(
    "Original per-study model: independent pi_i | data ~ Dirichlet(y_i1+a,...,y_i4+a).",
    "Categories: nonsignificant/BF<=1, nonsignificant/BF>1, significant/BF<=1, significant/BF>1.",
    "Significance: p_nhst<=.05, with p_nhst=p for unrestricted tests and 2*p1 for one-sided tests.",
    "Target: rho = sum_i m_i*pi_i4 / sum_i m_i*(pi_i3+pi_i4); m_i is each trial's eligible result count.",
    "Each of nine cells uses 50,000 independent normalized-gamma Dirichlet draws and a distinct fixed seed.",
    "Posterior-mean MCSE=SD/sqrt(draws); endpoint MCSE=SD(25 batch quantiles)/sqrt(25).",
    "Analytic weighted numerator and denominator expectations/variances/covariance audit simulation.",
    "ratio_of_expectations is E[numerator]/E[denominator], NOT E[rho]; second_order_ratio_mean is a delta-method approximation.",
    "Small-a symmetric Dirichlet priors concentrate near boundaries/vertices; results assess sensitivity within this family, not generic weak-prior robustness.",
    "Gamma draws for zero-count cells may underflow to zero; row sums remain strictly positive, normalization is checked, all phase ratios remain finite and positive.",
    "The supplied zero-count underflow counts describe floating-point values near zero; they are not posterior atoms at zero.",
    "No eligibility, Bayes-factor, bootstrap, original posterior, V2, main-figure or manuscript files are modified."
  ), file.path(output_dir, "prior_sensitivity_method.txt"))
  writeLines(trimws(capture.output(sessionInfo()), which = "right"),
             file.path(output_dir, "prior_sensitivity_session_info.txt"))
  invisible(combined)
}

if (sys.nframe() == 0L) run_prior_sensitivity()
