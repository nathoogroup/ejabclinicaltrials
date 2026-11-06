# Descriptive, result-weighted contradiction proportions. Resample whole trials
# within each phase; do not treat repeated analyses as independent observations.
# These are bootstrap confidence intervals, not posterior credible intervals.

summarise_clinical_phases <- function(eligible, draws = 20000L, seed = 20261105L) {
  required <- c("phase", "nctId", "statistically_significant", "ejab_01")
  stopifnot(all(required %in% names(eligible)), draws >= 100L,
            !anyNA(eligible[required]), is.logical(eligible$statistically_significant))
  phases <- c("Early Phase I", "Phase I", "Phase II")
  stopifnot(all(eligible$phase %in% phases))
  memberships <- dplyr::distinct(eligible, nctId, phase)
  stopifnot(!anyDuplicated(memberships$nctId))
  if ("p_nhst" %in% names(eligible)) {
    stopifnot(all(eligible$statistically_significant == (eligible$p_nhst <= 0.05)))
  }

  studies <- eligible |>
    dplyr::group_by(phase, nctId) |>
    dplyr::summarise(
      eligible_analyses = dplyr::n(),
      significant_analyses = sum(statistically_significant),
      contradictions = sum(statistically_significant & ejab_01 > 1),
      above_one_third = sum(statistically_significant & ejab_01 > 1 / 3),
      .groups = "drop"
    ) |>
    dplyr::arrange(phase, nctId)

  # Keep the legacy posterior sampler's RNG stream unchanged when sourced.
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

  summaries <- bootstrap <- leave_one_out <- vector("list", length(phases))
  for (j in seq_along(phases)) {
    phase <- phases[[j]]
    x <- studies[as.character(studies$phase) == phase, ]
    k <- nrow(x)
    n_sig <- sum(x$significant_analyses)
    n_contradiction <- sum(x$contradictions)
    stopifnot(k > 1L, n_sig > 0L)
    values <- numeric(draws)
    for (b in seq_len(draws)) {
      selected <- sample.int(k, size = k, replace = TRUE)
      denominator <- sum(x$significant_analyses[selected])
      values[[b]] <- if (denominator > 0) sum(x$contradictions[selected]) / denominator else NA_real_
    }
    invalid <- sum(is.na(values))
    if (invalid > 0L) {
      stop(sprintf("%s: %d bootstrap replicates had no significant results; review the interval method", phase, invalid))
    }
    loo_denominator <- n_sig - x$significant_analyses
    loo <- ifelse(loo_denominator > 0, (n_contradiction - x$contradictions) / loo_denominator, NA_real_)
    stopifnot(!anyNA(loo))
    limits <- unname(stats::quantile(values, c(0.025, 0.975), type = 7))
    summaries[[j]] <- data.frame(
      phase, eligible_trials = k, eligible_analyses = sum(x$eligible_analyses),
      significant_analyses = n_sig, contradictions = n_contradiction,
      observed = n_contradiction / n_sig, lower = limits[[1]], upper = limits[[2]],
      bootstrap_replicates = draws, seed, zero_denominator_draws = invalid,
      trials_no_significant_results = sum(x$significant_analyses == 0L),
      trials_with_contradictions = sum(x$contradictions > 0L),
      leave_one_trial_out_min = min(loo), leave_one_trial_out_max = max(loo)
    )
    bootstrap[[j]] <- data.frame(phase, draw = seq_len(draws), proportion = values)
    leave_one_out[[j]] <- data.frame(
      phase, omitted_trial = x$nctId, omitted_significant = x$significant_analyses,
      omitted_contradictions = x$contradictions, proportion = loo,
      change = loo - n_contradiction / n_sig
    )
  }
  raw_summary <- studies |>
    dplyr::group_by(phase) |>
    dplyr::summarise(
      significant_analyses = sum(significant_analyses), contradictions = sum(contradictions),
      above_one_third = sum(above_one_third), .groups = "drop"
    )
  raw_summary <- dplyr::bind_rows(raw_summary, data.frame(
    phase = "All included phases", significant_analyses = sum(studies$significant_analyses),
    contradictions = sum(studies$contradictions), above_one_third = sum(studies$above_one_third)
  )) |>
    dplyr::mutate(
      contradiction_proportion = contradictions / significant_analyses,
      above_one_third_proportion = above_one_third / significant_analyses
    )
  list(summary = dplyr::bind_rows(summaries), trial_counts = studies,
       draws = dplyr::bind_rows(bootstrap), leave_one_out = dplyr::bind_rows(leave_one_out),
       raw_evidence_summary = raw_summary)
}

make_clinical_phase_plot <- function(summary, base_size = 12) {
  phases <- c("Early Phase I", "Phase I", "Phase II")
  summary$phase <- factor(summary$phase, levels = rev(phases))
  labels <- stats::setNames(sprintf(
    "%s\n%s/%s (%.2f%%)", summary$phase,
    format(summary$contradictions, big.mark = ",", trim = TRUE),
    format(summary$significant_analyses, big.mark = ",", trim = TRUE),
    100 * summary$observed
  ), as.character(summary$phase))
  ggplot2::ggplot(summary, ggplot2::aes(x = observed, y = phase, colour = phase)) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = lower, xmax = upper),
                         width = 0.14, linewidth = 0.8, orientation = "y") +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::scale_color_manual(values = c(
      "Early Phase I" = "#0072B2", "Phase I" = "#D55E00", "Phase II" = "#009E73"
    ), guide = "none") +
    ggplot2::scale_y_discrete(labels = labels) +
    ggplot2::scale_x_continuous(labels = scales::label_percent(accuracy = 1),
                              limits = c(0, max(summary$upper) * 1.08),
                              expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(
      x = "Proportion with eJAB > 1\namong NHST-significant results", y = NULL,
      subtitle = "Observed proportions and 95% intervals\nfrom resampling whole trials"
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   axis.title = ggplot2::element_text(size = 11),
                   axis.text = ggplot2::element_text(size = 10),
                   plot.subtitle = ggplot2::element_text(size = 11))
}

make_clinical_prior_audit_plot <- function(draws, observed, base_size = 11) {
  phases <- c("Early Phase I", "Phase I", "Phase II")
  draws$phase <- factor(draws$phase, levels = phases)
  observed$phase <- factor(observed$phase, levels = phases)
  ggplot2::ggplot(draws, ggplot2::aes(value, fill = phase, colour = phase)) +
    ggplot2::geom_density(alpha = 0.18, linewidth = 0.85) +
    ggplot2::geom_vline(data = observed, ggplot2::aes(xintercept = observed, colour = phase),
                        linetype = "dashed", linewidth = 0.7, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = c("#0072B2", "#D55E00", "#009E73"), name = NULL) +
    ggplot2::scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73"), name = NULL) +
    ggplot2::scale_x_continuous(breaks = seq(0, 0.4, by = 0.1)) +
    ggplot2::coord_cartesian(xlim = c(0, 0.4)) +
    ggplot2::labs(
      title = "Previous per-trial Dirichlet model",
      subtitle = "eJAB > 1; shown for prior-sensitivity audit",
      x = "Proportion with eJAB > 1 among NHST-significant results", y = "Density",
      caption = "Dashed lines: observed proportions; curves: previous posterior densities."
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(legend.position = "top", panel.grid.minor = ggplot2::element_blank(),
                   plot.caption = ggplot2::element_text(hjust = 0))
}

write_clinical_phase_summary <- function(result, output_dir, figure_dirs) {
  outputs <- c(summary = "clinical_phase_bootstrap_summary.csv",
               trial_counts = "clinical_phase_trial_counts.csv",
               draws = "clinical_phase_bootstrap_draws.csv",
               leave_one_out = "clinical_phase_leave_one_trial_out.csv",
               raw_evidence_summary = "clinical_phase_raw_evidence_summary.csv")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  for (name in names(outputs)) readr::write_csv(result[[name]], file.path(output_dir, outputs[[name]]))
  plot <- make_clinical_phase_plot(result$summary, base_size = 11)
  for (directory in unique(figure_dirs)) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(file.path(directory, "early_phase_observed.png"), plot,
                    width = 10.5, height = 3.5, units = "in", dpi = 400, bg = "white")
  }
  invisible(result)
}
