# V3 screened Jeffreys-Lindley discordance panels. Existing V2 artifacts are
# deliberately untouched. Reported p-value bounds remain screening bounds.

jlp_phase_levels <- c("Early Phase I", "Phase I", "Phase II")
jlp_phase_colors <- stats::setNames(c("#0072B2", "#D55E00", "#009E73"), jlp_phase_levels)

compute_jlp_phase_panels <- function(eligible, contradiction_draws = NULL,
                                     contradiction_summary = NULL,
                                     draws = 20000L, seed = 20261106L) {
  required <- c("analysisId", "nctId", "phase", "p_nhst", "ejab_01", "p_left_censored")
  stopifnot(all(required %in% names(eligible)), !anyNA(eligible[required]),
            !anyDuplicated(eligible$analysisId), all(eligible$phase %in% jlp_phase_levels))
  membership <- dplyr::distinct(eligible, nctId, phase)
  stopifnot(!anyDuplicated(membership$nctId))
  if ("statistically_significant" %in% names(eligible)) {
    stopifnot(all(eligible$statistically_significant == (eligible$p_nhst <= 0.05)))
  }
  data <- eligible |>
    dplyr::mutate(significant = p_nhst <= 0.05,
                  screened = significant & ejab_01 > 3)
  trials <- data |>
    dplyr::group_by(phase, nctId) |>
    dplyr::summarise(
      eligible_analyses = dplyr::n(), significant_analyses = sum(significant),
      screened_jlp = sum(screened), numeric_screens = sum(screened & !p_left_censored),
      bound_screens = sum(screened & p_left_censored), .groups = "drop"
    ) |>
    dplyr::arrange(phase, nctId)

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

  summaries <- bootstrap <- list()
  threshold_grid <- sort(unique(c(
    exp(seq(log(1 / 300), log(300), length.out = 1000)),
    seq(3, 5, length.out = 1001), 1 / 3, 1, 3,
    data$ejab_01[data$significant & data$ejab_01 >= 1 / 300 & data$ejab_01 <= 300]
  )))
  sweeps <- list()
  for (phase in jlp_phase_levels) {
    x <- trials[as.character(trials$phase) == phase, ]
    d <- data[as.character(data$phase) == phase, ]
    k <- nrow(x)
    s <- sum(x$significant_analyses)
    j <- sum(x$screened_jlp)
    stopifnot(k > 1L, s > 0L)
    values <- numeric(draws)
    for (b in seq_len(draws)) {
      selected <- sample.int(k, k, replace = TRUE)
      denominator <- sum(x$significant_analyses[selected])
      values[[b]] <- if (denominator > 0) sum(x$screened_jlp[selected]) / denominator else NA_real_
    }
    if (anyNA(values)) stop(phase, ": bootstrap sample with no significant results")
    limits <- unname(stats::quantile(values, c(0.025, 0.975), type = 7))
    n_positive_trials <- sum(x$screened_jlp > 0L)
    summaries[[phase]] <- data.frame(
      phase, eligible_trials = k, eligible_analyses = nrow(d), significant_analyses = s,
      screened_jlp = j, numeric_screens = sum(x$numeric_screens), bound_screens = sum(x$bound_screens),
      observed_overall = j / nrow(d), observed_among_significant = j / s,
      lower = limits[[1]], upper = limits[[2]], zero_mass = mean(values == 0),
      zero_mass_theoretical = ((k - n_positive_trials) / k)^k,
      zero_replicates = sum(values == 0), positive_screen_trials = n_positive_trials,
      bootstrap_replicates = draws, seed, zero_denominator_draws = 0L
    )
    bootstrap[[phase]] <- data.frame(phase, draw = seq_len(draws), proportion = values)
    counts <- vapply(threshold_grid, function(threshold) {
      sum(d$significant & d$ejab_01 > threshold)
    }, integer(1))
    stopifnot(all(diff(counts) <= 0), counts[threshold_grid == 3] == j)
    sweeps[[phase]] <- data.frame(
      phase, bf_threshold = threshold_grid, screened_results = counts,
      eligible_analyses = nrow(d), proportion = counts / nrow(d)
    )
  }
  summary <- dplyr::bind_rows(summaries)
  bootstrap <- dplyr::bind_rows(bootstrap)

  result <- list(summary = summary, trials = trials, draws = bootstrap,
                 threshold_sweep = dplyr::bind_rows(sweeps))
  density <- phase_positive_kde(bootstrap, summary)
  result$kde <- density$kde
  result$kde_checks <- density$checks
  if (!is.null(contradiction_draws)) {
    result <- add_contradiction_density(result, contradiction_draws, contradiction_summary)
  }
  result
}

# A mixed empirical bootstrap law: a genuine point mass at zero and a
# separately smoothed positive component. Gaussian reflection at zero avoids
# negative support. Its area is 1 - the empirical point mass, not forcibly 1.
phase_positive_kde <- function(draws, summary) {
  required <- c("phase", "draw", "proportion")
  stopifnot(all(required %in% names(draws)), !anyNA(draws[required]),
            !anyDuplicated(draws[c("phase", "draw")]),
            all(draws$phase %in% jlp_phase_levels),
            all(draws$proportion >= 0 & draws$proportion <= 1))
  positive <- split(draws$proportion[draws$proportion > 0],
                    draws$phase[draws$proportion > 0])
  positive <- positive[lengths(positive) > 1L]
  if (!length(positive)) stop("No nondegenerate positive bootstrap component")
  bandwidth <- vapply(positive, stats::bw.nrd0, numeric(1))
  stopifnot(all(is.finite(bandwidth)), all(bandwidth > 0))
  xmax <- max(vapply(names(positive), function(phase) {
    max(positive[[phase]]) + 5 * bandwidth[[phase]]
  }, numeric(1)))
  curves <- checks <- list()
  for (phase in jlp_phase_levels) {
    phase_draws <- draws$proportion[draws$phase == phase]
    stopifnot(length(phase_draws) > 1L)
    zero_mass <- mean(phase_draws == 0)
    weight <- 1 - zero_mass
    if (!phase %in% names(positive)) {
      stopifnot(weight == 0)
      checks[[phase]] <- data.frame(phase, positive_draws = 0L,
        bandwidth = NA_real_, zero_mass, continuous_weight = 0,
        integrated_density = 0, integration_error = 0)
      next
    }
    kd <- stats::density(c(positive[[phase]], -positive[[phase]]),
                         bw = bandwidth[[phase]], kernel = "gaussian",
                         from = 0, to = xmax, n = 4096)
    y <- 2 * weight * kd$y
    area <- sum(diff(kd$x) * (head(y, -1) + tail(y, -1)) / 2)
    stopifnot(abs(area - weight) < 0.005)
    curves[[phase]] <- data.frame(phase, proportion = kd$x, density = y)
    checks[[phase]] <- data.frame(phase, positive_draws = length(positive[[phase]]),
      bandwidth = bandwidth[[phase]], zero_mass, continuous_weight = weight,
      integrated_density = area, integration_error = area - weight)
  }
  list(kde = dplyr::bind_rows(curves), checks = dplyr::bind_rows(checks))
}

add_contradiction_density <- function(result, draws, summary = NULL) {
  density <- phase_positive_kde(draws, summary)
  if (is.null(summary)) {
    # The exact empirical counts are also available at K = 1 in the sweep.
    summary <- result$threshold_sweep[result$threshold_sweep$bf_threshold == 1, ]
    names(summary)[names(summary) == "screened_results"] <- "contradictions"
    summary <- dplyr::left_join(summary,
      result$summary[c("phase", "significant_analyses")], by = "phase")
    summary$observed <- summary$contradictions / summary$significant_analyses
  }
  stopifnot(all(c("phase", "observed", "contradictions", "significant_analyses") %in% names(summary)))
  for (phase in jlp_phase_levels) {
    x <- summary[summary$phase == phase, ]
    k1 <- result$threshold_sweep[result$threshold_sweep$phase == phase &
                                 result$threshold_sweep$bf_threshold == 1, ]
    stopifnot(nrow(x) == 1L, x$contradictions == k1$screened_results,
              abs(x$observed - x$contradictions / x$significant_analyses) < 1e-12)
    # Original intervals must continue to agree with the supplied exact draws.
    if (all(c("lower", "upper") %in% names(x))) {
      ci <- unname(stats::quantile(draws$proportion[draws$phase == phase], c(.025, .975)))
      stopifnot(max(abs(ci - c(x$lower, x$upper))) < 1e-12)
    }
  }
  summary$zero_mass <- density$checks$zero_mass[match(summary$phase, density$checks$phase)]
  result$contradiction_draws <- draws
  result$contradiction_summary <- summary
  result$contradiction_kde <- density$kde
  result$contradiction_kde_checks <- density$checks
  result
}

make_jlp_threshold_plot <- function(sweep) {
  sweep$phase <- factor(sweep$phase, levels = jlp_phase_levels)
  base <- ggplot2::ggplot(sweep, ggplot2::aes(bf_threshold, proportion, colour = phase)) +
    ggplot2::geom_vline(xintercept = c(1, 3), linetype = "dashed",
                        linewidth = 0.5, colour = "grey45") +
    ggplot2::geom_step(linewidth = 0.85, direction = "hv") +
    ggplot2::annotate("text", x = 1, y = max(sweep$proportion) * 0.31,
                      label = "K = 1", size = 3.1, hjust = 1.15, colour = "grey25") +
    ggplot2::annotate("text", x = 3, y = max(sweep$proportion) * 0.31,
                      label = "K = 3", size = 3.1, hjust = -0.15, colour = "grey25") +
    ggplot2::scale_color_manual(values = jlp_phase_colors, drop = FALSE, name = NULL) +
    ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 1),
                              expand = ggplot2::expansion(mult = c(0, 0.04))) +
    ggplot2::scale_x_log10(breaks = c(1 / 300, 1 / 30, 1 / 3, 1, 3, 30, 300),
                         labels = c("1/300", "1/30", "1/3", "1", "3", "30", "300")) +
    ggplot2::labs(x = "Bayes-factor threshold K", y = "Proportion of all eligible results",
                  subtitle = "NHST significant and eJAB > K") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "top", legend.text = ggplot2::element_text(size = 9.5),
                   legend.key.width = grid::unit(0.42, "cm"),
                   panel.grid.minor = ggplot2::element_blank(),
                   plot.subtitle = ggplot2::element_text(size = 11),
                   axis.title = ggplot2::element_text(size = 11),
                   axis.text = ggplot2::element_text(size = 9))
  zoom <- sweep[sweep$bf_threshold >= 3 & sweep$bf_threshold <= 5, ]
  inset <- ggplot2::ggplot(zoom, ggplot2::aes(bf_threshold, proportion, colour = phase)) +
    ggplot2::geom_step(linewidth = 0.7, direction = "hv") +
    ggplot2::scale_color_manual(values = jlp_phase_colors, guide = "none") +
    ggplot2::scale_x_continuous(breaks = c(3, 4, 5)) +
    ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 0.01),
                              limits = c(0, max(zoom$proportion) * 1.08)) +
    ggplot2::labs(x = "K", y = NULL, title = "Screened JLP range") +
    ggplot2::theme_bw(base_size = 8) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 8.5),
                   plot.background = ggplot2::element_rect(fill = "white", colour = "black"),
                   plot.margin = ggplot2::margin(4, 5, 3, 4))
  cowplot::ggdraw(base) + cowplot::draw_plot(inset, x = 0.56, y = 0.42, width = 0.42, height = 0.42)
}

make_zero_mass_inset <- function(summary) {
  summary$phase <- factor(summary$phase, levels = jlp_phase_levels)
  summary$zero_label <- vapply(summary$zero_mass, function(x) {
    if (x == 0) "0%" else if (x == 1) "100%" else sprintf("%.2f%%", 100 * x)
  }, character(1))
  ggplot2::ggplot(summary, ggplot2::aes(phase, 0.6, colour = phase)) +
    ggplot2::geom_point(shape = 18, size = 2, show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(y = 0.25, label = zero_label), size = 2.7,
                       fontface = "bold", show.legend = FALSE) +
    ggplot2::scale_color_manual(values = jlp_phase_colors, guide = "none") +
    ggplot2::scale_x_discrete(labels = c("Early I", "I", "II")) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL, title = "Bootstrap mass at zero") +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size = 8),
                   plot.title = ggplot2::element_text(size = 8.5, hjust = 0.5),
                   plot.background = ggplot2::element_rect(fill = "white", colour = "grey65"),
                   plot.margin = ggplot2::margin(3, 3, 2, 3))
}

make_phase_density_row <- function(curves, summary, title, x_accuracy = 1,
                                   observed_column = "observed") {
  curves$phase <- factor(curves$phase, levels = jlp_phase_levels)
  summary$phase <- factor(summary$phase, levels = jlp_phase_levels)
  summary$empirical <- summary[[observed_column]]
  point_masses <- summary[summary$zero_mass == 1, ]
  base <- ggplot2::ggplot(curves, ggplot2::aes(proportion, density, colour = phase)) +
    ggplot2::geom_vline(data = summary, ggplot2::aes(xintercept = empirical, colour = phase),
                        linetype = "dashed", linewidth = 0.45, alpha = 0.8,
                        show.legend = FALSE) +
    ggplot2::geom_line(linewidth = 0.85) +
    ggplot2::geom_point(data = point_masses,
                         ggplot2::aes(x = 0, y = 0, colour = phase),
                         inherit.aes = FALSE, shape = 18, size = 2.8, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = jlp_phase_colors, drop = FALSE,
                               breaks = jlp_phase_levels, name = NULL) +
    ggplot2::scale_x_continuous(labels = scales::label_percent(accuracy = x_accuracy),
                              breaks = scales::breaks_width(if (x_accuracy >= 1) 0.05 else 0.0025),
                              limits = c(0, max(curves$proportion)),
                              expand = ggplot2::expansion(mult = c(0.005, 0.01))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.06))) +
    ggplot2::labs(x = "Proportion among NHST-significant results", y = "Bootstrap density",
                  subtitle = title) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "none", panel.grid.minor = ggplot2::element_blank(),
                   plot.subtitle = ggplot2::element_text(size = 10.5),
                   axis.title = ggplot2::element_text(size = 9.5),
                   axis.text = ggplot2::element_text(size = 8.5),
                   plot.margin = ggplot2::margin(4, 5, 3, 2))
  cowplot::ggdraw(base) + cowplot::draw_plot(make_zero_mass_inset(summary),
     x = 0.66, y = 0.49, width = 0.31, height = 0.31)
}

make_combined_phase_density_plot <- function(result, contradiction_draws = NULL,
                                             contradiction_summary = NULL) {
  if (!is.null(contradiction_draws)) {
    result <- add_contradiction_density(result, contradiction_draws, contradiction_summary)
  }
  if (is.null(result$contradiction_kde)) {
    stop("Supply existing contradiction draws to compute_jlp_phase_panels or this function")
  }
  legend_plot <- ggplot2::ggplot(result$contradiction_kde,
     ggplot2::aes(proportion, density, colour = factor(phase, levels = jlp_phase_levels))) +
    ggplot2::geom_line() +
    ggplot2::scale_color_manual(values = jlp_phase_colors, name = NULL, drop = FALSE) +
    ggplot2::theme_void() + ggplot2::theme(legend.position = "top",
      legend.text = ggplot2::element_text(size = 9.5),
      legend.key.width = grid::unit(0.5, "cm"), legend.margin = ggplot2::margin(0, 0, 0, 0))
  legend <- cowplot::get_legend(legend_plot)
  top <- make_phase_density_row(result$contradiction_kde, result$contradiction_summary,
    "Contradictions: eJAB > 1", x_accuracy = 1)
  bottom <- make_phase_density_row(result$kde, result$summary,
    "Screened JLPs: eJAB > 3", x_accuracy = 0.01,
    observed_column = "observed_among_significant")
  cowplot::plot_grid(legend, top, bottom, ncol = 1, rel_heights = c(0.18, 1, 1))
}

# Retained for callers that need only the stronger-evidence screen.
make_jlp_bootstrap_density_plot <- function(result) {
  make_phase_density_row(result$kde, result$summary,
    "Screened JLPs: eJAB > 3", x_accuracy = 0.01,
    observed_column = "observed_among_significant")
}

write_jlp_phase_panels <- function(result, output_dir, figure_dirs) {
  outputs <- c(summary = "jlp_phase_bootstrap_summary.csv", trials = "jlp_phase_trial_counts.csv",
             draws = "jlp_phase_bootstrap_draws.csv", threshold_sweep = "jlp_phase_threshold_sweep.csv",
             kde = "jlp_phase_bootstrap_kde.csv", kde_checks = "jlp_phase_kde_checks.csv")
  if (!is.null(result$contradiction_kde)) {
    outputs <- c(outputs, contradiction_kde = "clinical_phase_bootstrap_kde.csv",
                 contradiction_kde_checks = "clinical_phase_kde_checks.csv")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  for (name in names(outputs)) readr::write_csv(result[[name]], file.path(output_dir, outputs[[name]]))
  plots <- list(jlp_phase_thresholds = make_jlp_threshold_plot(result$threshold_sweep))
  if (!is.null(result$contradiction_kde)) {
    plots$phase_contradiction_jlp_density <- make_combined_phase_density_plot(result)
  }
  for (directory in unique(figure_dirs)) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    for (name in names(plots)) ggplot2::ggsave(file.path(directory, paste0(name, ".png")),
      plots[[name]], width = 6.6, height = 4.4, units = "in", dpi = 400, bg = "white")
  }
  invisible(plots)
}
