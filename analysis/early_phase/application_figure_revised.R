#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(cowplot)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
response_dir <- normalizePath(file.path(script_dir, ".."))
revision_dir <- normalizePath(file.path(script_dir, "figures"))
output_dir <- file.path(script_dir, "outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sample_levels <- c("[2, 30)", "[30, 900)", "[900, 30K)",
                   "[30K, 900K)", "[900K, 6030K]")
phase_levels <- c("Early Phase I", "Phase I", "Phase II")

eligible_file <- file.path(output_dir, "eligible_early_phase_analyses.csv")
if (!file.exists(eligible_file)) {
  stop("Run early_phase_analysis.R before application_figure_revised.R")
}

early_phase_data <- read_csv(eligible_file, show_col_types = FALSE) %>%
  transmute(
    analysisId,
    nctId,
    outcomeTitle,
    outcomeType,
    exact_phase = phase,
    pValue = p_two_sided,
    p_left_censored,
    statistically_significant,
    statisticalMethod = test_class,
    N = sample_size,
    JAB = ejab_01
  )

scatter_data <- early_phase_data %>%
  mutate(
    logJAB = log(JAB),
    sample_group = case_when(
      N < 30 ~ "[2, 30)",
      N < 900 ~ "[30, 900)",
      N < 30000 ~ "[900, 30K)",
      N < 900000 ~ "[30K, 900K)",
      TRUE ~ "[900K, 6030K]"
    ),
    sample_group = factor(sample_group, levels = sample_levels)
  )

test_counts <- scatter_data %>%
  count(statisticalMethod, name = "count")
test_labels <- setNames(
  paste0(test_counts$statisticalMethod, " (", test_counts$count, ")"),
  test_counts$statisticalMethod
)

evidence_bands <- tibble(
  ymin = c(-6, log(1 / 3), log(3)),
  ymax = c(log(1 / 3), log(3), 3),
  fill = c("#A6C7E5", "#B4B4B4", "#FFB84D")
)

scatter_base <- ggplot(scatter_data, aes(pValue, logJAB)) +
  geom_rect(
    data = evidence_bands,
    aes(xmin = 0, xmax = 1, ymin = ymin, ymax = ymax, fill = fill),
    alpha = 0.18,
    inherit.aes = FALSE
  ) +
  geom_point(
    aes(color = statisticalMethod, shape = sample_group),
    size = 0.55,
    alpha = 0.68
  ) +
  geom_vline(xintercept = 0.05, linetype = "dashed", linewidth = 0.55) +
  geom_hline(yintercept = c(log(1 / 3), log(3)),
             linetype = "dashed", linewidth = 0.55) +
  scale_fill_identity() +
  scale_color_brewer(
    name = "Test",
    palette = "Paired",
    labels = test_labels
  ) +
  scale_shape_manual(name = "Sample size", values = c(17, 15, 16, 1, 8)) +
  guides(
    color = guide_legend(override.aes = list(size = 2.3, alpha = 1), ncol = 1),
    shape = guide_legend(override.aes = list(size = 2.3), ncol = 1)
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.key.height = unit(0.36, "cm"),
    legend.text = element_text(size = 8.5),
    legend.title = element_text(size = 9.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

scatter_main <- scatter_base +
  coord_cartesian(xlim = c(0, 1), ylim = c(-6, 3), expand = FALSE) +
  scale_x_continuous(
    name = expression(italic(p) * "-value"),
    breaks = seq(0, 1, 0.25),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = expression("ln " * italic(eJAB)["01"]),
    breaks = seq(-6, 3, 3),
    expand = c(0, 0)
  )

scatter_inset <- scatter_base +
  coord_cartesian(xlim = c(0, 0.05), ylim = c(-3, 3), expand = FALSE) +
  scale_x_continuous(
    name = NULL,
    breaks = c(0, 0.025, 0.05),
    labels = c("0", ".025", ".05"),
    expand = c(0, 0)
  ) +
  scale_y_continuous(name = NULL, breaks = c(-3, 0, 3), expand = c(0, 0)) +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text = element_text(size = 7),
    axis.ticks = element_line(linewidth = 0.35),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = "black", linewidth = 0.7)
  )

scatter_plot <- ggdraw(scatter_main) +
  draw_plot(scatter_inset, x = 0.20, y = 0.10, width = 0.55, height = 0.60)

phase_data <- early_phase_data %>%
  mutate(
    phase = exact_phase,
    phase = factor(phase, levels = phase_levels)
  )

alpha_grid <- seq(0.001, 0.25, by = 0.001)
alpha_profile <- bind_rows(lapply(phase_levels, function(phase_name) {
  phase_frame <- filter(phase_data, phase == phase_name)
  tibble(
    phase = phase_name,
    alpha = alpha_grid,
    contradiction_rate = vapply(
      alpha_grid,
      function(alpha_level) {
        is_significant <- if_else(
          phase_frame$p_left_censored,
          phase_frame$pValue <= alpha_level,
          phase_frame$pValue < alpha_level
        )
        mean(is_significant & phase_frame$JAB > 1)
      },
      numeric(1)
    )
  )
})) %>%
  mutate(phase = factor(phase, levels = phase_levels))

alpha_plot <- ggplot(alpha_profile, aes(alpha, contradiction_rate, color = phase)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", linewidth = 0.6) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73")) +
  scale_x_continuous(limits = c(0, 0.25), breaks = seq(0, 0.25, 0.05), expand = expansion(mult = c(0, 0.015))) +
  scale_y_continuous(limits = c(0, 0.25), breaks = seq(0, 0.25, 0.05), expand = expansion(mult = c(0, 0.015))) +
  labs(
    x = expression("Significance level " * alpha),
    y = "Contradiction rate",
    color = "Trial phase"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

draw_phase_posterior <- function(phase_frame, n_sim = 5000, bf_threshold = 1) {
  study_counts <- phase_frame %>%
    group_by(nctId) %>%
    summarise(
      weight = n(),
      y1 = sum(!statistically_significant & JAB <= bf_threshold),
      y2 = sum(!statistically_significant & JAB > bf_threshold),
      y3 = sum(statistically_significant & JAB <= bf_threshold),
      y4 = sum(statistically_significant & JAB > bf_threshold),
      .groups = "drop"
    )

  numerator <- numeric(n_sim)
  denominator <- numeric(n_sim)
  for (index in seq_len(nrow(study_counts))) {
    alpha <- 1 + unlist(study_counts[index, c("y1", "y2", "y3", "y4")],
                        use.names = FALSE)
    draws <- vapply(alpha, function(shape) rgamma(n_sim, shape = shape),
                    numeric(n_sim))
    draws <- draws / rowSums(draws)
    weight <- study_counts$weight[index]
    numerator <- numerator + weight * draws[, 4]
    denominator <- denominator + weight * (draws[, 3] + draws[, 4])
  }

  tibble(draw = seq_len(n_sim), value = numerator / denominator)
}

set.seed(20260828)
posterior_long <- bind_rows(lapply(phase_levels, function(phase_name) {
  draw_phase_posterior(filter(phase_data, phase == phase_name)) %>%
    mutate(phase = phase_name)
})) %>%
  mutate(phase = factor(phase, levels = phase_levels))

posterior_summary <- posterior_long %>%
  group_by(phase) %>%
  summarise(
    mean = mean(value),
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975),
    .groups = "drop"
  ) %>%
  left_join(
    phase_data %>%
      group_by(phase) %>%
      summarise(analyses = n(), studies = n_distinct(nctId), .groups = "drop"),
    by = "phase"
  )

phase_labels <- posterior_summary %>%
  transmute(
    phase,
    label = sprintf("%s\nmean %.3f (%.3f, %.3f)", phase, mean, lower, upper)
  )
phase_labels <- setNames(phase_labels$label, phase_labels$phase)

write_csv(
  posterior_long,
  file.path(output_dir, "early_phase_contradiction_posterior_draws.csv")
)
write_csv(
  posterior_summary,
  file.path(output_dir, "early_phase_contradiction_posterior_summary.csv")
)

set.seed(20260903)
weak_evidence_posterior <- draw_phase_posterior(
  scatter_data,
  n_sim = 20000,
  bf_threshold = 1 / 3
)

application_summary <- tibble(
  metric = c(
    "eligible_results",
    "eligible_trials",
    "statistical_test_classes",
    "results_p_le_0_05",
    "contradictions_p_le_0_05_ejab_gt_1",
    "contradictions_p_le_0_01_ejab_gt_1",
    "jlp_p_le_0_05_ejab_gt_3",
    "weak_evidence_posterior_mean",
    "weak_evidence_posterior_sd"
  ),
  value = c(
    nrow(scatter_data),
    n_distinct(scatter_data$nctId),
    n_distinct(scatter_data$statisticalMethod),
    sum(scatter_data$statistically_significant),
    sum(scatter_data$statistically_significant & scatter_data$JAB > 1),
    sum(if_else(scatter_data$p_left_censored, scatter_data$pValue <= 0.01,
                scatter_data$pValue < 0.01) & scatter_data$JAB > 1),
    sum(scatter_data$statistically_significant & scatter_data$JAB > 3),
    mean(weak_evidence_posterior$value),
    sd(weak_evidence_posterior$value)
  )
)

write_csv(
  application_summary,
  file.path(output_dir, "early_phase_application_summary.csv")
)
write_csv(
  scatter_data %>%
    filter(statistically_significant, JAB > 1) %>%
    transmute(
      analysisId,
      nctId,
      phase = exact_phase,
      outcomeTitle,
      outcomeType,
      statisticalMethod,
      sample_size = N,
      p_value = pValue,
      ejab_01 = JAB,
      reported_as_inequality = p_left_censored
    ),
  file.path(output_dir, "early_phase_candidate_contradictions.csv")
)

posterior_plot <- ggplot(
  posterior_long,
  aes(value, fill = phase, color = phase)
) +
  geom_density(alpha = 0.30, linewidth = 0.75) +
  geom_vline(
    data = posterior_summary,
    aes(xintercept = mean, color = phase),
    linetype = "dashed",
    linewidth = 0.6,
    inherit.aes = FALSE
  ) +
  facet_wrap(
    ~phase,
    nrow = 1,
    scales = "free_y",
    labeller = as_labeller(phase_labels)
  ) +
  scale_fill_manual(values = c("#56B4E9", "#0072B2", "#009E73")) +
  scale_color_manual(values = c("#0072B2", "#004C80", "#007A55")) +
  scale_x_continuous(limits = c(0, 0.4), breaks = c(0, 0.1, 0.2, 0.3, 0.4)) +
  labs(x = "Posterior contradiction proportion", y = "Density") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(size = 10.5),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9.5)
  )

application_figure <- plot_grid(
  alpha_plot,
  scatter_plot,
  posterior_plot,
  labels = c("A", "B", "C"),
  label_size = 16,
  ncol = 1,
  rel_heights = c(1.0, 1.8, 1.0)
)

ggsave(
  file.path(revision_dir, "Figure3.png"),
  application_figure,
  width = 14,
  height = 14,
  dpi = 300,
  bg = "white"
)

print(application_summary)
print(posterior_summary)
