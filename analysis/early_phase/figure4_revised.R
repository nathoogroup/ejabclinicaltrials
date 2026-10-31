#!/usr/bin/env Rscript

# Build the outcome-distribution panels used in manuscript Figure 3D-G.
# Sourcing this file provides make_panel() without writing a figure. Running it
# directly also reproduces the legacy Figure 4 layout for comparison.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(lme4)
  library(readr)
  library(tidyr)
})

common_y_limits <- c(-18, 8)

make_panel <- function(data, outcome = "PRIMARY", adjusted = FALSE, compact = FALSE) {
  title_size <- if (compact) 12 else 11
  axis_size <- if (compact) 10 else 9
  plot_data <- data %>%
    filter(toupper(outcomeType) == outcome) %>%
    mutate(
      Phase = phase,
      Phase = factor(Phase, levels = c("Early Phase I", "Phase I", "Phase II")),
      logP = log(p_two_sided),
      logJAB = log(ejab_01),
      N = sample_size,
      Study = nctId
    ) %>%
    pivot_longer(cols = c(logJAB, logP), names_to = "Metric", values_to = "Value")

  average_n <- plot_data %>%
    group_by(outcomeType, Phase) %>%
    summarise(avg_N = mean(N, na.rm = TRUE), .groups = "drop") %>%
    mutate(label = paste0("italic(n) == ", round(avg_N)))

  if (adjusted) {
    adjust_metric <- function(metric_data) {
      metric_data <- filter(metric_data, !is.na(Study), !is.na(Value))
      model <- lmer(Value ~ 1 + (1 | Study), data = metric_data, REML = TRUE)
      metric_data$Value <- resid(model) + fixef(model)[["(Intercept)"]]
      metric_data
    }
    plot_data <- plot_data %>%
      group_by(Metric) %>%
      group_modify(~ adjust_metric(.x)) %>%
      ungroup()
  }

  ggplot(plot_data, aes(x = Phase, y = Value, fill = Metric)) +
    geom_boxplot(alpha = 0.7) +
    geom_hline(yintercept = log(1 / 3), linetype = "dashed") +
    geom_hline(yintercept = log(3), linetype = "dashed") +
    geom_hline(yintercept = log(0.05), linetype = "dashed", color = "red") +
    geom_text(
      data = average_n,
      aes(x = Phase, y = 7.1, label = label),
      parse = TRUE,
      inherit.aes = FALSE,
      size = if (compact) 3.6 else 3.4
    ) +
    coord_cartesian(ylim = common_y_limits) +
    labs(
      title = if (compact) {
        paste0(ifelse(outcome == "PRIMARY", "Primary", "Secondary"),
               " outcomes\n", ifelse(adjusted, "Adjusted", "Unadjusted"))
      } else {
        paste0(ifelse(adjusted, "Adjusted ", ""),
               ifelse(outcome == "PRIMARY", "Primary", "Secondary"), " Outcomes")
      },
      x = NULL,
      y = NULL
    ) +
    scale_x_discrete(labels = function(labels) {
      if (compact) sub("Early Phase I", "Early\nPhase I", labels, fixed = TRUE) else labels
    }) +
    scale_fill_manual(
      values = c(logJAB = "#1f77b4", logP = "#ff7f0e"),
      labels = c(expression("ln" ~ italic("eJAB")["01"]), expression("ln" ~ italic(p)))
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = axis_size),
      axis.text.x = element_text(size = if (compact) 10 else 8),
      plot.title = element_text(size = title_size),
      plot.margin = margin(6, 14, 6, 14),
      legend.key.spacing.x = unit(20, "pt"),
      legend.title = element_text(margin = margin(r = 25))
    )
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = FALSE)
  script_arg <- grep("^--file=", args, value = TRUE)
  script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else "figure4_revised.R"
  script_dir <- dirname(normalizePath(script_file, mustWork = FALSE))
  repo_dir <- normalizePath(file.path(script_dir, "../.."), mustWork = TRUE)
  input_file <- file.path(script_dir, "outputs", "eligible_early_phase_analyses.csv")
  output_file <- file.path(script_dir, "figures", "Figure4.png")

  if (!file.exists(input_file)) {
    stop("Run early_phase_analysis.R before figure4_revised.R")
  }
  analysis_data <- read_csv(input_file, show_col_types = FALSE)

  panel_a <- make_panel(analysis_data, "PRIMARY", adjusted = FALSE)
  panel_b <- make_panel(analysis_data, "SECONDARY", adjusted = FALSE)
  panel_c <- make_panel(analysis_data, "PRIMARY", adjusted = TRUE)
  panel_d <- make_panel(analysis_data, "SECONDARY", adjusted = TRUE)

  figure4 <- ggarrange(
    panel_a, panel_c, panel_b, panel_d,
    nrow = 2,
    ncol = 2,
    labels = c("A", "C", "B", "D"),
    align = "hv",
    common.legend = TRUE,
    legend = "bottom",
    font.label = list(size = 12)
  )

  ggsave(output_file, figure4, width = 10.83, height = 7.17, units = "in", dpi = 300, bg = "white")
  cat("Saved", output_file, "\n")
}
