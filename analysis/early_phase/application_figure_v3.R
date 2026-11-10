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
revision_dir <- normalizePath(file.path(script_dir, "figures"))
output_dir <- file.path(script_dir, "outputs")

phase_levels <- c("Early Phase I", "Phase I", "Phase II")

eligible_file <- file.path(output_dir, "eligible_early_phase_analyses.csv")
if (!file.exists(eligible_file)) {
  stop("Run early_phase_analysis.R before application_figure_v3.R")
}

source(file.path(script_dir, "clinical_scatter.R"))
early_phase_data <- read_clinical_scatter_data(eligible_file)
scatter_objects <- make_clinical_scatter(early_phase_data)
scatter_data <- scatter_objects$scatter_data
scatter_main <- scatter_objects$scatter_main
scatter_inset <- scatter_objects$scatter_inset
scatter_plot <- scatter_objects$scatter_plot

# This renderer writes the six-panel figure and a large standalone scatterplot.
# Small-a posterior densities
# are generated separately; historical a=1 and bootstrap artifacts are retained
# for reproducibility but are not recomputed or plotted here.
source(file.path(script_dir, "jlp_phase_panels.R"))
threshold_grid <- sort(unique(c(
  exp(seq(log(1 / 300), log(300), length.out = 1000)),
  seq(3, 5, length.out = 1001), 1 / 3, 1, 3,
  scatter_data$JAB[scatter_data$statistically_significant &
                   scatter_data$JAB >= 1 / 300 & scatter_data$JAB <= 300]
)))
threshold_sweep <- bind_rows(lapply(phase_levels, function(phase_name) {
  phase_frame <- filter(scatter_data, exact_phase == phase_name)
  counts <- vapply(threshold_grid, function(threshold) {
    sum(phase_frame$statistically_significant & phase_frame$JAB > threshold)
  }, integer(1))
  stopifnot(all(diff(counts) <= 0))
  tibble(
    phase = phase_name,
    bf_threshold = threshold_grid,
    screened_results = counts,
    eligible_analyses = nrow(phase_frame),
    proportion = counts / nrow(phase_frame)
  )
}))
threshold_plot <- make_jlp_threshold_plot(threshold_sweep)

# Outcome-distribution panels C-F display folded two-sided p-values; significance
# membership and eJAB remain those of the original directional analysis.
evidence_environment <- new.env(parent = environment())
sys.source(file.path(script_dir, "figure4_revised.R"), envir = evidence_environment)
evidence_data <- read_csv(eligible_file, show_col_types = FALSE)
outcome_panels <- list(
  evidence_environment$make_panel(evidence_data, "PRIMARY", FALSE, compact = TRUE),
  evidence_environment$make_panel(evidence_data, "PRIMARY", TRUE, compact = TRUE),
  evidence_environment$make_panel(evidence_data, "SECONDARY", FALSE, compact = TRUE),
  evidence_environment$make_panel(evidence_data, "SECONDARY", TRUE, compact = TRUE)
)
outcome_legend <- get_legend(
  outcome_panels[[1]] + theme(legend.position = "bottom", legend.text = element_text(size = 12))
)
outcome_panels <- lapply(outcome_panels, function(panel) {
  panel + theme(legend.position = "none", plot.margin = margin(6, 5, 6, 5))
})
summary_row <- plot_grid(
  threshold_plot, scatter_plot,
  labels = c("A", "B"), label_size = 16,
  nrow = 1, rel_widths = c(1, 1.3)
)
outcome_row <- plot_grid(
  plotlist = outcome_panels,
  labels = c("C", "D", "E", "F"), label_size = 16,
  nrow = 1, align = "h", axis = "tb"
)
outcome_section <- plot_grid(
  outcome_row, outcome_legend, ncol = 1, rel_heights = c(1, 0.09)
)
application_figure <- plot_grid(
  summary_row, outcome_section,
  ncol = 1, rel_heights = c(1.55, 1.3)
)

ggsave(
  file.path(revision_dir, "Figure3_V3.png"),
  application_figure,
  width = 12,
  height = 8.5,
  dpi = 300,
  bg = "white"
)

# A full-page version retains the same records, axes and evidence bands as
# panel B, with larger text and a spacious inset for standalone inspection.
scatter_main_large <- scatter_main +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 11.5),
    legend.title = element_text(size = 12.5),
    legend.key.height = unit(0.48, "cm"),
    plot.margin = margin(8, 12, 8, 8)
  )
scatter_inset_large <- scatter_inset +
  theme(
    axis.text = element_text(size = 11),
    axis.ticks = element_line(linewidth = 0.45),
    plot.margin = margin(6, 6, 6, 6)
  )
scatter_large <- ggdraw(scatter_main_large) +
  draw_plot(scatter_inset_large, x = 0.20, y = 0.12, width = 0.47, height = 0.56)
ggsave(
  file.path(revision_dir, "Clinical_scatter_large.png"),
  scatter_large,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)
