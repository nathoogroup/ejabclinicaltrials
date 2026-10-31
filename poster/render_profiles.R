# Rebuild plot labels from existing CSV results; no simulation or risk fitting.
# Run: Rscript poster/render_profiles.R [repository-root] [--revision-only]
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(patchwork)
  library(knitr)
})

args <- commandArgs(trailingOnly = TRUE)
revision_only <- "--revision-only" %in% args
args <- args[args != "--revision-only"]
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
repo <- if (length(args)) normalizePath(args[1]) else
  normalizePath(file.path(dirname(script_path), ".."))
source_lines <- readLines(file.path(repo, "bayes_risk_analysis.Rmd"))
output_dir <- file.path(repo, "bayes_risk_outputs/revised_three_comparisons")

chunk <- function(name) {
  opening <- grep(paste0("^```\\{r ", name, "([,}]|$)"), source_lines)
  stopifnot(length(opening) == 1)
  closing <- which(seq_along(source_lines) > opening & source_lines == "```")[1]
  parse(text = source_lines[seq.int(opening + 1L, closing - 1L)])
}
run_assignments <- function(code, names) {
  for (expr in code) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        as.character(expr[[2]]) %in% names) eval(expr, envir = .GlobalEnv)
  }
}
run_plot <- function(name) {
  # Evaluate assignments and ggsave calls, skipping print-only expressions
  # and CSV exports. The plotting chunks do not compute decision risks.
  for (expr in chunk(name)) {
    if (is.call(expr) && as.character(expr[[1]]) %in% c("<-", "ggsave")) {
      invisible(eval(expr, envir = .GlobalEnv))
    }
  }
  message("Rendered ", name)
}

dense_risk <- read.csv(file.path(output_dir, "dense_bayes_risk_grid.csv"))
bic_risk <- read.csv(file.path(output_dir, "large_sample_empirical_bic_risk_grid.csv"))
chisq_dimension_risk <- read.csv(file.path(output_dir, "chisq_dimensions_risk_grid.csv"))
dense_n_grid <- sort(unique(dense_risk$n))
dense_theta_grid <- sort(unique(dense_risk$theta))
bic_n_grid <- sort(unique(bic_risk$n))
bic_theta_grid <- sort(unique(bic_risk$theta))

run_assignments(chunk("simulation-models"), c("test_labels", "bic_test_keys"))
run_assignments(chunk("plotting-functions"), c(
  "facet_labels", "local_alternative_c_values", "local_alternative_contours",
  "plot_criticality_map", "cell_boundaries", "x_bounds", "y_bounds",
  "empirical_x_bounds", "empirical_y_bounds", "plot_data",
  "empirical_bic_plot_data", "chisq_dimensions", "chisq_w_max", "chisq_w_grid",
  "tsbf_dimension_plot_data", "tsbf_estimated_omega_dimension_plot_data"
))
if (revision_only) {
  # Component plots are temporary when updating only the journal revision.
  # The cached simulation results and existing poster assets are not rewritten.
  output_dir <- tempfile("revision-ewde-")
  dir.create(output_dir)
}
plot_chunks <- c("plot-p", "plot-p05", "plot-bic", "plot-tsbf",
                 "plot-tsbf-estimated-omega", "manuscript-figure")
if (!revision_only) {
  plot_chunks <- c(plot_chunks, "combined-figure", "chi-square-dimension-risk-plot")
}
for (name in plot_chunks) run_plot(name)

# The revision uses p <= 0.005 in Panel C and retains p <= 0.05 as Figure S3.
# The existing poster and companion methods draft retain their p <= 0.05 panel.
poster_a <- manuscript_a +
  facet_wrap(~ facet_label, ncol = 3, labeller = label_parsed)
poster_b <- manuscript_b +
  facet_wrap(~ facet_label, ncol = 3, labeller = label_parsed)
poster_d <- bic_vs_ejab + labs(title = expression(
  "D. " * BF["01"]^plain(BIC) <= 1 ~ "versus" ~ eJAB["01"] <= 1
))
four_panel_profiles <- function(panel_c) {
  (poster_a / poster_b / panel_c / poster_d) +
    plot_layout(heights = c(1, 1, 2, 2), guides = "collect") +
    plot_annotation(
      title = "Equally Weighted Decision Error Profiles",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 18))
    )
}
revision_c <- p_vs_ejab + labs(title = expression(
  "C. " * p <= 0.005 ~ "versus" ~ eJAB["01"] <= 1
))
revision_figure <- file.path(
  repo, "bayes_risk_outputs/revised_three_comparisons/06_revision_ewde.png"
)
ggsave(revision_figure, four_panel_profiles(revision_c),
       width = 12, height = 13.2, dpi = 300)
message("Rendered revision Figure 2 with p <= 0.005 in Panel C.")
supplement_profiles <- p05_vs_ejab +
  facet_wrap(~ facet_label, ncol = 2, labeller = label_parsed) +
  labs(title = expression(atop(
    "Equally Weighted Decision Error Profiles",
    p <= 0.05 ~ "versus" ~ eJAB["01"] <= 1
  )))
supplement_figure <- file.path(
  repo, "bayes_risk_outputs/revised_three_comparisons/07_supplement_ewde_p05.png"
)
ggsave(supplement_figure, supplement_profiles,
       width = 8, height = 9.4, dpi = 300)
message("Rendered standalone p <= 0.05 comparison for Figure S3.")

if (!revision_only) {
  ggsave(file.path(repo, "poster/figures/decision_error_profiles.png"),
         four_panel_profiles(manuscript_c), width = 12, height = 13.2, dpi = 300)
  copies <- c(
    "05_manuscript_ewde.png" = "ewde_main.png",
    "01_ejab_vs_p005_all_10_tests.png" = "ewde_p005.png",
    "02_ejab_vs_bic_bf01_six_parametric_tests.png" = "ewde_bic.png"
  )
  for (name in names(copies)) {
    stopifnot(file.copy(file.path(output_dir, name),
                       file.path(repo, "poster/figures", copies[[name]]),
                       overwrite = TRUE))
  }
}
# Keep the active manuscript and reviewer response synchronized. Archived
# submissions are left unchanged.
paper <- file.path(repo, paste0(
  "Reassessing_the_Statistical_Evidence_in_Cancer_Research_",
  "Using_Approximate_Objective_Bayes_Factors"
))
if (dir.exists(paper)) {
  stopifnot(file.copy(
    revision_figure,
    file.path(paper, "REVISION/Figure2.png"), overwrite = TRUE
  ))
  stopifnot(file.copy(
    supplement_figure,
    file.path(paper, "REVISION/FigureS3.png"), overwrite = TRUE
  ))
  response_copies <- c(
    "01_ejab_vs_p005_all_10_tests.png" = "BayesRiskP005.png",
    "02_ejab_vs_p05_all_10_tests.png" = "BayesRiskP05.png",
    "02_ejab_vs_bic_bf01_six_parametric_tests.png" = "BayesRiskBIC.png",
    "03_ejab_vs_tsbf01_three_table_sizes.png" = "BayesRiskTSBF.png",
    "03b_ejab_vs_tsbf01_estimated_omega_three_table_sizes.png" = "BayesRiskTSBFSensitivity.png",
    "04_chisq_dimensions_ejab_vs_p05_and_tsbf.png" = "BayesRiskChiSquare.png"
  )
  if (!revision_only) {
    for (name in names(response_copies)) {
      stopifnot(file.copy(
        file.path(output_dir, name),
        file.path(paper, "response_to_reviewers/Figures", response_copies[[name]]),
        overwrite = TRUE
      ))
    }
  }
}
if (revision_only) unlink(output_dir, recursive = TRUE)
message("Updated titles and legends from cached results; risk CSV files were read only.")
