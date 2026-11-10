#!/usr/bin/env Rscript

# One-column version of Supplementary Figure S7B, including its inset.
# Uses every eligible record and the existing folded plotting p-value.
suppressPackageStartupMessages({
  library(cowplot)
  library(dplyr)
  library(ggplot2)
  library(readr)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
eligible_file <- file.path(script_dir, "outputs", "eligible_early_phase_analyses.csv")
output_file <- file.path(script_dir, "figures", "Clinical_scatter_column.png")
input_hash <- tools::md5sum(eligible_file)
source(file.path(script_dir, "clinical_scatter.R"))
plot_data <- read_clinical_scatter_data(eligible_file)
objects <- make_clinical_scatter(plot_data)

short_names <- c(
  "Kruskal-Wallis test" = "Kruskal-Wallis",
  "Linear regression" = "Linear regression",
  "Logistic regression" = "Logistic regression",
  "Mann-Whitney test" = "Mann-Whitney",
  "One-sample/paired t-test" = "One-sample/paired t",
  "One-way ANOVA" = "One-way ANOVA",
  "Two-sample t-test" = "Two-sample t",
  "Wilcoxon signed-rank test" = "Wilcoxon signed-rank"
)
test_labels <- setNames(
  paste0(short_names[objects$test_counts$statisticalMethod],
         " (", objects$test_counts$count, ")"),
  objects$test_counts$statisticalMethod
)
stopifnot(!anyNA(test_labels))

main <- objects$scatter_main +
  scale_color_brewer(name = "Test (records)", palette = "Paired", labels = test_labels) +
  guides(
    color = guide_legend(ncol = 2, byrow = TRUE, order = 1,
                         title.position = "top",
                         override.aes = list(size = 1.7, alpha = 1)),
    shape = guide_legend(nrow = 1, order = 2, title.position = "left",
                         override.aes = list(size = 1.7))
  ) +
  theme(
    axis.title = element_text(size = 8.5),
    axis.text = element_text(size = 7.5),
    axis.title.x = element_text(margin = margin(t = 4)),
    axis.title.y = element_text(margin = margin(r = 4)),
    legend.position = "bottom", legend.box = "vertical",
    legend.title = element_text(size = 7, face = "bold"),
    legend.text = element_text(size = 6.5),
    legend.key.width = unit(0.14, "in"),
    legend.key.height = unit(0.11, "in"),
    legend.spacing.y = unit(2, "pt"),
    legend.box.spacing = unit(1, "pt"),
    legend.margin = margin(1, 0, 1, 0),
    plot.margin = margin(5, 10, 2, 1)
  )
# Finer marks avoid overplotting when the same records occupy one column.
main$layers[[2]]$aes_params$size <- 0.32
main$layers[[3]]$aes_params$linewidth <- 0.35
main$layers[[4]]$aes_params$linewidth <- 0.35
legend <- get_legend(main)
main <- main + theme(legend.position = "none")
inset <- objects$scatter_inset +
  theme(
    axis.text = element_text(size = 7),
    axis.ticks = element_line(linewidth = 0.25),
    panel.border = element_rect(linewidth = 0.4),
    plot.background = element_rect(fill = "white", color = "black", linewidth = 0.45),
    plot.margin = margin(3, 4, 3, 3)
  )
panel <- ggdraw(main) +
  draw_plot(inset, x = 0.35, y = 0.17, width = 0.61, height = 0.50)
column_plot <- plot_grid(panel, legend, ncol = 1, rel_heights = c(1, 0.40))
ggsave(output_file, column_plot, width = 3.35, height = 4.1,
       units = "in", dpi = 600, bg = "white")

stopifnot(identical(input_hash, tools::md5sum(eligible_file)))
cat("Rendered", nrow(plot_data), "records from", n_distinct(plot_data$nctId),
    "trials to", output_file, "\n")
