#!/usr/bin/env Rscript

# Reproduce the application-related supplementary figures from the same
# screened primary- and secondary-outcome set used in the manuscript.

suppressPackageStartupMessages({
  library(cowplot)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(httr)
  library(lme4)
  library(readr)
  library(stringr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else "supplementary_application_figures.R"
script_dir <- dirname(normalizePath(script_file, mustWork = FALSE))
repo_dir <- normalizePath(file.path(script_dir, "../.."), mustWork = TRUE)
input_file <- file.path(script_dir, "outputs", "eligible_early_phase_analyses.csv")
mesh_cache <- file.path(script_dir, "outputs", "mesh_tree_numbers.csv")
figure_dir <- file.path(script_dir, "figures")

if (!file.exists(input_file)) {
  stop("Run early_phase_analysis.R before supplementary_application_figures.R")
}

sample_levels <- c("[2, 30)", "[30, 900)", "[900, 30K)",
                   "[30K, 900K)", "[900K, 6030K]")

analysis_data <- read_csv(input_file, show_col_types = FALSE) %>%
  mutate(
    phase = factor(phase, levels = c("Early Phase I", "Phase I", "Phase II")),
    log_p = log(p_two_sided),
    log_ejab = log(ejab_01),
    statisticalMethod = test_class,
    sample_group = case_when(
      sample_size < 30 ~ "[2, 30)",
      sample_size < 900 ~ "[30, 900)",
      sample_size < 30000 ~ "[900, 30K)",
      sample_size < 900000 ~ "[30K, 900K)",
      TRUE ~ "[900K, 6030K]"
    ),
    sample_group = factor(sample_group, levels = sample_levels)
  )

test_counts <- analysis_data %>%
  count(statisticalMethod, name = "count")
test_levels <- test_counts$statisticalMethod
test_labels <- setNames(
  paste0(test_counts$statisticalMethod, " (", test_counts$count, ")"),
  test_counts$statisticalMethod
)
sample_levels_present <- sample_levels[
  sample_levels %in% as.character(analysis_data$sample_group)
]

query_mesh_batch <- function(id_batch) {
  formatted_ids <- paste0("mesh:", id_batch, collapse = " ")
  query_string <- paste0(
    "PREFIX meshv: <http://id.nlm.nih.gov/mesh/vocab#> ",
    "PREFIX mesh: <http://id.nlm.nih.gov/mesh/> ",
    "SELECT ?meshID ?treeNumber WHERE { VALUES ?meshID { ",
    formatted_ids,
    " } ?meshID meshv:treeNumber ?treeNumber . }"
  )
  query_url <- paste0(
    "https://id.nlm.nih.gov/mesh/sparql?query=",
    URLencode(query_string, reserved = TRUE),
    "&format=CSV&inference=false"
  )
  response <- GET(query_url)
  stop_for_status(response)
  read_csv(I(content(response, as = "text", encoding = "UTF-8")), show_col_types = FALSE)
}

condition_data <- analysis_data %>%
  select(analysisId, conditionIds) %>%
  separate_rows(conditionIds, sep = ";") %>%
  mutate(conditionIds = str_trim(conditionIds)) %>%
  filter(!is.na(conditionIds), conditionIds != "")

mesh_ids <- sort(unique(condition_data$conditionIds))
if (file.exists(mesh_cache)) {
  mesh_mapping <- read_csv(mesh_cache, show_col_types = FALSE)
} else {
  batches <- split(mesh_ids, ceiling(seq_along(mesh_ids) / 100))
  mesh_mapping <- bind_rows(lapply(batches, query_mesh_batch)) %>%
    transmute(
      meshID = sub(".*/", "", meshID),
      tree_code = sub("^([A-Z]\\d+).*", "\\1", sub(".*/", "", treeNumber))
    ) %>%
    distinct()
  write_csv(mesh_mapping, mesh_cache)
}

mesh_data <- condition_data %>%
  left_join(
    mesh_mapping,
    by = c("conditionIds" = "meshID"),
    relationship = "many-to-many"
  ) %>%
  filter(!is.na(tree_code)) %>%
  select(analysisId, tree_code) %>%
  distinct() %>%
  left_join(analysis_data, by = "analysisId")

make_mesh_plot <- function(tree_data, tree_code) {
  n_findings <- n_distinct(tree_data$analysisId)
  base <- ggplot(tree_data, aes(p_two_sided, log_ejab)) +
    annotate("rect", xmin = 0, xmax = 1, ymin = -Inf, ymax = log(1 / 3), fill = "#A6C7E5", alpha = 0.22) +
    annotate("rect", xmin = 0, xmax = 1, ymin = log(1 / 3), ymax = log(3), fill = "grey70", alpha = 0.22) +
    annotate("rect", xmin = 0, xmax = 1, ymin = log(3), ymax = Inf, fill = "#FFB84D", alpha = 0.22) +
    geom_point(
      aes(color = statisticalMethod, shape = sample_group),
      size = 0.38,
      alpha = 0.42
    ) +
    geom_vline(xintercept = 0.05, linetype = "dashed", linewidth = 0.35) +
    geom_hline(yintercept = c(log(1 / 3), log(3)), linetype = "dashed", linewidth = 0.35) +
    scale_color_brewer(
      name = "Test",
      palette = "Paired",
      limits = test_levels,
      labels = test_labels
    ) +
    scale_shape_manual(
      name = "Sample size",
      values = c(17, 15, 16, 1, 8),
      limits = sample_levels_present
    ) +
    theme_bw(base_size = 7) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 6)
    )

  main <- base +
    coord_cartesian(xlim = c(0, 1), ylim = c(-12, 7), expand = FALSE) +
    scale_x_continuous(name = expression(italic(p) * "-value"), breaks = c(0, 0.5, 1)) +
    scale_y_continuous(name = expression("ln " * italic(eJAB)["01"]), breaks = c(-10, 0, 5)) +
    ggtitle(sprintf("%s (n = %s)", tree_code, format(n_findings, big.mark = ",")))

  inset <- base +
    coord_cartesian(xlim = c(0, 0.05), ylim = c(-3, 3), expand = FALSE) +
    scale_x_continuous(name = NULL, breaks = c(0, 0.05), labels = c("0", ".05")) +
    scale_y_continuous(name = NULL, breaks = c(-3, 0, 3)) +
    theme(
      plot.title = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(size = 4.5),
      plot.background = element_rect(fill = "white", color = "black", linewidth = 0.35)
    )

  ggdraw(main) + draw_plot(inset, x = 0.29, y = 0.14, width = 0.60, height = 0.48)
}

shared_legend <- get_legend(
  ggplot(
    analysis_data,
    aes(p_two_sided, log_ejab, color = statisticalMethod, shape = sample_group)
  ) +
    geom_point(size = 1.8, alpha = 1) +
    scale_color_brewer(
      name = "Test",
      palette = "Paired",
      limits = test_levels,
      labels = test_labels
    ) +
    scale_shape_manual(
      name = "Sample size",
      values = c(17, 15, 16, 1, 8),
      limits = sample_levels_present
    ) +
    guides(
      color = guide_legend(
        override.aes = list(size = 2.3, alpha = 1),
        nrow = 2,
        byrow = TRUE,
        order = 2
      ),
      shape = guide_legend(
        override.aes = list(size = 2.3),
        nrow = 1,
        order = 1
      )
    ) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.justification = "center",
      legend.text = element_text(size = 7.5),
      legend.title = element_text(size = 8.5),
      legend.key.width = unit(0.38, "cm"),
      legend.spacing.y = unit(0.04, "cm")
    )
)

tree_codes <- sort(unique(mesh_data$tree_code))
tree_plots <- lapply(tree_codes, function(code) {
  make_mesh_plot(filter(mesh_data, tree_code == code), code)
})
page_membership <- cut(seq_along(tree_plots), breaks = 3, labels = FALSE)
plot_pages <- lapply(split(tree_plots, page_membership), function(plots) {
  panel_grid <- ggarrange(
    plotlist = plots,
    ncol = 3,
    nrow = ceiling(length(plots) / 3),
    align = "hv"
  )
  plot_grid(panel_grid, shared_legend, ncol = 1, rel_heights = c(1, 0.13))
})

for (page_index in seq_along(plot_pages)) {
  ggsave(
    file.path(figure_dir, paste0("FigureS", page_index + 3, ".png")),
    plot_pages[[page_index]],
    width = 10,
    height = 13.2,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}

fit_residuals <- function(metric_data) {
  fit <- lmer(value ~ 1 + (1 | nctId), data = metric_data, REML = TRUE)
  metric_data %>%
    mutate(residual = resid(fit), fitted = value - residual)
}

residual_data <- analysis_data %>%
  select(analysisId, nctId, phase, log_ejab, log_p) %>%
  pivot_longer(c(log_ejab, log_p), names_to = "metric", values_to = "value") %>%
  group_by(metric) %>%
  group_modify(~ fit_residuals(.x)) %>%
  ungroup() %>%
  mutate(metric = factor(
    metric,
    levels = c("log_ejab", "log_p"),
    labels = c("ln~italic(eJAB)[0*1]", "ln~italic(p)")
  ))

residual_plot <- ggplot(residual_data, aes(fitted, residual)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(aes(color = phase), alpha = 0.20, size = 0.5) +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 0.65) +
  facet_wrap(~metric, labeller = label_parsed, scales = "free") +
  scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73")) +
  labs(x = "Fitted value", y = "Residual", color = "Trial phase") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

ggsave(file.path(figure_dir, "FigureS7.png"), residual_plot,
       width = 9.5, height = 5.8, units = "in", dpi = 300, bg = "white")

p_plot <- ggplot(analysis_data, aes(p_two_sided, color = phase, fill = phase)) +
  geom_density(alpha = 0.16, linewidth = 0.75) +
  geom_vline(xintercept = 0.05, color = "black", linetype = "dashed") +
  scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73")) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00", "#009E73")) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = expression(italic(p) * "-value"), y = "Density", color = "Trial phase", fill = "Trial phase") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

n_plot <- ggplot(analysis_data, aes(log(sample_size), color = phase, fill = phase)) +
  geom_density(alpha = 0.16, linewidth = 0.75) +
  scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73")) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00", "#009E73")) +
  scale_x_continuous(
    breaks = log(c(10, 100, 1000, 10000, 100000)),
    labels = c("10", "100", "1,000", "10,000", "100,000")
  ) +
  labs(x = "Sample size (log scale)", y = "Density", color = "Trial phase", fill = "Trial phase") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

phase_distribution_plot <- ggarrange(
  p_plot, n_plot,
  nrow = 1,
  labels = c("A", "B"),
  common.legend = TRUE,
  legend = "top"
)

ggsave(file.path(figure_dir, "FigureS8.png"), phase_distribution_plot,
       width = 9.5, height = 5.2, units = "in", dpi = 300, bg = "white")

cat("Saved Figures S4--S8 from", nrow(analysis_data), "analyses in",
    n_distinct(analysis_data$nctId), "trials across", length(tree_codes),
    "MeSH tree categories.\n")
