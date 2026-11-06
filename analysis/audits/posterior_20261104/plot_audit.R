#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(readr); library(dplyr); library(tidyr); library(ggplot2)})
audit_dir <- "analysis/audits/posterior_20261104"
d <- read_csv(file.path(audit_dir, "results/posterior_comparison.csv"), show_col_types = FALSE) %>%
  filter(p_rule == "<= .05", prior_per_cell == 1) %>%
  mutate(label = case_when(dataset == "Original" ~ "Original: all trials",
                          group == "All" ~ "Revision: all early phases",
                          TRUE ~ paste("Revision:", group)),
         label = factor(label, levels = rev(c("Original: all trials", "Revision: all early phases",
                             "Revision: Early Phase I", "Revision: Phase I", "Revision: Phase II"))),
         threshold = factor(if_else(threshold == 1, "eJAB > 1", "eJAB > 1/3"),
                            levels = c("eJAB > 1/3", "eJAB > 1")))
p <- ggplot(d, aes(y = label)) +
  geom_segment(aes(x = empirical, xend = posterior_mean, yend = label), color = "grey72", linewidth = .6) +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = .14, linewidth = .7, color = "#0072B2") +
  geom_point(aes(x = posterior_mean, color = "Posterior mean and 95% credible interval", shape = "Posterior mean and 95% credible interval"), size = 2.6) +
  geom_point(aes(x = empirical, color = "Empirical proportion", shape = "Empirical proportion"), size = 2.8) +
  facet_grid(. ~ threshold) +
  scale_x_continuous(limits = c(0, .6), breaks = seq(0, .6, .1), labels = scales::label_percent(accuracy = 1)) +
  scale_color_manual(values = c("Empirical proportion" = "#111111", "Posterior mean and 95% credible interval" = "#0072B2")) +
  scale_shape_manual(values = c("Empirical proportion" = 18, "Posterior mean and 95% credible interval" = 16)) +
  labs(x = "Proportion among NHST-significant results", y = NULL, color = NULL, shape = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        strip.background = element_rect(fill = "#EDF2F5"), strip.text = element_text(face = "bold"),
        legend.position = "bottom", legend.box = "vertical",
        legend.text = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(margin = margin(t = 8)), panel.spacing = grid::unit(.7, "cm")) +
  guides(color = guide_legend(ncol = 1), shape = guide_legend(ncol = 1))
ggsave(file.path(audit_dir, "results/empirical_posterior_comparison.png"), p,
       width = 9, height = 4.2, units = "in", dpi = 300)
