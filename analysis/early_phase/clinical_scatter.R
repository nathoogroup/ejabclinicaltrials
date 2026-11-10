# Shared clinical scatter objects used by the supplementary composite and the
# standalone manuscript panel. Callers load dplyr, readr, ggplot2 and cowplot.
# No records or significance classifications are modified here.

read_clinical_scatter_data <- function(eligible_file) {
  read_csv(eligible_file, show_col_types = FALSE) %>%
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
}

make_clinical_scatter <- function(early_phase_data) {
  sample_levels <- c("[2, 30)", "[30, 900)", "[900, 30K)",
                     "[30K, 900K)", "[900K, 6030K]")

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
      name = "Two-sided p-value",
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
    draw_plot(scatter_inset, x = 0.17, y = 0.13, width = 0.42, height = 0.49)

  list(scatter_data = scatter_data, test_counts = test_counts,
       scatter_base = scatter_base, scatter_main = scatter_main,
       scatter_inset = scatter_inset, scatter_plot = scatter_plot)
}
