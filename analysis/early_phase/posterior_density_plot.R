# Shared presentation for phase-specific posterior contradiction proportions.
# The densities retain their probability-density scale on common axes.
make_contradiction_density_plot <- function(draws, summaries, base_size = 11) {
  phase_levels <- c("Early Phase I", "Phase I", "Phase II")
  phase_colors <- setNames(c("#0072B2", "#D55E00", "#009E73"), phase_levels)
  draws$phase <- factor(draws$phase, levels = phase_levels)
  summaries$phase <- factor(summaries$phase, levels = phase_levels)
  legend_labels <- setNames(
    sprintf("%s\n%.3f [%.3f, %.3f]", summaries$phase,
            summaries$mean, summaries$lower, summaries$upper),
    summaries$phase
  )

  ggplot2::ggplot(draws, ggplot2::aes(value, fill = phase, color = phase)) +
    ggplot2::geom_density(alpha = 0.18, linewidth = 0.85) +
    ggplot2::geom_vline(
      data = summaries,
      ggplot2::aes(xintercept = mean, color = phase),
      linetype = "dashed", linewidth = 0.6, show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = phase_colors, breaks = phase_levels, labels = legend_labels, name = NULL
    ) +
    ggplot2::scale_color_manual(
      values = phase_colors, breaks = phase_levels, labels = legend_labels, name = NULL
    ) +
    ggplot2::scale_x_continuous(breaks = seq(0, 0.4, by = 0.1)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.06))) +
    ggplot2::coord_cartesian(xlim = c(0, 0.4)) +
    ggplot2::labs(
      x = "Proportion with eJAB > 1\namong NHST-significant results", y = "Density"
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(nrow = 1),
      color = ggplot2::guide_legend(nrow = 1)
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      legend.position = "top",
      legend.text = ggplot2::element_text(size = 10.5),
      legend.key.width = grid::unit(0.3, "cm"),
      legend.spacing.x = grid::unit(0.1, "cm"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 11),
      axis.text = ggplot2::element_text(size = 9.5)
    )
}
