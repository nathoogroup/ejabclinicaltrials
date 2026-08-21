############################
##### Data Preparation #####
############################

load("~/Downloads/osfstorage-archive/Figure 3/setup.RData")

myTitle <- 22
myAxis <- myLabel <- myTitle - 2 # font size (Figure 3A-C & E-F, excluding Figure 3D)

library(dplyr)    # v1.1.4
library(ggplot2)  # v3.5.2
library(cowplot)  # v1.2.0


## Figure 3A -------------------------------------------------------------------

RERUN <- function(filtered_df, byP = TRUE) {
  
  if (byP) {
    fine <- seq(0.001, 1, length.out = 1000) # define fine alpha grid
    proportions <- sapply(fine, function(a) {
      mean(filtered_df$JAB > 1 & filtered_df$pValue <= a)
    }) # compute proportions over fine grid
    
  } else {
    bf_breaks <- c(1/300, 1/100, 1/30, 1/10, 1/3, 1, 3, 10, 30, 100, 300) # extended canonical BF breaks
    fine <- exp(seq(log(min(bf_breaks)), log(max(bf_breaks)), length.out = 1000)) # define fine grid (log scale)
    proportions <- sapply(fine, function(b) {
      mean(filtered_df$pValue <= 0.05 & filtered_df$JAB >= b)
    })
  }
  
  n <- nrow(filtered_df)
  low <- proportions - 1.96 * sqrt((proportions * (1 - proportions)) / n)
  high <- proportions + 1.96 * sqrt((proportions * (1 - proportions)) / n) # CI
  
  proportion_data <- data.frame(
    x = fine,
    proportion = proportions,
    n = n,
    low = low,
    high = high
  ) # create the data frame
  
  if (byP) {
    alphas_main <- seq(0.01, 1, length.out = 100) # downsample for main plot
    plot_main_data <- proportion_data %>% filter(x %in% alphas_main)
    
    main_plot <- ggplot(plot_main_data, aes(x = x, y = proportion)) +
      geom_line(color = "#1f77b4", linewidth = 2) +
      ylim(0, .8) +
      geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black", linewidth = 1.5) +
      geom_abline(intercept = mean(filtered_df$JAB > 1), slope = 0, linetype = "dotted", color = "red", linewidth = 1.5) +
      labs(
        title = expression("Prop.  " * italic(p) * " \u2264 " * italic(α) * "  &  " * italic("eJAB")["01"] > 1),
        x = expression(italic(α)),
        y = "Observed  Prop."
      )
    
    inset_plot <- ggplot(proportion_data %>% filter(x <= 0.4), aes(x = x, y = proportion)) +
      geom_line(color = "#1f77b4", linewidth = 1) +
      ylim(0, .4) +
      geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black", linewidth = 1)
      
  } else {
    main_plot <- ggplot(proportion_data, aes(x = x, y = proportion)) +
      geom_line(color = "black", linewidth = 1) +
      labs(
        title = expression("Prop.  " * italic(p) * " \u2264 .05  &  " * italic("eJAB")["01"] * " \u2265 " * italic(K)),
        x = expression(italic(K)),
        y = "Observed  Prop."
      ) +
      scale_x_log10(
        breaks = bf_breaks,
        expand = c(0.005, 0.067),
        labels = c(0, 0.01, 0.03, 0.1, 0.33, 1, 3, 10, 30, 100, 300) # scales::label_number(accuracy = 0.01, trim = TRUE)
      )
    
    inset_plot <- ggplot(proportion_data %>% filter(x >= 3), aes(x = x, y = proportion)) +
      geom_line(color = "black", linewidth = 0.8) +
      scale_x_log10(
        breaks = bf_breaks[bf_breaks >= 3],
        labels = c(3, 10, 30, 100, 300) # scales::label_number(accuracy = 0.01, trim = TRUE)
      )
  }
  
  main_plot <- main_plot +
    geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2) +
    theme_minimal() + 
    theme(
      axis.title = element_text(size = myLabel),
      axis.text = element_text(size = myAxis),
      plot.title = element_text(size = myTitle),
      axis.title.y = element_text(margin = margin(r = 15)),
      axis.title.x = element_text(margin = margin(t = 20))
    ) # create main plot
  
  inset_plot <- inset_plot +
    geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title.y = element_blank(),
      axis.text = element_text(size = 12),
      plot.background = element_rect(fill = "white", color = "black")
    ) # create inset plot
  
  ggdraw() +
    draw_plot(main_plot) +
    draw_plot(inset_plot, 
              x = ifelse(byP, 0.6, 0.47), 
              y = ifelse(byP, 0.15, 0.5), 
              width = ifelse(byP, 0.35, 0.5), 
              height = ifelse(byP, 0.35, 0.4)) # combine with inset
}


df1 <- study3_summary %>%
  mutate(
    Significant05 = if_else(pValue <= 0.05, "Significant", "Non-significant"),
    Significant005 = if_else(pValue <= 0.005, "Significant", "Non-significant"),
    BF = case_when(
      JAB < 0.3 ~ " (0, 1/3)",
      JAB >= 0.3 & JAB <= 3 ~ " [1/3, 3]",
      JAB > 3 ~ " (3, \u221E)"
    ),
    BF = factor(BF, levels = c(" (0, 1/3)", " [1/3, 3]", " (3, \u221E)")) # intentional empty space
  )

FigA <- RERUN(df1)


## Figure 3B -------------------------------------------------------------------
## Filter out NAs and pValues exactly equal to 0.05

df2 <- study3_summary %>% filter(has_less_than == FALSE, pValue != 0.05)
FigB <- RERUN(df2)


## Figure 3C -------------------------------------------------------------------

FigC <- RERUN(df1, byP = FALSE)

study3_summary %>% 
  filter(JAB > 3, pValue <= 0.05) %>% 
  nrow() # 487


## Figure 3D -------------------------------------------------------------------

plot_data_summary <- function(data, extra_title = "") {
  
  data_summary2 <- data %>%
    mutate(
      logJAB = log(JAB),
      N_shape = case_when(
        N < 30 ~ "[2, 30)",
        N >= 30 & N < 900 ~ "[30, 900)",
        N >= 900 & N < 30000 ~ "[900, 30K)",
        N >= 30000 & N < 900000 ~ "[30K, 900K)",
        N >= 900000 ~ "[900K, 6030K]"
      ),
      N_shape = factor(
        N_shape,
        levels = c("[2, 30)", "[30, 900)", "[900, 30K)", "[30K, 900K)", "[900K, 6030K]")
      )
    )
  
  total_analyses <- data_summary2 %>%
    pull(analysisId) %>% unique() %>% length()
  
  test_counts <- data_summary2 %>%
    group_by(statisticalMethod) %>%
    summarise(n = n_distinct(analysisId), .groups = "drop")
  
  test_labels <- setNames(
    paste0(test_counts$statisticalMethod, " (", test_counts$n, ")"),
    test_counts$statisticalMethod
  )
  
  data_summary2 <- data_summary2 %>%
    mutate(statisticalMethod = factor(statisticalMethod, levels = names(test_labels)))
  
  # ---- base plot: no x/y scales set here ----
  p_base <- ggplot(data_summary2, aes(x = pValue, y = logJAB, shape = N_shape)) +
    geom_rect(
      data = data.frame(
        ybottom = c(-15, log(1/3), log(3)),
        ytop    = c(log(1/3), log(3), 9),
        color   = c("#A6C7E5", "#B4B4B4", "#FFB84D")
      ),
      aes(xmin = 0, xmax = 1, ymin = ybottom, ymax = ytop, fill = color),
      alpha = 0.2, inherit.aes = FALSE
    ) +
    geom_point(aes(color = statisticalMethod), size = 1) +
    scale_fill_identity() +
    scale_color_brewer(name = "Test", palette = "Paired", labels = test_labels) +
    scale_shape_manual(
      name = "Sample Size",
      values = c(
        "[2, 30)" = 17,
        "[30, 900)" = 15,
        "[900, 30K)" = 16,
        "[30K, 900K)" = 1,
        "[900K, 6030K]" = 8
      )
    ) +
    geom_vline(xintercept = 0.05, linetype = "dashed", color = "black") +
    geom_hline(yintercept = log(1/3), linetype = "dashed", color = "black") +
    geom_hline(yintercept = log(3), linetype = "dashed", color = "black") +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 14),
      plot.title = element_text(size = 14),
      axis.title.x = element_text(margin = margin(t = 15)),
      legend.title = element_text(size = 14),
      legend.key.size = unit(0.8, "cm"),
      legend.text = element_text(size = 14)
    ) +
    guides(
      color = guide_legend(override.aes = list(size = 3, alpha = 1)),
      shape = guide_legend(override.aes = list(size = 3))
    )
  
  # ---- main plot: add scales/limits here (once) ----
  p <- p_base +
    coord_cartesian(xlim = c(-0.02, 1), ylim = c(-15, 9)) +
    scale_x_continuous(name = expression(italic(p) * "-Value"), expand = c(0, 0)) +
    scale_y_continuous(
      name = expression("ln"~italic("eJAB")["01"]),
      breaks = seq(-15, 9, 2),
      expand = c(0, 0)
    ) +
    ggtitle(paste("Findings:", total_analyses, extra_title))
  
  # ---- inset: use base and set different scales (once) ----
  inset_plot <- p_base +
    coord_cartesian(xlim = c(0, 0.05), ylim = c(-3, 3)) +
    scale_x_continuous(name = "", breaks = c(0, 0.025, 0.05), expand = c(0, 0.005)) +
    scale_y_continuous(name = "", breaks = c(-3, 0, 3), expand = c(0, 0)) +
    theme(
      legend.position = "none",
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.background = element_rect(fill = "white", color = "black", linewidth = 1.5),
      axis.text = element_text(size = 8),
      axis.ticks = element_line(linewidth = 0.5)
    )
  
  ggdraw(p) +
    draw_plot(inset_plot, x = 0.25, y = 0.15, width = 0.4, height = 0.4)
}

FigD <- plot_data_summary(study3_summary)


## Figure 3E -------------------------------------------------------------------
## Create a density plot using kernel density estimation

alpha_0.05_results <- posterior_analysis(study3_summary, 5000, 0.05) # roughly 29 seconds of run time

(mean_val <- mean(alpha_0.05_results$data, na.rm = TRUE)) # posterior mean = 35.5%
(sd_val <- sd(alpha_0.05_results$data, na.rm = TRUE)) # posterior SD = 0.22%

FigE <- ggplot(data.frame(value = alpha_0.05_results$data), aes(x = value)) +
  geom_density(fill = "#B4B4B4", alpha = 0.5) +  # Gaussian kernel is default
  geom_vline(xintercept = mean_val, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    y = "",
    x = expression("Prop.  " * italic(p) * " \u2264 .05  &  " * italic("eJAB")["01"] * " > 1/3")
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(
    axis.title = element_text(size = myLabel),
    axis.text = element_text(size = myAxis),
    plot.title = element_text(size = myTitle),
    axis.title.x = element_text(margin = margin(t = 20)),
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  ) # create density plot


## Figure 3F -------------------------------------------------------------------
## Density plot for the significant group

nice_Ns <- c(20, 200, 2000, 20000, 200000)
FigF <- df1 %>%
  filter(Significant05 == "Significant", !is.na(BF)) %>%
  ggplot(aes(x = log(N), fill = BF)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c("#A6C7E5",
                               "#B4B4B4",
                               "#FFB84D")
  ) +
  scale_x_continuous(
    breaks = log(nice_Ns),
    labels = c(expression(2 %*% 10^1), expression(2 %*% 10^2), expression(2 %*% 10^3), expression(2 %*% 10^4), expression(2 %*% 10^5)),
    name = expression("Sample Size (" * italic(p) * " \u2264 .05)")
  ) +
  labs(fill = expression(italic("eJAB")["01"])) +
  theme_minimal() + 
  theme(
    axis.title = element_text(size = myLabel),
    axis.text = element_text(size = myAxis),
    plot.title = element_text(size = myTitle),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(margin = margin(t = 15)),
    legend.title = element_text(size = myTitle),
    legend.key.size = unit(0.8, "cm"),
    legend.text = element_text(size = myTitle),
    legend.position=c(0.85, 0.6) # `legend.position.inside`
  )


##############################
##### Data Visualization #####
##############################

ggsave("~/Downloads/osfstorage-archive/Figure 3/Fig 3A.png", FigA, width=3000, height=2000, units="px", device="jpg", dpi=300)
ggsave("~/Downloads/osfstorage-archive/Figure 3/Fig 3B.png", FigB, width=3000, height=2000, units="px", device="jpg", dpi=300)
ggsave("~/Downloads/osfstorage-archive/Figure 3/Fig 3C.png", FigC, width=3000, height=2000, units="px", device="jpg", dpi=300)
ggsave("~/Downloads/osfstorage-archive/Figure 3/Fig 3D.png", FigD, width=3375 * .8, height=2200, units="px", device="jpg", dpi=300)
ggsave("~/Downloads/osfstorage-archive/Figure 3/Fig 3E.png", FigE, width=2700, height=1200, units="px", device="jpg", dpi=300)
ggsave("~/Downloads/osfstorage-archive/Figure 3/Fig 3F.png", FigF, width=3375, height=1200, units="px", device="jpg", dpi=300)

