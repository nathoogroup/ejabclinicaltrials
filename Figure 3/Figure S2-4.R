############################
##### Data Preparation #####
############################

load("~/Downloads/osfstorage-archive/Figure 3/setup.RData")

library(dplyr)    # v1.1.4
library(httr)     # v1.4.7
library(readr)    # v2.1.5
library(tidyr)    # v1.3.1
library(stringr)  # v1.5.1


all_ids <- unlist(strsplit(study3_summary$conditionIds, ";\\s*")) # parse and clean unique MeSH IDs
unique_ids <- unique(all_ids[!is.na(all_ids)])

batches <- split(unique_ids, ceiling(seq_along(unique_ids) / 100)) # split into batches of 100 and query
results_list <- lapply(batches, query_mesh_batch)
 
mesh_df <- bind_rows(results_list) %>%
  mutate(
    meshID = sub(".*/", "", meshID),
    treeNumber = sub(".*/", "", treeNumber)
  ) %>%
  distinct() # clean URIs to extract ID and tree number codes

mesh_df <- mesh_df %>%
  mutate(treeNumber = sub("^([A-Z]\\d+).*", "\\1", treeNumber)) %>%
  distinct()

study3_summary_expanded <- study3_summary %>%
  separate_rows(conditionIds, sep = ";") %>% # split the conditionIds into separate rows
  mutate(conditionIds = str_trim(conditionIds)) %>% # trim whitespace from conditionIds
  filter(!is.na(conditionIds) & conditionIds != "") %>%
  left_join(mesh_df, by = c("conditionIds" = "meshID"))

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
    scale_color_manual(
      name   = "Test",
      values = setNames(brewer.pal(9, "Paired"), 
                        c("1-ttest", "2-ttest", "ANOVA", "Htest", "Utest", "clogit", "glm", "lm", "wilcox")),
      labels = test_labels,
      drop   = FALSE
    ) +
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
      color = guide_legend(override.aes = list(size = 3, alpha = 1), order = 2),
      shape = guide_legend(override.aes = list(size = 3), order = 1)
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

plot_data_per_tree <- function(data) {
  unique_trees <- unique(data$treeNumber)
  
  lapply(unique_trees, function(tree) {
    tree_data <- data[data$treeNumber == tree, ]
    plot_data_summary(tree_data, tree)  # call the abstracted function for each tree
  })
}


##############################
##### Data Visualization #####
##############################

library(ggplot2)       # v3.5.2
library(RColorBrewer)  # v1.1-3
library(cowplot)       # v1.2.0

tree_plots <- plot_data_per_tree(study3_summary_expanded) # 41 plots


library(ggpubr)    # v0.6.1

p1 <- ggarrange(plotlist=tree_plots[1:15], nrow=5, ncol=3, align="v")
p2 <- ggarrange(plotlist=tree_plots[16:30], nrow=5, ncol=3, align="v")
p3 <- ggarrange(plotlist=tree_plots[31:41], nrow=5, ncol=3, align="v")

ggsave("~/Downloads/osfstorage-archive/Figure 3/Figure S2.png", p1, width=2700 * 3, height=2200 * 5, units="px", device="jpg", dpi=300)
ggsave("~/Downloads/osfstorage-archive/Figure 3/Figure S3.png", p2, width=2700 * 3, height=2200 * 5, units="px", device="jpg", dpi=300)
ggsave("~/Downloads/osfstorage-archive/Figure 3/Figure S4.png", p3, width=2700 * 3, height=2200 * 5, units="px", device="jpg", dpi=300)

