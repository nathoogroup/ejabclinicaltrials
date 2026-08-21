############################
##### Data Preparation #####
############################

load("~/Downloads/osfstorage-archive/Figure 3/setup.RData")

library(dplyr)    # v1.1.4
library(tidyr)    # v1.3.1


study_col <- "nctId"
test_col  <- "analysisId"

df_long <- study3_summary %>%
  pivot_longer(
    cols = c(EARLY_PHASE1, PHASE1, PHASE2, PHASE3, PHASE4),
    names_to = "Phase",
    values_to = "Indicator"
  ) %>%
  filter(Indicator == 1,
         !is.na(JAB), !is.na(pValue), JAB > 0, pValue > 0) %>%
  mutate(
    Phase = case_when(
      Phase == "EARLY_PHASE1" ~ "Early Phase 1",
      Phase == "PHASE1" ~ "Phase 1",
      Phase == "PHASE2" ~ "Phase 2",
      Phase == "PHASE3" ~ "Phase 3",
      Phase == "PHASE4" ~ "Phase 4"
    ),
    Phase = factor(Phase, levels = c("Early Phase 1","Phase 1","Phase 2","Phase 3","Phase 4")),
    logP = log(pValue),
    logJAB = log(JAB),
    Study = !!sym("nctId"),
    Test  = !!sym("analysisId")
  ) %>%
  pivot_longer(
    cols = c(logJAB, logP),
    names_to = "Metric",
    values_to = "Value"
  )


##############################
##### Data Visualization #####
##############################

library(ggplot2)  # v3.5.2
library(ggridges) # v0.5.7


## by p-value ------------------------------------------------------------------

kw_result <- kruskal.test(log(pValue) ~ Phase, data = df_long)
(p_value <- kw_result$p.value)
# p_text <- ifelse(p_value < 0.001, "p < 0.001", paste("p =", round(p_value, 3)))

figP <- ggplot(df_long, aes(x = pValue, y = Phase, fill = Phase)) +
  geom_density_ridges(aes(height = after_stat(density)), stat = "density",
                      scale = 1, rel_min_height = 0, alpha = 0.5, color = "black") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), name = expression(italic(p) * "-Value")) +
  geom_vline(xintercept = 0.05, color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = 0.2, y = 0.7, label = "italic(α) == .05", parse = TRUE, color = "red", size = 4.5) +
  labs(y="") +
  theme_minimal() + 
  theme(legend.position = "none")  + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title.x=element_text(margin=margin(t = 12)))

## by sample size --------------------------------------------------------------

kw_result <- kruskal.test(log(N) ~ Phase, data = df_long)
(p_value <- kw_result$p.value)
# p_text <- ifelse(p_value < 0.001, "p < 0.001", paste("p =", round(p_value, 3)))

figN <- ggplot(df_long, aes(x = log(N), y = Phase, fill = Phase)) +
  geom_density_ridges(aes(height = after_stat(density)), stat = "density",
                      scale = 1, rel_min_height = 0, alpha = 0.5, color = "black") +
  scale_x_continuous(breaks = log(c(20, 200, 2000, 20000)), 
                     labels = c(expression(2 %*% 10^1), 
                                expression(2 %*% 10^2), 
                                expression(2 %*% 10^3), 
                                expression(2 %*% 10^4)), 
                     name = "Sample Size") +
  labs(y="") +
  theme_minimal() + 
  theme(legend.position = "none")  + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x=element_text(margin=margin(t = 12)))


library(ggpubr)   # v0.6.1

ggarrange(figP, figN, nrow=1, ncol=2, align="h")

ggsave("~/Downloads/osfstorage-archive/Figure 4/Figure S6.png", width=2500, height=1300, units="px", device="jpg", dpi=300)

