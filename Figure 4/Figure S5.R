############################
##### Data Preparation #####
############################

load("~/Downloads/osfstorage-archive/Figure 3/setup.RData")

library(dplyr)    # v1.1.4
library(tidyr)    # v1.3.1
library(lme4)     # v1.1-37


study_col <- "nctId"
test_col  <- "analysisId"

df_long2 <- study3_summary %>%
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

# fit RE per outcomeType + Metric, 
# add residuals and mu-hat, then compute Adjusted = mu_hat + residual
fit_add_adjusted <- function(df) {
  df <- df %>% filter(!is.na(Study), !is.na(Value))
  m <- lmer(Value ~ 1 + (1 | Study), data = df, REML = TRUE)
  mu_hat <- fixef(m)[["(Intercept)"]]
  df$Residual <- resid(m)
  df$mu_hat  <- mu_hat
  df$Adjusted <- df$Residual + df$mu_hat
  df
}

df_adj_long2 <- df_long2 %>%
  group_by(Metric) %>%
  group_modify(~ fit_add_adjusted(.x)) %>% 
  mutate(Fitted = Value - Residual)


##############################
##### Data Visualization #####
##############################

library(ggplot2)  # v3.5.2


df_adj_long2$Metric <- factor(df_adj_long2$Metric,
                              levels = c("logJAB", "logP"),
                              labels = c("ln~italic('eJAB')['01']", "ln~italic(p)"))

p <- df_adj_long2 %>%
  filter(!is.na(Fitted), !is.na(Residual)) %>%
  ggplot(aes(x = Fitted, y = Residual)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = Phase, shape = Phase), alpha = 0.5, size = 1.4) +
  scale_color_brewer(palette = "Paired") +
  geom_smooth(method = "loess", se = FALSE, color = "black") +
  facet_wrap(~ Metric, labeller = label_parsed, scales = "free") +
  labs(
    title = "",
    x = "Fitted Values",
    y = "Residuals",
    color = "Phase"
    ) +
  theme_minimal() +
  theme(text=element_text(size = 14),
        strip.text = element_text(size = 14),
        axis.title.x=element_text(margin=margin(t = 12)),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14)) +
  guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1)))

# roughly 33 seconds of run time
ggsave("~/Downloads/osfstorage-archive/Figure 4/Figure S5.png", p, width=3250, height=2000, units="px", device="jpg", dpi=300)

