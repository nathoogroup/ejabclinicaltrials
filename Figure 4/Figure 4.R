############################
##### Data Preparation #####
############################

load("~/Downloads/osfstorage-archive/Figure 3/setup.RData")

myTitle <- 12
myAxis <- myLabel <- myTitle - 2 # font size

library(dplyr)    # v1.1.4
library(tidyr)    # v1.3.1
library(ggplot2)  # v3.5.2
library(lme4)     # v1.1-37


myPlot <- function(study3_summary, outcome = "PRIMARY", adjusted = FALSE) {
  
  df_long2 <- study3_summary %>%
    filter(outcomeType == outcome) %>%
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
  
  avg_sample_size <- df_long2 %>%
    group_by(outcomeType, Phase) %>%
    summarize(avg_N = mean(N, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      label = paste0("italic(n) == ", round(avg_N, 0))
      ) # calculate average sample size by phase
  
  if (adjusted) {
    fit_add_adjusted <- function(df) {
      df <- df %>% filter(!is.na(Study), !is.na(Value))
      m <- lmer(Value ~ 1 + (1 | Study), data = df, REML = TRUE)
      mu_hat <- fixef(m)[["(Intercept)"]]
      df$Residual <- resid(m)
      df$mu_hat  <- mu_hat
      df$Value <- df$Residual + df$mu_hat
      df
    }
    
    df_long2 <- df_long2 %>%
      group_by(Metric) %>%
      group_modify(~ fit_add_adjusted(.x))
  }
  
  ggplot(df_long2, aes(x = Phase, y = Value, fill = Metric)) +
    geom_boxplot(alpha = 0.7) +
    geom_hline(yintercept = log(1/3), linetype = "dashed") +
    geom_hline(yintercept = log(3), linetype = "dashed") +
    geom_hline(yintercept = log(0.05), linetype = "dashed", color = "red") +
    geom_text(
      data = avg_sample_size,
      aes(x = Phase, y = ifelse(adjusted, 8, 6), label = label),
      parse = TRUE,
      inherit.aes = FALSE,
      size = 4
    ) +
    labs(title = paste0(ifelse(adjusted, "Adjusted ", ""), 
                        ifelse(outcome=="PRIMARY", "Primary", "Secondary"), 
                        " Outcomes"),
         x = "", y = "") +
    scale_fill_manual(
      values = c("logJAB" = "#1f77b4", "logP" = "#ff7f0e"),
      labels = c(expression("ln"~italic("eJAB")["01"]), expression("ln"~italic(p)))
    ) +
    theme_minimal() + 
    theme(
      axis.title = element_text(size = myLabel),
      axis.text = element_text(size = myAxis),
      plot.title = element_text(size = myTitle),
      legend.key.spacing.x = unit(20, "pt"),
      legend.title = element_text(margin = margin(r = 25))
    )
}


##############################
##### Data Visualization #####
##############################

fig1 <- myPlot(study3_summary)
fig2 <- myPlot(study3_summary, "SECONDARY")
fig3 <- myPlot(study3_summary, adjusted = TRUE)
fig4 <- myPlot(study3_summary, "SECONDARY", adjusted = TRUE)


library(ggpubr)   # v0.6.1

ggarrange(fig1, fig2, fig3, fig4,
          nrow=2, ncol=2, labels=c("A", "B", "C", "D"), align="h",
          common.legend=T, legend="bottom", font.label=list(size=myLabel+3))

ggsave("~/Downloads/osfstorage-archive/Figure 4/Figure 4.png", width=3250, height=2152, units="px", device="jpg", dpi=300)

