{
  ##### Chi-Squared Test for Independence (3 × 3 Contingency Table)
  load("~/Downloads/osfstorage-archive/Simulations/Data/chisqSim_jointMulti_R3C3.RData")
  index <- 1; sub6J <- data.frame()
  for (n in nObs) {
    for (omega in OMEGA) {
      temp <- as.data.frame(reportBF[[index]])
      temp$N <- n # grand total
      temp$ES <- omega
      sub6J <- rbind(sub6J, temp)
      index <- index + 1
    }
  }
  sub6J$p <- pchisq(-2 * log(sub6J$JAB / (sub6J$N^2)) * sub6J$N / (sub6J$N-1), 4, lower.tail=F)
  sub6J$test <- "chisq multinom" # fixed total
  
  
  load("~/Downloads/osfstorage-archive/Simulations/Data/chisqSim_indepMulti_R3C3.RData")
  rm(list=setdiff(ls(), c("sub6J", "reportBF", "nObs", "OMEGA")))
  index <- 1; sub6I <- data.frame()
  for (n in nObs) {
    for (omega in OMEGA) {
      temp <- as.data.frame(reportBF[[index]])
      temp$N <- n # grand total
      temp$ES <- omega
      sub6I <- rbind(sub6I, temp)
      index <- index + 1
    }
  }
  sub6I$p <- pchisq(-2 * log(sub6I$JAB / (sub6I$N^2)) * sub6I$N / (sub6I$N-1), 4, lower.tail=F)
  sub6I$test <- "chisq prod-multinom" # fixed row sums
  
}

df <- rbind(sub6J, sub6I)
df$test <- factor(df$test)
set.seed(277); df <- df[sample(1:nrow(df), nrow(df)),] # shuffle to prevent overlapping


library(ggplot2) # v 3.5.2

df$H0 <- factor(ifelse(df$ES==0, "Yes", "No"))
df$ES2 <- factor(df$ES, 
                 levels=c(0.0, 0.1, 0.3, 0.5),
                 labels=c("null", "small", "medium", "large"))

fig <- ggplot(subset(df, N %in% c(45, 300, 1500)), # N = 45, 90, 300, 750, 1500; total number of objects in a 3 × 3 contingency table
              aes(-log(eJAB), -log(ctBF), color=ES2)) +
  geom_point(size=2, alpha=.1, na.rm=T) +
  facet_grid(rows=vars(test), cols=vars(N), scales="free",
             labeller=label_bquote(cols=italic(n)==.(N))) +
  scale_color_manual(values=c("black", "#5494cb", "#ff7411", "#8263d9")) +
  labs(x=expression("ln"~italic("eJAB")[10]),
       y="Dirichlet Bayes Factor",
       color="True effect size is ") +
  scale_x_continuous(limits=c(-10, 15), expand=c(0, 0),
                     breaks=c(-5,-log(3), log(3), 5, 10),
                     labels=c(-5, "-ln3  ", "  ln3", 5, 10)) +
  scale_y_continuous(limits=c(-10, 15), expand=c(0, 0),
                     breaks=c(-5,-log(3), log(3), 5, 10),
                     labels=c(-5, "-ln3", "ln3", 5, 10)) +
  geom_abline(slope=1, color="gray") +
  geom_vline(xintercept=c(-log(3),log(3)), linetype="dashed", color="gray") +
  geom_hline(yintercept=c(-log(3),log(3)), linetype="dashed", color="gray") +
  theme_classic() +
  theme(text=element_text(size=12),
        axis.title.x=element_text(margin=margin(t=12)),
        legend.position="bottom",
        legend.title=element_text(margin=margin(0,0.4,0,0,"cm")),
        legend.text=element_text(margin=margin(0,0.4,0,0.3,"cm")),
        legend.direction="horizontal") +
  guides(colour=guide_legend(override.aes=list(size=4, alpha=1)))

ggsave("~/Downloads/osfstorage-archive/Figure 2/Figure S1.png", fig, width=2100, height=1700, units="px", device="png", dpi=300)

