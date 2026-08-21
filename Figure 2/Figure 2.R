############################
##### Data Preparation #####
############################

{
  ##### Two-Sample t-Test
  load("~/Downloads/osfstorage-archive/Simulations/Data/ttestSim_totTRUE.RData") # sample size is the total number of observations
  rm(list=setdiff(ls(), c("reportBF", "nSubj", "SMD")))
  index <- 1; sub1 <- data.frame()
  for (n1 in nSubj) {
    for (delta in names(SMD)) {
      temp <- as.data.frame(reportBF[[index]])
      temp$N <- 2.2 * n1 # n2 = 1.2 * n1, unbalanced group
      temp$ES <- delta
      sub1 <- rbind(sub1, temp)
      index <- index + 1
    }
  }
  # The original simulations didn't save the p-values, but the p-values can be back-transformed from JAB.
  sub1$p <- pt(sqrt(-2 * log(sub1$JAB / sqrt(sub1$N))), sub1$N-2, lower.tail=F) * 2
  sub1$eJAB <- sqrt(sub1$N) * exp(-0.5 * qchisq(sub1$p, 1, lower.tail=F) * (sub1$N-1) / sub1$N)
  sub1$test <- "ttest"
  
  
  ##### Simple Linear Regression
  load("~/Downloads/osfstorage-archive/Simulations/Data/lmSim.RData")
  rm(list=setdiff(ls(), c("sub1", 
                          "reportBF", "nObs", "BETA1")))
  index <- 1; sub2 <- data.frame()
  for (n in nObs) {
    for (slope in names(BETA1)) {
      temp <- as.data.frame(reportBF[[index]])
      temp$N <- n
      temp$ES <- slope
      sub2 <- rbind(sub2, temp)
      index <- index + 1
    }
  }
  # The original simulations didn't save the p-values, but the p-values can be back-transformed from JAB.
  sub2$p <- pt(sqrt(-2 * log(sub2$JAB / sqrt(sub2$N))), sub2$N-2, lower.tail=F) * 2
  sub2$eJAB <- sqrt(sub2$N) * exp(-0.5 * qchisq(sub2$p, 1, lower.tail=F) * (sub2$N-1) / sub2$N)
  sub2$test <- "lm"
  
  
  ##### Simple Logistic Regression
  load("~/Downloads/osfstorage-archive/Simulations/Data/glmSim.RData")
  rm(list=setdiff(ls(), c("sub1", "sub2", 
                          "reportBF", "nObs", "BETA1")))
  index <- 1; sub3 <- data.frame()
  for (n in nObs) {
    for (slope in names(BETA1)) {
      temp <- as.data.frame(reportBF[[index]])
      temp$N <- n
      temp$ES <- slope
      sub3 <- rbind(sub3, temp)
      index <- index + 1
    }
  }
  # The original simulations didn't save the p-values, but the p-values can be back-transformed from JAB.
  sub3$p <- pchisq(-2 * log(sub3$JAB / sqrt(sub3$N)), 1, lower.tail=F)
  sub3$eJAB <- sqrt(sub3$N) * exp(-0.5 * qchisq(sub3$p, 1, lower.tail=F) * (sub3$N-1) / sub3$N)
  sub3$test <- "glm"
  
  
  ##### One-Way Analysis of Variance (ANOVA)
  load("~/Downloads/osfstorage-archive/Simulations/Data/anovaSim.RData")
  rm(list=setdiff(ls(), c("sub1", "sub2", "sub3", 
                          "reportBF", "nSubj", "SMD")))
  index <- 1; sub4 <- data.frame()
  for (n in nSubj) {
    for (delta in names(SMD)) {
      temp <- as.data.frame(reportBF[[index]])
      temp$N <- 3 * n # three groups
      temp$ES <- delta
      sub4 <- rbind(sub4, temp)
      index <- index + 1
    }
  }
  # The original simulations didn't save the p-values, but the p-values can be back-transformed from JAB.
  sub4$p <- pf(-log(sub4$JAB_OG / sub4$N), 2, sub4$N-3, lower.tail=F)
  sub4$eJAB <- sqrt(sub4$N) * exp(-0.5 * qchisq(sub4$p, 2, lower.tail=F) * (sqrt(sub4$N)-1) / sqrt(sub4$N))
  sub4$test <- "anova"
  
  
  ##### One-Way Repeated-Measures ANOVA (rANOVA)
  load("~/Downloads/osfstorage-archive/Simulations/Data/rmanovaSim_rho0.9.RData")
  rm(list=setdiff(ls(), c("sub1", "sub2", "sub3", "sub4", 
                          "reportBF", "nSubj", "SMD")))
  index <- 1; sub5H <- data.frame()
  for (n in nSubj) {
    for (delta in names(SMD)) {
      temp <- as.data.frame(reportBF[[index]])
      temp$N <- 2 * n # three conditions
      temp$ES <- delta
      sub5H <- rbind(sub5H, temp)
      index <- index + 1
    }
  }
  # The original simulations didn't save the p-values, but the p-values can be back-transformed from JAB.
  sub5H$p <- pf(-log(sub5H$JAB / sqrt(sub5H$N)) * sqrt(sub5H$N) / (sqrt(sub5H$N)-1),
                2, sub5H$N-2, lower.tail=F)
  sub5H$eJAB <- sqrt(sub5H$N) * exp(-0.5 * qchisq(sub5H$p, 2, lower.tail=F) * (sqrt(sub5H$N)-1) / sqrt(sub5H$N))
  sub5H$test <- "ranova high" # high correlation
  
  
  load("~/Downloads/osfstorage-archive/Simulations/Data/rmanovaSim_rho0.2.RData")
  rm(list=setdiff(ls(), c("sub1", "sub2", "sub3", "sub4", "sub5H", 
                          "reportBF", "nSubj", "SMD")))
  index <- 1; sub5L <- data.frame()
  for (n in nSubj) {
    for (delta in names(SMD)) {
      temp <- as.data.frame(reportBF[[index]])
      temp$N <- 2 * n # three conditions
      temp$ES <- delta
      sub5L <- rbind(sub5L, temp)
      index <- index + 1
    }
  }
  # The original simulations didn't save the p-values, but the p-values can be back-transformed from JAB.
  sub5L$p <- pf(-log(sub5L$JAB / sqrt(sub5L$N)) * sqrt(sub5L$N) / (sqrt(sub5L$N)-1),
                2, sub5L$N-2, lower.tail=F)
  sub5L$eJAB <- sqrt(sub5L$N) * exp(-0.5 * qchisq(sub5L$p, 2, lower.tail=F) * (sqrt(sub5L$N)-1) / sqrt(sub5L$N))
  sub5L$test <- "ranova low" # low correlation
  
  
  ##### Chi-Squared Test for Independence (3 × 3 Contingency Table)
  load("~/Downloads/osfstorage-archive/Simulations/Data/chisqSim_jointMulti_R3C3.RData")
  rm(list=setdiff(ls(), c("sub1", "sub2", "sub3", "sub4", "sub5H", "sub5L", 
                          "reportBF", "nObs", "OMEGA")))
  index <- 1; sub6J <- data.frame()
  for (n in nObs) {
    for (omega in names(OMEGA)) {
      temp <- as.data.frame(reportBF[[index]])
      temp$N <- n # grand total
      temp$ES <- omega
      sub6J <- rbind(sub6J, temp)
      index <- index + 1
    }
  }
  # The original simulations didn't save the p-values, but the p-values can be back-transformed from JAB.
  sub6J$p <- pchisq(-2 * log(sub6J$JAB / (sub6J$N^2)) * sub6J$N / (sub6J$N-1), 4, lower.tail=F)
  sub6J$test <- "chisq multinom" # fixed total
  
  
  load("~/Downloads/osfstorage-archive/Simulations/Data/chisqSim_indepMulti_R3C3.RData")
  rm(list=setdiff(ls(), c("sub1", "sub2", "sub3", "sub4", "sub5H", "sub5L", "sub6J", 
                          "reportBF", "nObs", "OMEGA")))
  index <- 1; sub6I <- data.frame()
  for (n in nObs) {
    for (omega in names(OMEGA)) {
      temp <- as.data.frame(reportBF[[index]])
      temp$N <- n # grand total
      temp$ES <- omega
      sub6I <- rbind(sub6I, temp)
      index <- index + 1
    }
  }
  # The original simulations didn't save the p-values, but the p-values can be back-transformed from JAB.
  sub6I$p <- pchisq(-2 * log(sub6I$JAB / (sub6I$N^2)) * sub6I$N / (sub6I$N-1), 4, lower.tail=F)
  sub6I$test <- "chisq prod-multinom" # fixed row sums
  
  
  ##### Cox Proportional-Hazards Regression
  load("~/Downloads/osfstorage-archive/Simulations/Data/coxphSim.RData")
  rm(list=setdiff(ls(), c("sub1", "sub2", "sub3", "sub4", "sub5H", "sub5L", "sub6J", "sub6I", 
                          "reportBF", "reportP", "nObs", "BETA1")))
  index <- 1; sub7 <- data.frame()
  for (n in nObs) {
    for (slope in names(BETA1)) {
      temp <- as.data.frame(reportBF[[index]])
      temp$N <- n #  <== number of observations
      temp$ES <- slope
      temp$p <- reportP[[index]][,1]
      temp$nEvent <- reportP[[index]][,2]
      sub7 <- rbind(sub7, temp)
      index <- index + 1
    }
  }
  sub7$eJAB <- sqrt(sub7$nEvent) * exp(-0.5 * qchisq(sub7$p, 1, lower.tail=F) * (sub7$nEvent-1) / sub7$nEvent)
  sub7$test <- "cox"
  
}


# [1] "22-ttest"         "44-ttest"         "66-ttest"         "220-ttest"        "1100-ttest"       "10-lm"           
# [7] "20-lm"            "30-lm"            "100-lm"           "500-lm"           "10-glm"           "20-glm"          
# [13] "30-glm"           "100-glm"          "500-glm"          "30-anova"         "60-anova"         "90-anova"        
# [19] "300-anova"        "1500-anova"       "20-ranova hi"     "40-ranova hi"     "60-ranova hi"     "20-ranova lo"    
# [25] "40-ranova lo"     "60-ranova lo"     "45-chisq joint"   "90-chisq joint"   "300-chisq joint"  "750-chisq joint" 
# [31] "1500-chisq joint" "45-chisq indep"   "90-chisq indep"   "300-chisq indep"  "750-chisq indep"  "1500-chisq indep"
# [37] "10-cox"         "20-cox"         "30-cox"         "100-cox"        "500-cox"       





##############################
##### Data Visualization #####
############################## scatter plot
library(ggplot2) # v 3.5.2

## Panel A, eJAB ---------------------------------------------------------------

select <- c("SD", "eJAB", "N", "ES", "p", "test")
df <- rbind(subset(sub1, N==66)[select], 
            subset(sub2, N==30)[select], 
            subset(sub3, N==30)[select], 
            subset(sub4, N==90)[select], 
            subset(sub5H, N==40)[select], 
            subset(sub5L, N==40)[select], 
            setNames(subset(sub6J, N==300)[c("TSBF", "eJAB", "N", "ES", "p", "test")], select), #  <== no SDdr gold standard
            setNames(subset(sub6I, N==300)[c("TSBF", "eJAB", "N", "ES", "p", "test")], select), #  <== no SDdr gold standard
            subset(sub7, N==30)[select])
df$test <- factor(df$test,
                  levels=c("ttest", "lm", "glm", "anova", "ranova high", "ranova low", "chisq multinom", "chisq prod-multinom", "cox"))
df$label <- paste0(ifelse(df$test == "cox", "italic(n)^'*' == ", "italic(n) == "), df$N)
set.seed(277); df2 <- df[sample(1:nrow(df), nrow(df)),] # shuffle to prevent overlapping

figA <- ggplot(subset(na.omit(df2), SD > 0), #  <== negative SDdr in ANOVA and rANOVA due to estimation issues
              aes(-log(eJAB), -log(SD),
                  color=factor(cut(p, breaks=c(-1E-5, .01, .05, .1, 1),
                                   labels=c("[0, .01]", "(.01, .05]", "(.05, .1]", "(.1, 1]"))),
                  shape=factor(cut(p, breaks=c(-1E-5, .01, .05, .1, 1),
                                   labels=c("[0, .01]", "(.01, .05]", "(.05, .1]", "(.1, 1]"))))) +
  geom_point(size=2, alpha=.1, na.rm=T) +
  facet_wrap(.~test, nrow=2) +
  scale_color_manual(values=c("#8263d9", "#ff7411", "#5494cb", "black")) +
  labs(x=expression("ln"~italic("eJAB")[10]),
       y="MCMC and Test-Statistic Bayes Factors\n",
       color=expression(paste(italic("p-"), "Value")),
       shape=expression(paste(italic("p-"), "Value"))) + #title="Accuracy of the Bayes-Factor Approximation", subtitle="Fixed Sample Size"
  scale_x_continuous(limits=c(-3, 5), expand=c(0, 0),
                     breaks=c(-log(3), 0, log(3)),
                     labels=c("-ln3", 0, "ln3")) +
  scale_y_continuous(limits=c(-3, 5), expand=c(0, 0),
                     breaks=c(-log(3), 0, log(3)),
                     labels=c("-ln3", 0, "ln3")) +
  geom_abline(slope=1, color="gray") +
  geom_vline(xintercept=c(-log(3),log(3)), linetype="dashed", color="gray") +
  geom_hline(yintercept=c(-log(3),log(3)), linetype="dashed", color="gray") +
  geom_text(data=df, aes(x=-1.2, y=1.5, label=label), alpha=.01, parse=T, col="red") +
  theme_classic() +
  theme(text=element_text(size=12),
        axis.title.x=element_text(margin=margin(t=12)),
        legend.position=c(0.9, 0.25)) + #`legend.position.inside`
  guides(colour=guide_legend(override.aes=list(size=2, alpha=1)))


## Panel B, BIC approx ---------------------------------------------------------

select <- c("BIC_approx", "SD", "N", "ES", "p", "test")
df <- rbind(subset(sub1, N==66)[select], 
            subset(sub2, N==30)[select], 
            subset(sub3, N==30)[select], 
            subset(sub4, N==90)[select], 
            setNames(subset(sub5H, N==40)[c("SBC_approx", "SD", "N", "ES", "p", "test")], select), # or "NM16"
            setNames(subset(sub5L, N==40)[c("SBC_approx", "SD", "N", "ES", "p", "test")], select),
            setNames(subset(sub6J, N==300)[c("BIC_approx", "TSBF", "N", "ES", "p", "test")], select), #  <== no SDdr gold standard
            setNames(subset(sub6I, N==300)[c("BIC_approx", "TSBF", "N", "ES", "p", "test")], select), #  <== no SDdr gold standard
            subset(sub7, N==30)[select])
df$test <- factor(df$test,
                  levels=c("ttest", "lm", "glm", "anova", "ranova high", "ranova low", "chisq multinom", "chisq prod-multinom", "cox"))
df$label <- paste0(ifelse(df$test == "cox", "italic(n)^'*' == ", "italic(n) == "), df$N)
set.seed(277); df2 <- df[sample(1:nrow(df), nrow(df)),] # shuffle to prevent overlapping

figB <- ggplot(subset(na.omit(df2), SD > 0), #  <== negative SDdr in ANOVA and rANOVA due to estimation issues
              aes(-log(BIC_approx), -log(SD),
                  color=factor(cut(p, breaks=c(-1E-5, .01, .05, .1, 1),
                                   labels=c("[0, .01]", "(.01, .05]", "(.05, .1]", "(.1, 1]"))),
                  shape=factor(cut(p, breaks=c(-1E-5, .01, .05, .1, 1),
                                   labels=c("[0, .01]", "(.01, .05]", "(.05, .1]", "(.1, 1]"))))) +
  geom_point(size=2, alpha=.1, na.rm=T) +
  facet_wrap(.~test, nrow=2) +
  scale_color_manual(values=c("#8263d9", "#ff7411", "#5494cb", "black")) +
  labs(x=expression("ln"~italic("BF")[10]^"(BIC)"),
       y="MCMC and Test-Statistic Bayes Factors\n",
       color=expression(paste(italic("p-"), "Value")),
       shape=expression(paste(italic("p-"), "Value"))) + #title="Accuracy of the Bayes-Factor Approximation", subtitle="Fixed Sample Size"
  scale_x_continuous(limits=c(-3, 5), expand=c(0, 0),
                     breaks=c(-log(3), 0, log(3)),
                     labels=c("-ln3", 0, "ln3")) +
  scale_y_continuous(limits=c(-3, 5), expand=c(0, 0),
                     breaks=c(-log(3), 0, log(3)),
                     labels=c("-ln3", 0, "ln3")) +
  geom_abline(slope=1, color="gray") +
  geom_vline(xintercept=c(-log(3),log(3)), linetype="dashed", color="gray") +
  geom_hline(yintercept=c(-log(3),log(3)), linetype="dashed", color="gray") +
  geom_text(data=df, aes(x=2.2, y=-1.5, label=label), alpha=.01, parse=T, col="red") +
  theme_classic() +
  theme(text=element_text(size=12),
        axis.title.x=element_text(margin=margin(t=12)),
        legend.position=c(0.9, 0.25)) + #`legend.position.inside`
  guides(colour=guide_legend(override.aes=list(size=2, alpha=1)))


library(ggpubr)   # v0.6.1

ggarrange(figA, figB,
          nrow=2, ncol=1, labels=c("A", "B"), align="v",
          font.label=list(size=13))

ggsave("~/Downloads/osfstorage-archive/Figure 2/Figure 2.png", width=3250, height=3400, units="px", device="jpg", dpi=300)

