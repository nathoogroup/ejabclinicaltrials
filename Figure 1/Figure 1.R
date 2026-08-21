####################
##### Function #####
####################

computeDist <- function(n, ES=0.5, sd=1, loc.alt=0, rand=F) {
  #' Input -
  #' n:           number of observations
  #' ES:          population effect size
  #' sd:          population standard deviation
  #' loc.alt:     whether to simulate under a local alternative;
  #'              0 -- θ₀ + c;            (fixed);
  #'              1 -- θ₀ + c / n;
  #'              2 -- θ₀ + c / √n;       (local);
  #'              3 -- θ₀ + c / n^{1/4};
  #'              4 -- θ₀ + c / ln{n};
  #'              5 -- θ₀ + c ln{n} / √n
  #' rand:        whether to draw the effect (c) from a prior or fix it
  #' 
  #' Output - 
  #' a list of the p-value and $D_{n}=\sqrt{n}\cdot p\cdot(-\ln{p})$
  
  if(rand) ES <- rnorm(1, ES, sd) # random effect
  if (loc.alt==1) {
    multiplier <- 1 / n
  } else if (loc.alt==2) {
    multiplier <- 1 / sqrt(n)
  } else if (loc.alt==3) {
    multiplier <- n^-0.25
  } else if (loc.alt==4) {
    multiplier <- 1 / log(n)
  } else if (loc.alt==5) {
    multiplier <- log(n) / sqrt(n)
  } else {
    multiplier <- 1
  }
  mu <- ES * multiplier
  
  pVal <- t.test(rnorm(n, mu, sd))$p.value # p-value of the one-sample t-test
  list("p"=pVal,
       "D"=sqrt(n) * pVal * -log(pVal) )
}


#######################
##### Simulations #####
#######################

nSim <- 5000 # number of simulation runs for each setting
n <- c(50, 250, 1000, 5000, 25000, 50000) # numbers of observations

set.seed(277)
distAlt2 <- sapply(n, function(x) replicate(nSim, computeDist(x, loc.alt=2)$D)) # θ₀ + c / √n

set.seed(277)
distAlt3 <- sapply(n, function(x) replicate(nSim, computeDist(x, loc.alt=3)$D)) # θ₀ + c / n^{1/4}

loc <- data.frame("dist"=c(distAlt2, distAlt3),
                  "n"=rep(n, each=nSim),
                  "Alt"=factor(rep(c("italic(θ)[0] + italic(c) / sqrt(italic(n))", "italic(θ)[0] + italic(c) / italic(n)^0.25"), each=nSim*length(n)),
                               levels=c("italic(θ)[0] + italic(c) / sqrt(italic(n))", "italic(θ)[0] + italic(c) / italic(n)^0.25")))


############################
##### Data Preparation #####
############################

{
  load("~/Downloads/osfstorage-archive/Theorem/regularity_conditions.RData")
  
  ##### Two-Sample t-Test (balanced group sizes with equal variance)
  sub1 <- data.frame("dist"=c(distMat1),
                     "N"=rep(2 * c(10, 25, 50, 100, 250, 500), each=nrow(distMat1)), # total number of observations
                     "test"="ttest")
  
  
  ##### Simple Linear Regression
  sub2 <- data.frame("dist"=c(distMat2),
                     "N"=rep(c(10, 25, 50, 100, 250, 500), each=nrow(distMat2)),
                     "test"="lm")
  
  
  ##### Simple Logistic Regression (binomial family)
  sub3 <- data.frame("dist"=c(distMat3),
                     "N"=rep(c(10, 25, 50, 100, 250, 500), each=nrow(distMat3)),
                     "test"="glm")
  
  
  ##### One-Way Analysis of Variance (ANOVA) with Four Groups
  sub4 <- data.frame("dist"=c(distMat4),
                     "N"=rep(4 * c(10, 25, 50, 100, 250, 500), each=nrow(distMat4)),
                     "test"="anova")
  
  
  ##### One-Way Repeated-Measures ANOVA (rANOVA) with Four Conditions
  sub5 <- data.frame("dist"=c(distMat5),
                     "N"=rep(3 * c(10, 25, 50, 100, 250, 500), each=nrow(distMat5)), # total number of independent observations
                     "test"="ranova")
  
  
  ##### Chi-Squared Test for Independence (3 × 3 Contingency Table)
  # The smallest sample size is 3 × 100 to ensure accuracy of the chi-squared test with no cells having an expected count lower than 5
  sub6 <- data.frame("dist"=c(distMat6),
                     "N"=rep(3 * c(100, 250, 500, 1000, 2000, 3000), each=nrow(distMat6)),
                     "test"="chisq")
  
  
  ##### Cox Proportional-Hazards Regression (random right-censoring)
  # Note that the smallest sample size is N=25 , because, with right censoring, the Cox model cannot estimate parameters when N=10
  sub7 <- data.frame("dist"=c(distMat7),
                     "N"=rep(c(25, 50, 100, 250, 500), each=nrow(distMat7)), # number of observations (> number of events)
                     "test"="coxph")
  
  
  ##### Wilcoxon Signed-Rank Test
  sub8 <- data.frame("dist"=c(distMat8),
                     "N"=rep(c(10, 25, 50, 100, 250, 500), each=nrow(distMat8)),
                     "test"="wilcox")
  
  
  ##### Mann–Whitney U-Test
  sub9 <- data.frame("dist"=c(distMat9),
                     "N"=rep(2 * c(10, 25, 50, 100, 250, 500), each=nrow(distMat9)), # total number of observations
                     "test"="Utest")
  
  
  ##### Kruskal-Wallis H-Test with Five Groups
  sub10 <- data.frame("dist"=c(distMat10),
                      "N"=rep(5 * c(10, 25, 50, 100, 250, 500), each=nrow(distMat10)),
                      "test"="Htest")
  
  ##### Conditional Logistic Regression (1:1 matched case-control study with two predictors)
  sub11 <- data.frame("dist"=c(distMat11),
                      "N"=rep(c(100, 220, 340, 460, 580, 700), each=nrow(distMat11)), # number of pairs (>= number of matched pairs)
                      "test"="clogit")
}

df <- rbind(sub1, sub2, sub3, sub4, sub5, sub6, sub7, sub8, sub9, sub10)
rm(list=setdiff(ls(), c("df", "sub11", "loc")))
df$test <- factor(df$test,
                  levels=c("ttest", "lm", "glm", "anova", "ranova", "chisq", "coxph",
                           "wilcox", "Utest", "Htest"))





##############################
##### Data Visualization #####
############################## boxplot
library(ggplot2) # v 3.5.2

## Panel A ---------------------------------------------------------------------

fA <- ggplot(df, aes(x=factor(N), y=dist)) +
  geom_boxplot(alpha=.05, na.rm=T) +
  facet_wrap(.~test, scales="free", nrow=2) +
  #labs(x="Sample Size", y=expression(italic(D)[italic(n)])) +
  theme_classic() +
  theme(text=element_text(size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(angle=45, hjust=0.9))


## Panel B ---------------------------------------------------------------------

fB <- ggplot(sub11, aes(x=factor(N), y=dist)) +
  geom_boxplot(alpha=.05, na.rm=T) +
  facet_wrap(.~test, scales="free") +
  theme_classic() +
  theme(text=element_text(size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.margin=margin(t=5, r=1, b=5, l=15),
        #strip.background=element_rect(linewidth=0.55),
        axis.text.x=element_text(angle=45, hjust=0.9))


## Panel C ---------------------------------------------------------------------

fC <- ggplot(loc, aes(x=factor(n), y=dist)) +
  geom_boxplot(alpha=.05, na.rm=T) +
  facet_wrap(.~Alt, scales="free", labeller=label_parsed, nrow=1) +
  theme_classic() +
  theme(text=element_text(size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.margin=margin(t=5, r=5, b=5, l=10),
        axis.text.x=element_text(angle=45, hjust=0.9))


library(ggpubr)   # v0.6.1

f1 <- ggarrange(fA,
                ggarrange(ggarrange(ggplot()+theme_classic(), fB, ggplot()+theme_classic(),
                                    nrow=1, ncol=3, labels=c("", "B", ""),
                                    font.label=list(size=13), widths=c(0.5, 1, 0.5)),
                          fC,
                          nrow=2, ncol=1, labels=c("", "C"),
                          font.label=list(size=13)),
                nrow=1, ncol=2, labels="A", font.label=list(size=13), widths=c(5, 2))

f1 <- annotate_figure(f1,
                      left=text_grob(expression(italic(D)[italic(n)]), rot=90, size=13),
                      bottom=text_grob("Sample Size", size=13))

ggsave("~/Downloads/osfstorage-archive/Figure 1/Figure 1.png", f1, width=3500, height=1500, units="px", device="jpg", dpi=300)

