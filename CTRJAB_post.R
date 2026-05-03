# CTRJAB_post.R --- post-processing only.
# Skips the network pull and operates on whatever is already in
# ../jtc/euctr.sqlite. Identical extract/categorize/JAB logic to CTRJAB.R.

library(ctrdata)
library(nodbi)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)

source("../rfuncs/JAB.R")

dbc <- nodbi::src_sqlite(dbname = "../jtc/euctr.sqlite", collection = "euctr")

needed <- unique(c(
  unlist(f.primaryEndpointResults()),
  unlist(f.sampleSize()),
  unlist(f.hasResults()),
  unlist(f.trialPhase()),
  "a3_full_title_of_the_trial",
  "a2_eudract_number"
))
df <- dbGetFieldsIntoDf(fields = needed, con = dbc)
for (col in setdiff(needed, names(df))) df[[col]] <- NA
cat("rows in DB:", nrow(df), "\n")

pe <- f.primaryEndpointResults(df)
ss <- f.sampleSize(df)
hr <- f.hasResults(df)
ph <- f.trialPhase(df)

unique_ids <- dbFindIdsUniqueTrials(con = dbc, prefermemberstate = "FR",
                                    include3rdcountrytrials = TRUE)

per_trial <- df %>%
  select(`_id`, a2_eudract_number, a3_full_title_of_the_trial,
         e71_human_pharmacology_phase_i, e72_therapeutic_exploratory_phase_ii,
         e73_therapeutic_confirmatory_phase_iii, e74_therapeutic_use_phase_iv) %>%
  filter(`_id` %in% unique_ids) %>%
  left_join(pe, by = "_id") %>%
  left_join(ss, by = "_id") %>%
  left_join(hr, by = "_id") %>%
  left_join(ph, by = "_id")

ctr_clean <- per_trial %>%
  filter(.hasResults,
         !is.na(.primaryEndpointFirstPvalue),
         .primaryEndpointFirstPvalue > 0,
         .primaryEndpointFirstPvalue < 1,
         !is.na(.primaryEndpointFirstPsize),
         .primaryEndpointFirstPsize > 1) %>%
  rename(pValue         = .primaryEndpointFirstPvalue,
         originalMethod = .primaryEndpointFirstPmethod,
         N              = .primaryEndpointFirstPsize)

categorize_method <- function(m) {
  case_when(
    str_detect(m, "ttest") & str_detect(m, "1sample|paired") ~ "One-sample t-test",
    str_detect(m, "ttest") ~ "Two-sample t-test",
    str_detect(m, "logistic") & str_detect(m, "regression") & str_detect(m, "conditional") ~ "Conditional logistic regression",
    str_detect(m, "logistic") & str_detect(m, "regression") ~ "Logistic regression",
    str_detect(m, "cox") & !str_detect(m, "wilcoxon") ~ "Cox",
    str_detect(m, "logrank|mantelcox") & !str_detect(m, "wilcoxon|signed") ~ "Logrank",
    str_detect(m, "mannwhitney") |
      (str_detect(m, "wilcoxon") & str_detect(m, "rank") & !str_detect(m, "signed")) ~ "Mann-Whitney",
    str_detect(m, "wilcoxon") & str_detect(m, "signed") ~ "Wilcoxon",
    str_detect(m, "kruskal|wallis") ~ "Kruskal-Wallis",
    str_detect(m, "chisquared|chisq|cochranmantelhaenszel|mantelhaenszel") ~ "Chi-square test",
    str_detect(m, "mmrm|mixedeffectmodelrepeatedmeasurement|repeatedmeasures|mixedmodel") ~ "Repeated measures analysis",
    str_detect(m, "linear") & str_detect(m, "regression") ~ "Linear regression",
    str_detect(m, "ancova|anova|analysisofvariance") ~ "ANOVA",
    TRUE ~ NA_character_
  )
}

ctr_cleaned <- ctr_clean %>%
  mutate(statisticalMethod = categorize_method(originalMethod)) %>%
  filter(!is.na(statisticalMethod)) %>%
  mutate(pValue = if_else(str_detect(originalMethod, "1sided"), pmin(pValue * 2, 1), pValue))

euctr_summary <- ctr_cleaned %>%
  mutate(
    model = case_when(
      statisticalMethod %in% c("Two-sample t-test", "One-sample t-test")        ~ "t-test",
      statisticalMethod == "Linear regression"                                   ~ "linear_regression",
      statisticalMethod %in% c("Logistic regression", "Conditional logistic regression") ~ "logistic_regression",
      statisticalMethod == "ANOVA"                                               ~ "anova",
      statisticalMethod == "Kruskal-Wallis"                                      ~ "kruskal_wallis",
      statisticalMethod == "Mann-Whitney"                                        ~ "mann_whitney",
      statisticalMethod == "Wilcoxon"                                            ~ "wilcoxon",
      TRUE ~ NA_character_
    ),
    R = NA_integer_, C = NA_integer_, I = NA_integer_
  ) %>%
  filter(!is.na(model), N > 1) %>%
  rowwise() %>%
  mutate(JAB = tryCatch(JAB01(N, pValue, model, R = R, C = C, I = I),
                        error = function(e) NA_real_)) %>%
  ungroup() %>%
  filter(!is.na(JAB))

dir.create("../csv", showWarnings = FALSE, recursive = TRUE)
write.csv(per_trial,    "../csv/CTR_per_trial.csv",    row.names = FALSE)
write.csv(ctr_cleaned,  "../csv/CTR_cleaned.csv",      row.names = FALSE)
write.csv(euctr_summary,"../csv/CTR_summary.csv",      row.names = FALSE)

cat("\n--- counts at each stage ---\n")
cat("unique trials:           ", length(unique_ids), "\n")
cat("per_trial rows:          ", nrow(per_trial), "\n")
cat("ctr_clean (pval+psize):  ", nrow(ctr_clean), "\n")
cat("ctr_cleaned (categorized):", nrow(ctr_cleaned), "\n")
cat("euctr_summary (with JAB):", nrow(euctr_summary), "\n")

cat("\n--- categorized methods ---\n")
print(ctr_cleaned %>% count(statisticalMethod, sort = TRUE))

cat("\n--- model distribution in summary ---\n")
print(euctr_summary %>% count(model, sort = TRUE))
