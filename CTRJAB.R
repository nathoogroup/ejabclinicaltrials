# CTRJAB.R --- EUCTR analogue of CTGJAB.R.
#
# Pulls EU CTR records (with results posted) into a local SQLite-backed ctrdata
# collection, then uses ctrdata's built-in helper f.primaryEndpointResults() to
# extract the first primary endpoint's p-value, statistical method, and
# per-analysis sample size for each unique trial. Categorization mirrors
# CTGJAB.R, adapted for the normalized method strings ctrdata emits
# (lower-cased, punctuation/whitespace stripped, e.g. "regressioncox").

library(ctrdata)
library(nodbi)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)

source("../rfuncs/JAB.R")

# ---- 1. SQLite-backed ctrdata collection -----------------------------------
db_path <- "../jtc/euctr.sqlite"
dbc <- nodbi::src_sqlite(dbname = db_path, collection = "euctr")

# ---- 2. Populate the DB (idempotent: ctrLoadQueryIntoDb skips already-loaded records)
# Scope: every EUCTR trial with results posted, all phases. EUCTR caps a single
# query at 10,000 trials, so we iterate over year windows on the trial's start date.
year_windows <- seq(2004, as.integer(format(Sys.Date(), "%Y")))
for (yr in year_windows) {
  qt <- sprintf(
    "https://www.clinicaltrialsregister.eu/ctr-search/search?query=&resultsstatus=trials-with-results&dateFrom=%d-01-01&dateTo=%d-12-31",
    yr, yr
  )
  message(sprintf("[%s] EUCTR pull, year %d", format(Sys.time(), "%H:%M:%S"), yr))
  ctrLoadQueryIntoDb(
    queryterm    = qt,
    euctrresults = TRUE,
    con          = dbc
  )
}

# ---- 3. Pull only the fields the helpers need ------------------------------
needed <- unique(c(
  unlist(f.primaryEndpointResults()),
  unlist(f.sampleSize()),
  unlist(f.hasResults()),
  unlist(f.trialPhase()),
  "a3_full_title_of_the_trial",
  "a2_eudract_number"
))
df <- dbGetFieldsIntoDf(fields = needed, con = dbc)

# Helpers internally check for fields belonging to all four registers (ctgov,
# ctgov2, ctis, isrctn). Stub absent columns with NA so the helpers do not error.
for (col in setdiff(needed, names(df))) df[[col]] <- NA

# ---- 4. Run the structured-extraction helpers ------------------------------
pe <- f.primaryEndpointResults(df)   # .primaryEndpointFirstPvalue/Pmethod/Psize
ss <- f.sampleSize(df)               # .sampleSize  (trial-level enrolment)
hr <- f.hasResults(df)               # .hasResults
ph <- f.trialPhase(df)               # .trialPhase

# ---- 5. Deduplicate to one row per unique trial ----------------------------
# EUCTR holds one record per (EUDRACT x reporting member state); keep one per trial.
unique_ids <- dbFindIdsUniqueTrials(
  con                     = dbc,
  prefermemberstate       = "FR",
  include3rdcountrytrials = TRUE
)

per_trial <- df %>%
  select(`_id`, a2_eudract_number, a3_full_title_of_the_trial,
         e71_human_pharmacology_phase_i, e72_therapeutic_exploratory_phase_ii,
         e73_therapeutic_confirmatory_phase_iii, e74_therapeutic_use_phase_iv) %>%
  filter(`_id` %in% unique_ids) %>%
  left_join(pe, by = "_id") %>%
  left_join(ss, by = "_id") %>%
  left_join(hr, by = "_id") %>%
  left_join(ph, by = "_id")

# ---- 6. Clean p-value handling and categorize method -----------------------
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

# Categorize ctrdata's normalized method strings (lower-case, alpha-only).
categorize_method <- function(m) {
  case_when(
    str_detect(m, "ttest") &
      str_detect(m, "1sample|paired")                           ~ "One-sample t-test",
    str_detect(m, "ttest")                                      ~ "Two-sample t-test",

    str_detect(m, "logistic") & str_detect(m, "regression") &
      str_detect(m, "conditional")                              ~ "Conditional logistic regression",
    str_detect(m, "logistic") & str_detect(m, "regression")     ~ "Logistic regression",

    str_detect(m, "cox") & !str_detect(m, "wilcoxon")           ~ "Cox",

    str_detect(m, "logrank|mantelcox") &
      !str_detect(m, "wilcoxon|signed")                         ~ "Logrank",

    str_detect(m, "mannwhitney") |
      (str_detect(m, "wilcoxon") & str_detect(m, "rank") &
         !str_detect(m, "signed"))                              ~ "Mann-Whitney",

    str_detect(m, "wilcoxon") & str_detect(m, "signed")         ~ "Wilcoxon",

    str_detect(m, "kruskal|wallis")                             ~ "Kruskal-Wallis",
    str_detect(m, "chisquared|chisq|cochranmantelhaenszel|mantelhaenszel") ~ "Chi-square test",

    str_detect(m, "mmrm|mixedeffectmodelrepeatedmeasurement|repeatedmeasures|mixedmodel") ~ "Repeated measures analysis",

    str_detect(m, "linear") & str_detect(m, "regression")       ~ "Linear regression",
    str_detect(m, "ancova|anova|analysisofvariance")            ~ "ANOVA",

    TRUE ~ NA_character_
  )
}

ctr_cleaned <- ctr_clean %>%
  mutate(statisticalMethod = categorize_method(originalMethod)) %>%
  filter(!is.na(statisticalMethod)) %>%
  # Double 1-sided p-values when the method string indicates so.
  mutate(pValue = if_else(str_detect(originalMethod, "1sided"),
                          pmin(pValue * 2, 1),
                          pValue))

# ---- 7. Map to JAB models and compute --------------------------------------
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
  mutate(
    JAB = tryCatch(JAB01(N, pValue, model, R = R, C = C, I = I),
                   error = function(e) NA_real_)
  ) %>%
  ungroup() %>%
  filter(!is.na(JAB))

# ---- 8. Persist intermediate + final tables --------------------------------
dir.create("../csv", showWarnings = FALSE, recursive = TRUE)
write.csv(per_trial,    "../csv/CTR_per_trial.csv",    row.names = FALSE)
write.csv(ctr_cleaned,  "../csv/CTR_cleaned.csv",      row.names = FALSE)
write.csv(euctr_summary,"../csv/CTR_summary.csv",      row.names = FALSE)

cat("\nrows: per_trial=", nrow(per_trial),
    " ctr_clean=", nrow(ctr_clean),
    " ctr_cleaned=", nrow(ctr_cleaned),
    " euctr_summary=", nrow(euctr_summary), "\n", sep = "")
