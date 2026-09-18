#!/usr/bin/env Rscript

# Protocol-aligned reanalysis of early-phase ClinicalTrials.gov results.
#
# The archived input used here is the August 7, 2025 study3 snapshot. Because
# that snapshot predates the July 2026 preregistration, this script implements
# the protocol's eligibility rules as a transparent retrospective reanalysis;
# it is not described as a prospective confirmatory analysis.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else "early_phase_analysis.R"
script_dir <- dirname(normalizePath(script_file, mustWork = FALSE))
response_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
project_dir <- normalizePath(file.path(script_dir, "../.."), mustWork = TRUE)
input_zip <- file.path(project_dir, "jtc", "study3.zip")
output_dir <- file.path(script_dir, "outputs")
figure_dir <- file.path(script_dir, "figures")
response_figure_dir <- figure_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(response_figure_dir, recursive = TRUE, showWarnings = FALSE)

bad_analysis_ids <- c(
  "9f08f829a4d81c7b436b818dd972e1b4",
  "01bf66ca01ef01f21b9e3bfd18d4c81f",
  "25638018d5a7c1c3cd22fcddd0cf44d8",
  "b39984376c029bac02dc6de78af004cd"
)

classify_method <- function(method) {
  m <- str_to_lower(coalesce(method, ""))
  case_when(
    str_detect(m, "\\bt[- ]?test") & str_detect(m, "one[- ]?sample|paired") ~ "One-sample/paired t-test",
    str_detect(m, "\\bt[- ]?test") ~ "Two-sample t-test",
    str_detect(m, "logistic") & str_detect(m, "regression") & str_detect(m, "conditional") ~ "Conditional logistic regression",
    str_detect(m, "logistic") & str_detect(m, "regression") ~ "Logistic regression",
    str_detect(m, "mann[- ]?whitney") |
      (str_detect(m, "wilcoxon") & str_detect(m, "rank") & !str_detect(m, "signed")) ~ "Mann-Whitney test",
    str_detect(m, "wilcoxon") & str_detect(m, "signed") & str_detect(m, "rank") ~ "Wilcoxon signed-rank test",
    str_detect(m, "kruskal|wallis") ~ "Kruskal-Wallis test",
    str_detect(m, "linear") & str_detect(m, "regression") ~ "Linear regression",
    !str_detect(m, "repeated") & str_detect(m, "(^|[^a-z])anova([^a-z]|$)|analysis of variance") ~ "One-way ANOVA",
    TRUE ~ NA_character_
  )
}

raw <- read_csv(
  unz(input_zip, "study3.csv"),
  col_types = cols(.default = col_character()),
  progress = FALSE,
  show_col_types = FALSE
) %>%
  mutate(
    across(c(EARLY_PHASE1, PHASE1, PHASE2, PHASE3, PHASE4), ~ suppressWarnings(as.integer(.x))),
    value_num = parse_number(value, locale = locale(grouping_mark = ","), na = c("", "NA")),
    p_text = str_squish(coalesce(pValue, "")),
    p_bound = suppressWarnings(as.numeric(str_extract(p_text, "[0-9]+[.]?[0-9]*(?:[eE][-+]?[0-9]+)?"))),
    p_left_censored = str_detect(p_text, "<"),
    p_right_censored = str_detect(p_text, ">"),
    interventional = str_to_upper(coalesce(studyType, "")) == "INTERVENTIONAL",
    eligible_outcome = str_to_upper(coalesce(outcomeType, "")) %in% c("PRIMARY", "SECONDARY"),
    exact_early_phase1 = coalesce(EARLY_PHASE1, 0L) == 1L & coalesce(PHASE1, 0L) == 0L &
      coalesce(PHASE2, 0L) == 0L & coalesce(PHASE3, 0L) == 0L & coalesce(PHASE4, 0L) == 0L,
    exact_phase1 = coalesce(PHASE1, 0L) == 1L & coalesce(EARLY_PHASE1, 0L) == 0L &
      coalesce(PHASE2, 0L) == 0L & coalesce(PHASE3, 0L) == 0L & coalesce(PHASE4, 0L) == 0L,
    exact_phase2 = coalesce(PHASE2, 0L) == 1L & coalesce(EARLY_PHASE1, 0L) == 0L &
      coalesce(PHASE1, 0L) == 0L & coalesce(PHASE3, 0L) == 0L & coalesce(PHASE4, 0L) == 0L
  )

# Rows sharing analysisId represent group denominators for a single statistical
# analysis; they are collapsed rather than treated as duplicate analyses.
analysis_level <- raw %>%
  group_by(analysisId) %>%
  summarise(
    nctId = first(nctId),
    briefTitle = first(briefTitle),
    outcomeTitle = first(outcomeTitle),
    outcomeType = first(outcomeType),
    conditionIds = first(condition_ids),
    originalMethod = first(statisticalMethod),
    p_reported = first(p_text),
    p_bound = first(p_bound),
    p_left_censored = first(p_left_censored),
    p_right_censored = first(p_right_censored),
    studyType = first(studyType),
    exact_early_phase1 = first(exact_early_phase1),
    exact_phase1 = first(exact_phase1),
    exact_phase2 = first(exact_phase2),
    interventional = first(interventional),
    eligible_outcome = first(eligible_outcome),
    group_text = paste(unique(na.omit(groupDescription)), collapse = " | "),
    values_complete = all(is.finite(value_num) & value_num > 0),
    group_rows = n(),
    value_first = first(value_num),
    value_sum = sum(value_num),
    value_max = max(value_num),
    .groups = "drop"
  ) %>%
  mutate(
    phase = case_when(
      exact_early_phase1 ~ "Early Phase I",
      exact_phase1 ~ "Phase I",
      exact_phase2 ~ "Phase II",
      TRUE ~ NA_character_
    ),
    test_class = classify_method(originalMethod),
    sample_size = case_when(
      test_class %in% c("One-sample/paired t-test", "Wilcoxon signed-rank test") ~ value_first,
      test_class == "Conditional logistic regression" ~ value_max,
      TRUE ~ value_sum
    ),
    dimension = case_when(
      test_class %in% c("One-way ANOVA", "Kruskal-Wallis test") ~ group_rows - 1L,
      !is.na(test_class) ~ 1L,
      TRUE ~ NA_integer_
    ),
    direction_text = str_to_lower(paste(originalMethod, group_text, outcomeTitle)),
    explicit_one_sided_method = str_detect(str_to_lower(originalMethod), "one[ -]?sided|1[ -]?sided|one[ -]?tailed|1[ -]?tailed"),
    explicit_two_sided_method = str_detect(str_to_lower(originalMethod), "two[ -]?sided|2[ -]?sided|two[ -]?tailed|2[ -]?tailed"),
    contextual_one_sided = str_detect(direction_text, "one[ -]?sided|1[ -]?sided|one[ -]?tailed|1[ -]?tailed"),
    # The Wagenmakers correction is directional and scalar. We require an
    # explicit one-sided label in the statistical-method field. Contextual
    # mentions often describe a power calculation rather than the reported
    # analysis and are therefore not used to reclassify its direction.
    one_sided = dimension == 1L & explicit_one_sided_method,
    test_direction = if_else(one_sided, "One-sided", "Two-sided/unrestricted"),
    row_structure_valid = case_when(
      test_class %in% c("One-sample/paired t-test", "Wilcoxon signed-rank test") ~ group_rows == 1L,
      test_class %in% c("One-way ANOVA", "Kruskal-Wallis test") ~ group_rows >= 2L,
      !is.na(test_class) ~ group_rows >= 1L,
      TRUE ~ FALSE
    ),
    p_eligible = is.finite(p_bound) & p_bound > 0 & p_bound < 1 & !p_right_censored,
    sample_eligible = values_complete & row_structure_valid & is.finite(sample_size) & sample_size > 1 &
      is.finite(dimension) & dimension >= 1
  )

flow_steps <- list(
  "Unique statistical analyses in archived extraction" = rep(TRUE, nrow(analysis_level)),
  "After removing four prespecified malformed analysis identifiers" = !analysis_level$analysisId %in% bad_analysis_ids,
  "Interventional studies" = analysis_level$interventional,
  "Exact Early Phase I, Phase I, or Phase II classification" = !is.na(analysis_level$phase),
  "Primary- and secondary-outcome analyses" = analysis_level$eligible_outcome,
  "Parseable p-value in (0,1), excluding right-censored reports" = analysis_level$p_eligible,
  "Method mapped to a prespecified eJAB test class" = !is.na(analysis_level$test_class),
  "Required denominator counts and test dimension available" = analysis_level$sample_eligible
)

keep <- rep(TRUE, nrow(analysis_level))
flow <- bind_rows(lapply(seq_along(flow_steps), function(i) {
  before <- sum(keep)
  keep <<- keep & flow_steps[[i]]
  tibble(
    step = i,
    criterion = names(flow_steps)[[i]],
    excluded_at_step = before - sum(keep),
    retained = sum(keep),
    studies_retained = n_distinct(analysis_level$nctId[keep])
  )
}))

eligible <- analysis_level[keep, ] %>%
  mutate(
    p_two_sided = if_else(one_sided, 2 * pmin(p_bound, 1 - p_bound), p_bound),
    ejab_two_sided_01 = sqrt(sample_size) * exp(
      -0.5 * ((sample_size^(1 / dimension) - 1) / sample_size^(1 / dimension)) *
        qchisq(1 - p_two_sided, df = dimension)
    ),
    # The scalar one-sided approximation omits the (n - 1) / n correction
    # in the exponent, before applying the directional mass correction.
    ejab_01 = if_else(
      one_sided,
      sqrt(sample_size) * exp(-0.5 * qchisq(1 - p_two_sided, df = 1)) /
        (2 * (1 - p_bound)),
      ejab_two_sided_01
    ),
    significance_threshold = if_else(one_sided, 0.025, 0.05),
    # For a report such as p < .05, the stored numerical bound is .05 but the
    # underlying p-value is strictly below it. The bound is retained for the
    # conservative eJAB calculation while the inequality determines whether
    # the direction-specific significance rule is met.
    statistically_significant = if_else(
      p_left_censored,
      p_bound <= significance_threshold,
      p_bound < significance_threshold
    ),
    favors_null = ejab_01 > 1,
    contradiction = statistically_significant & favors_null,
    evidence_category = case_when(
      !statistically_significant & !favors_null ~ "not significant / favors alternative",
      !statistically_significant & favors_null ~ "not significant / favors null",
      statistically_significant & !favors_null ~ "significant / favors alternative",
      statistically_significant & favors_null ~ "significant / favors null"
    )
  ) %>%
  filter(is.finite(ejab_01), ejab_01 > 0)

if (nrow(eligible) != tail(flow$retained, 1)) {
  stop("At least one record passed screening but did not yield a finite eJAB value.")
}

phase_summary <- eligible %>%
  group_by(phase) %>%
  summarise(
    eligible_studies = n_distinct(nctId),
    eligible_analyses = n(),
    one_sided_analyses = sum(one_sided),
    left_censored_p_values = sum(p_left_censored),
    statistically_significant = sum(statistically_significant),
    significant_favoring_null = sum(contradiction),
    raw_contradiction_proportion = significant_favoring_null / statistically_significant,
    .groups = "drop"
  )

test_summary <- eligible %>%
  group_by(phase, test_class, test_direction) %>%
  summarise(
    eligible_analyses = n(),
    statistically_significant = sum(statistically_significant),
    significant_favoring_null = sum(contradiction),
    .groups = "drop"
  )

study_categories <- eligible %>%
  count(phase, nctId, evidence_category) %>%
  complete(
    nesting(phase, nctId),
    evidence_category = c(
      "not significant / favors alternative",
      "not significant / favors null",
      "significant / favors alternative",
      "significant / favors null"
    ),
    fill = list(n = 0)
  ) %>%
  pivot_wider(names_from = evidence_category, values_from = n) %>%
  mutate(cluster_size = rowSums(across(where(is.numeric))))

draw_phase <- function(studies, draws = 20000L) {
  numerator <- numeric(draws)
  denominator <- numeric(draws)
  category_columns <- c(
    "not significant / favors alternative",
    "not significant / favors null",
    "significant / favors alternative",
    "significant / favors null"
  )
  for (i in seq_len(nrow(studies))) {
    alpha <- 1 + as.numeric(studies[i, category_columns])
    gamma_draws <- matrix(rgamma(draws * 4L, shape = rep(alpha, each = draws)), ncol = 4L)
    probability_draws <- gamma_draws / rowSums(gamma_draws)
    cluster_size <- studies$cluster_size[[i]]
    numerator <- numerator + cluster_size * probability_draws[, 4]
    denominator <- denominator + cluster_size * (probability_draws[, 3] + probability_draws[, 4])
  }
  numerator / denominator
}

set.seed(20260828)
phase_levels <- c("Early Phase I", "Phase I", "Phase II")
phase_draws <- setNames(
  lapply(phase_levels, function(phase_name) {
    draw_phase(filter(study_categories, phase == phase_name))
  }),
  phase_levels
)

summarise_draws <- function(x, parameter) {
  tibble(
    parameter = parameter,
    posterior_mean = mean(x),
    posterior_median = median(x),
    posterior_sd = sd(x),
    credible_lower_95 = unname(quantile(x, 0.025)),
    credible_upper_95 = unname(quantile(x, 0.975))
  )
}

posterior_summary <- bind_rows(lapply(phase_levels, function(phase_name) {
  summarise_draws(phase_draws[[phase_name]], paste(phase_name, "proportion"))
}))

posterior_draws <- tibble(draw = seq_along(phase_draws[[1]]))
for (phase_name in phase_levels) {
  posterior_draws[[paste(phase_name, "proportion")]] <- phase_draws[[phase_name]]
}

mixed_model_summary <- tibble(
  term = character(), estimate = double(), std_error = double(),
  odds_ratio = double(), conf_low_95 = double(), conf_high_95 = double(),
  n_analyses = integer(), n_studies = integer()
)
if (requireNamespace("lme4", quietly = TRUE)) {
  significant_data <- eligible %>%
    filter(statistically_significant) %>%
    mutate(phase = relevel(factor(phase, levels = phase_levels), ref = "Phase I"))
  mixed_fit <- lme4::glmer(
    contradiction ~ phase + (1 | nctId),
    data = significant_data,
    family = binomial(),
    control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  )
  coefficient_table <- coef(summary(mixed_fit))
  phase_terms <- grep("^phase", rownames(coefficient_table), value = TRUE)
  mixed_model_summary <- bind_rows(lapply(phase_terms, function(phase_term) {
    estimate <- coefficient_table[phase_term, "Estimate"]
    standard_error <- coefficient_table[phase_term, "Std. Error"]
    tibble(
      term = phase_term,
      estimate = estimate,
      std_error = standard_error,
      odds_ratio = exp(estimate),
      conf_low_95 = exp(estimate - 1.96 * standard_error),
      conf_high_95 = exp(estimate + 1.96 * standard_error),
      n_analyses = nrow(significant_data),
      n_studies = n_distinct(significant_data$nctId)
    )
  }))
}

write_csv(flow, file.path(output_dir, "inclusion_exclusion_flow.csv"))
write_csv(phase_summary, file.path(output_dir, "phase_summary.csv"))
write_csv(test_summary, file.path(output_dir, "test_summary.csv"))
write_csv(posterior_summary, file.path(output_dir, "posterior_summary.csv"))
write_csv(posterior_draws, file.path(output_dir, "posterior_draws.csv"))
write_csv(mixed_model_summary, file.path(output_dir, "mixed_model_summary.csv"))
write_csv(
  eligible %>% select(
    analysisId, nctId, phase, outcomeTitle, outcomeType, conditionIds, originalMethod, test_class,
    test_direction, one_sided, p_reported, p_bound, p_two_sided, p_left_censored, sample_size,
    dimension, statistically_significant, ejab_01, favors_null, contradiction
  ),
  file.path(output_dir, "eligible_early_phase_analyses.csv")
)

posterior_long <- posterior_draws %>%
  pivot_longer(-draw, names_to = "parameter", values_to = "value") %>%
  mutate(parameter = factor(
    parameter,
    levels = paste(phase_levels, "proportion")
  ))

density_figure <- ggplot(posterior_long, aes(x = value, fill = parameter, color = parameter)) +
  geom_density(alpha = 0.30, linewidth = 0.7) +
  geom_vline(data = posterior_summary,
             aes(xintercept = posterior_mean, color = parameter),
             linetype = "dashed", linewidth = 0.6, show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free", ncol = 3) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#009E73"), guide = "none") +
  scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73"), guide = "none") +
  labs(x = "Posterior value", y = "Density") +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey95"), panel.grid.minor = element_blank())

ggsave(file.path(figure_dir, "early_phase_posterior.png"), density_figure,
       width = 10.5, height = 3.5, units = "in", dpi = 400)
ggsave(file.path(response_figure_dir, "early_phase_posterior.png"), density_figure,
       width = 10.5, height = 3.5, units = "in", dpi = 400)

alpha_grid <- seq(0.001, 0.25, by = 0.001)
alpha_profile <- bind_rows(lapply(phase_levels, function(phase_name) {
  phase_frame <- filter(eligible, phase == phase_name)
  bind_rows(lapply(alpha_grid, function(alpha_level) {
    is_significant <- if_else(
      phase_frame$p_left_censored,
      phase_frame$p_two_sided <= alpha_level,
      phase_frame$p_two_sided < alpha_level
    )
    tibble(
      phase = phase_name,
      alpha = alpha_level,
      contradiction_rate = mean(is_significant & phase_frame$ejab_01 > 1)
    )
  }))
})) %>%
  mutate(phase = factor(phase, levels = phase_levels))

alpha_figure <- ggplot(
  alpha_profile,
  aes(alpha, contradiction_rate, color = phase)
) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", linewidth = 0.6) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73")) +
  scale_x_continuous(limits = c(0, 0.25), breaks = seq(0, 0.25, 0.05), expand = expansion(mult = c(0, 0.015))) +
  scale_y_continuous(limits = c(0, 0.25), breaks = seq(0, 0.25, 0.05), expand = expansion(mult = c(0, 0.015))) +
  labs(
    x = expression("Significance level " * alpha),
    y = "Contradiction rate",
    color = "Trial phase"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

ggsave(file.path(response_figure_dir, "early_phase_alpha_contradictions.png"), alpha_figure,
       width = 7.2, height = 4.5, units = "in", dpi = 400)

interval_figure <- posterior_summary %>%
  mutate(parameter = factor(parameter, levels = rev(parameter))) %>%
  ggplot(aes(x = posterior_median, y = parameter, color = parameter)) +
  geom_errorbar(aes(xmin = credible_lower_95, xmax = credible_upper_95),
                width = 0.15, linewidth = 0.8, orientation = "y") +
  geom_point(size = 2.4) +
  scale_color_manual(values = c("#009E73", "#D55E00", "#0072B2"), guide = "none") +
  labs(x = "Posterior median and 95% credible interval", y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(figure_dir, "early_phase_summary.png"), interval_figure,
       width = 7.2, height = 3.2, units = "in", dpi = 400)

cat("Analysis complete\n")
print(flow)
print(phase_summary)
print(posterior_summary)
print(mixed_model_summary)
