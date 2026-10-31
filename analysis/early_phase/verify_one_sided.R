#!/usr/bin/env Rscript
# Check generated values using an independent normal-quantile expression.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
data <- read.csv(file.path(script_dir, "outputs", "eligible_early_phase_analyses.csv"))
one <- data$one_sided
stopifnot(all(data$dimension[one] == 1), sum(one) == 427, sum(!one) == 16136)
expected_one <- sqrt(data$sample_size[one]) * exp(
  -qnorm(data$p_bound[one])^2 / 2
) / (2 * (1 - data$p_bound[one]))
expected_two <- sqrt(data$sample_size[!one]) * exp(
  -0.5 * (1 - data$sample_size[!one]^(-1 / data$dimension[!one])) *
    qchisq(1 - data$p_bound[!one], df = data$dimension[!one])
)
stopifnot(isTRUE(all.equal(data$ejab_01[one], expected_one, tolerance = 1e-10)))
stopifnot(isTRUE(all.equal(data$ejab_01[!one], expected_two, tolerance = 1e-10)))
stopifnot(identical(data$contradiction, data$statistically_significant & data$ejab_01 > 1))
two_at_reconstructed_p <- sqrt(data$sample_size[one]) * exp(
  -0.5 * (1 - 1 / data$sample_size[one]) * qchisq(1 - data$p_two_sided[one], df = 1)
)
sample_size <- data$sample_size[one]
identity_log <- log(sample_size) / (2 * (1 - sample_size)) +
  (sample_size / (sample_size - 1)) * log(two_at_reconstructed_p) -
  log(2 * (1 - data$p_bound[one]))
stopifnot(isTRUE(all.equal(log(data$ejab_01[one]), identity_log, tolerance = 1e-10)))
direct_correction <- two_at_reconstructed_p / (2 * (1 - data$p_bound[one]))
stopifnot(all(data$ejab_01[one] <= direct_correction * (1 + 1e-10)))
cat("Verified", sum(one), "one-sided and", sum(!one), "two-sided Bayes factors.\n")

# Boundary reports must be included in the manuscript's directional NHST rule.
threshold <- ifelse(one, 0.025, 0.05)
boundary <- !data$p_left_censored & data$p_bound == threshold
stopifnot(sum(boundary & one) == 2, sum(boundary & !one) == 63)
stopifnot(all(data$statistically_significant[boundary]))
stopifnot(!any(data$statistically_significant[data$p_bound > threshold]))
stopifnot(all(data$statistically_significant[data$p_bound < threshold]))
stopifnot(sum(data$statistically_significant) == 6474, sum(data$contradiction) == 493)
stopifnot(sum(data$statistically_significant & data$ejab_01 > 3) == 20)
cat("Verified inclusive thresholds: 65 boundary reports, 6474 significant results, 493 contradictions.\n")
