#!/usr/bin/env Rscript
# Check generated values using an independent normal-quantile expression.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
data <- read.csv(file.path(script_dir, "outputs", "eligible_early_phase_analyses.csv"))
one <- data$one_sided
stopifnot(all(data$dimension[one] == 1), sum(one) == 427, sum(!one) == 16136)
expected_one <- sqrt(data$sample_size[one]) * exp(-qnorm(data$p_bound[one])^2 / 2) /
  (2 * (1 - data$p_bound[one]))
expected_two <- sqrt(data$sample_size[!one]) * exp(
  -0.5 * (1 - data$sample_size[!one]^(-1 / data$dimension[!one])) *
    qchisq(1 - data$p_bound[!one], df = data$dimension[!one])
)
stopifnot(isTRUE(all.equal(data$ejab_01[one], expected_one, tolerance = 1e-10)))
stopifnot(isTRUE(all.equal(data$ejab_01[!one], expected_two, tolerance = 1e-10)))
stopifnot(identical(data$contradiction, data$statistically_significant & data$ejab_01 > 1))
cat("Verified", sum(one), "one-sided and", sum(!one), "two-sided Bayes factors.\n")
