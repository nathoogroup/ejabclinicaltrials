args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: fir_bic_merge_render.R SHARD_DIR WORK_DIR")
}

shard_dir <- normalizePath(args[[1]], mustWork = TRUE)
work_dir <- normalizePath(args[[2]], mustWork = TRUE)
setwd(work_dir)

bic_test_keys <- c("ttest", "lm", "logistic", "anova", "ranova", "chisq")
bic_n_grid <- unique(round(exp(seq(log(30), log(1e7), length.out = 70))))
bic_theta_grid <- unique(c(
  exp(seq(log(0.001), log(0.02), length.out = 14)),
  seq(0.025, 0.9, length.out = 46),
  seq(0.925, 2, by = 0.025)
))
task_grid <- tidyr::crossing(test_key = bic_test_keys, n = bic_n_grid) |>
  dplyr::mutate(task_id = dplyr::row_number(), .before = 1)

expected_files <- file.path(
  shard_dir, sprintf("bic_task_%04d.rds", task_grid$task_id)
)
missing_files <- expected_files[!file.exists(expected_files)]
if (length(missing_files)) {
  stop("Missing ", length(missing_files), " BIC shards; first missing: ",
       basename(missing_files[[1]]))
}

bic_risk <- dplyr::bind_rows(lapply(expected_files, readRDS))
expected_keys <- tidyr::crossing(
  test_key = bic_test_keys,
  n = bic_n_grid,
  theta = bic_theta_grid
)

if (nrow(bic_risk) != nrow(expected_keys)) {
  stop("Merged row count is ", nrow(bic_risk),
       "; expected ", nrow(expected_keys))
}
if (anyDuplicated(bic_risk[c("test_key", "n", "theta")])) {
  stop("Merged BIC keys are not unique")
}
observed_keys <- dplyr::arrange(bic_risk, test_key, n, theta) |>
  dplyr::select(test_key, n, theta)
expected_keys <- dplyr::arrange(expected_keys, test_key, n, theta)
if (!isTRUE(all.equal(observed_keys, expected_keys, check.attributes = FALSE))) {
  stop("Merged (test_key, n, theta) keys are incomplete or unexpected")
}
if (any(!is.finite(bic_risk$risk_bic))) {
  stop("Merged BIC risks contain non-finite values")
}

bic_risk <- dplyr::select(bic_risk, -task_id)
output_dir <- file.path("bayes_risk_outputs", "revised_three_comparisons")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cache_file <- file.path(output_dir, "large_sample_empirical_bic_simulation.rds")
temporary_file <- paste0(cache_file, ".tmp-", Sys.getpid())
saveRDS(bic_risk, temporary_file)
if (!file.rename(temporary_file, cache_file)) {
  unlink(temporary_file)
  stop("Could not atomically publish merged BIC cache")
}

rmarkdown::render(
  "bayes_risk_analysis.Rmd",
  output_file = "bayes_risk_analysis.html",
  params = list(refresh = TRUE, dense_refresh = TRUE, bic_refresh = FALSE),
  envir = new.env(parent = globalenv()),
  quiet = FALSE
)
