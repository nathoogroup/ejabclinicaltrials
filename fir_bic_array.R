args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: fir_bic_array.R TASK_ID SHARD_DIR WORK_DIR")
}

task_id <- as.integer(args[[1]])
shard_dir <- normalizePath(args[[2]], mustWork = FALSE)
work_dir <- normalizePath(args[[3]], mustWork = TRUE)
if (is.na(task_id)) stop("TASK_ID must be an integer")

setwd(work_dir)
render_dir <- tempfile(pattern = sprintf("bic-task-%04d-", task_id))
dir.create(render_dir)
on.exit(unlink(render_dir, recursive = TRUE), add = TRUE)

rmarkdown::render(
  file.path(work_dir, "bayes_risk_analysis.Rmd"),
  output_file = "task.html",
  output_dir = render_dir,
  intermediates_dir = render_dir,
  knit_root_dir = work_dir,
  params = list(
    bic_task_id = task_id,
    bic_shard_dir = shard_dir
  ),
  envir = new.env(parent = globalenv()),
  quiet = TRUE
)

shard_file <- file.path(shard_dir, sprintf("bic_task_%04d.rds", task_id))
if (!file.exists(shard_file)) stop("Expected shard was not created: ", shard_file)
message("Wrote ", shard_file)
