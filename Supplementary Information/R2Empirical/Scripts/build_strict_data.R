## Build ORIGINAL and REPLICATION (p, n) series for the D_n figure under the
## strict rule. Adapted from R2Empirical/Scripts/build_replication_data.R
## (3v4n-5/documents).
##
## Rule:
##   1. original p <= ALPHA
##   2. EVERY replication attempt of the result with a recorded p has p <= ALPHA
##      (all attempts, not any); a result with no replication p fails
##   3. only p and n are required, for original and replication. No q / df
##      condition: originals are read from the raw files, not the *_cleaned.csv
##      outputs, which carry the q filters (rpp.R complete cases on T_df1..O.,
##      score.R q >= 1).
##
## p-values recorded as 0 count as p <= ALPHA. For plotting, log p is recovered
## from the test statistic (same recovery as the original script, now applied to
## originals too); a zero p with no usable statistic is dropped as a point.
##
##   Rscript Scripts/build_strict_data.R [RAW_DIR]
## Output: Data/{EERP,RPP,RPCB,SCORE}_orig.csv, *_repl.csv
##         (id, p, logp, n, censored), Data/counts.csv
##
## censored = the recorded p is a reporting bound ("p < .001"), not a value:
##   SCORE originals: orig_p_value_type is not "exact" ("less-than",
##     "threshold", or no type).
##   SCORE replications (no type field): p is a power of ten and its own test
##     statistic gives less than half of it, or there is no statistic to check.
##     A t statistic with no df is checked against the normal tail, which is a
##     lower bound on the t tail for every df.
##   EERP: never. The four powers of ten among the plotted values are values,
##     per the project's replication reports: studies 6 and 11 original
##     "The original p-value is 0.001"; study 8 original "the exact p-value
##     based on the posted data which is p = 0.0001"; study 16 replication
##     "a p-value of 0.001".
##   RPP and RPCB: never; their p is computed from the test statistic.
##
## Excluded from the H1 subset: RPP study 69 (Tracy & Robins). Its p-values of
## 1e-06 (T_pval_USE, original and replication) have no test value behind them;
## the reported values are "<.05" (original) and 0.05 (replication).
EXCLUDE <- list(RPP = 69)

ALPHA <- 0.005

args    <- commandArgs(trailingOnly = TRUE)
RAW_DIR <- if (length(args) >= 1) args[1] else "/Volumes/STORAGE/Research/AhabsHarpoon/Data"
HERE    <- dirname(normalizePath(sub("^--file=", "",
             grep("^--file=", commandArgs(FALSE), value = TRUE))))
OUT_DIR <- file.path(dirname(HERE), "Data")
num <- function(x) suppressWarnings(as.numeric(x))

## --------------------------------------------------------------------------
## p-value recovery (unchanged from build_replication_data.R)
## --------------------------------------------------------------------------
recover_logp <- function(stat_type, stat, df1, df2) {
  st <- tolower(trimws(as.character(stat_type)))
  x  <- num(stat); d1 <- num(df1); d2 <- num(df2)
  out <- rep(NA_real_, length(x))
  for (i in seq_along(x)) {
    if (is.na(x[i]) || is.na(st[i])) next
    out[i] <- tryCatch({
      if (st[i] == "z") {
        log(2) + pnorm(-abs(x[i]), log.p = TRUE)
      } else if (st[i] == "t") {
        if (is.na(d1[i]) || d1[i] < 1) NA_real_
        else log(2) + pt(-abs(x[i]), df = d1[i], log.p = TRUE)
      } else if (st[i] %in% c("chi_squared", "chi-squared", "chisq",
                              "chi", "chi2")) {
        if (is.na(d1[i]) || d1[i] < 1) NA_real_
        else pchisq(x[i], df = d1[i], lower.tail = FALSE, log.p = TRUE)
      } else if (st[i] == "f") {
        if (is.na(d1[i]) || d1[i] < 1) NA_real_
        else {
          lp <- if (!is.na(d2[i]) && d2[i] >= 1)
            suppressWarnings(pf(x[i], df1 = d1[i], df2 = d2[i],
                                lower.tail = FALSE, log.p = TRUE)) else NA_real_
          ## pf underflows for a large F on a large df2; df1*F -> chi2_df1
          if (!is.finite(lp))
            lp <- pchisq(d1[i] * x[i], df = d1[i], lower.tail = FALSE, log.p = TRUE)
          lp
        }
      } else NA_real_
    }, error = function(e) NA_real_)
  }
  out
}

is_pow10 <- function(p) !is.na(p) & p > 0 & abs(log10(p) - round(log10(p))) < 1e-9

## log p for plotting: log(p) where 0 < p < 1, recovered where p == 0
with_logp <- function(p, stat_type = NA, stat = NA, df1 = NA, df2 = NA) {
  lp <- ifelse(p > 0 & p < 1, log(p), NA_real_)
  zero <- which(p == 0)
  if (length(zero)) {
    rec <- function(v) if (length(v) == 1) rep(v, length(p)) else v
    lp[zero] <- recover_logp(rec(stat_type)[zero], rec(stat)[zero],
                             rec(df1)[zero], rec(df2)[zero])
  }
  lp
}

## --------------------------------------------------------------------------
## One corpus. o: one row per result (id, p, n, logp).
##             r: one row per replication attempt (id, p, n, logp).
## --------------------------------------------------------------------------
counts <- list()
build <- function(name, o, r) {
  if (is.null(o$censored)) o$censored <- FALSE
  if (is.null(r$censored)) r$censored <- FALSE
  o <- o[!is.na(o$p) & !is.na(o$n), ]                    # denominator
  r <- r[r$id %in% o$id, ]
  rp <- r[!is.na(r$p), ]                                  # attempts with a p
  all_pass <- tapply(rp$p <= ALPHA, rp$id, all)
  has_rep  <- o$id %in% rp$id
  pass_ids <- o$id[o$p <= ALPHA & has_rep & o$id %in% names(all_pass)[all_pass]]
  excluded <- intersect(pass_ids, EXCLUDE[[name]])
  pass_ids <- setdiff(pass_ids, excluded)

  op <- o[o$id %in% pass_ids, ]
  rr <- rp[rp$id %in% pass_ids, ]
  no_n   <- sum(is.na(rr$n))
  n_le1  <- sum(!is.na(rr$n) & rr$n <= 1)
  r_zero <- sum(rr$p == 0); r_unrec <- sum(rr$p == 0 & !is.finite(rr$logp))
  o_zero <- sum(op$p == 0); o_unrec <- sum(op$p == 0 & !is.finite(op$logp))
  rr_pts <- rr[!is.na(rr$n) & rr$n > 1 & is.finite(rr$logp), ]
  op_pts <- op[is.finite(op$logp), ]

  write.csv(op_pts[order(op_pts$id), c("id", "p", "logp", "n", "censored")],
            file.path(OUT_DIR, paste0(name, "_orig.csv")), row.names = FALSE)
  write.csv(rr_pts[order(rr_pts$id), c("id", "p", "logp", "n", "censored")],
            file.path(OUT_DIR, paste0(name, "_repl.csv")), row.names = FALSE)

  ## attempts with no p, attached to results that pass on the others
  no_p_att <- sum(is.na(r$p) & r$id %in% pass_ids)
  counts[[name]] <<- data.frame(
    project = name,
    results_p_n = nrow(o),
    orig_p_le_alpha = sum(o$p <= ALPHA),
    no_repl_p = sum(o$p <= ALPHA & !has_rep),
    excluded = length(excluded),
    originals_passing = length(pass_ids),
    attempts_of_passing = nrow(rr),
    repl_dropped_no_n = no_n,
    repl_dropped_n_le1 = n_le1,
    repl_zero_p = r_zero,
    repl_dropped_zero_unrecoverable = r_unrec,
    replication_points = nrow(rr_pts),
    orig_zero_p = o_zero,
    orig_dropped_zero_unrecoverable = o_unrec,
    original_points = nrow(op_pts),
    passing_with_attempt_missing_p = no_p_att,
    orig_points_censored = sum(op_pts$censored),
    repl_points_censored = sum(rr_pts$censored))
  cat(sprintf("%-5s results %3d | orig p<=%.3f %3d | passing %3d | repl attempts %3d (no n %d, n<=1 %d, zero p %d unrecov %d) -> %3d points | orig zero p %d unrecov %d\n",
              name, nrow(o), ALPHA, sum(o$p <= ALPHA), length(pass_ids), nrow(rr),
              no_n, n_le1, r_zero, r_unrec, nrow(rr_pts), o_zero, o_unrec))
}

## --------------------------------------------------------------------------
## EERP -- one replication per study (camerer2016_studydetails.csv)
## --------------------------------------------------------------------------
e <- read.csv(file.path(RAW_DIR, "EERP", "camerer2016_studydetails.csv"))
build("EERP",
      data.frame(id = e$study, p = num(e$porig), n = num(e$norig),
                 logp = with_logp(num(e$porig))),
      data.frame(id = e$study, p = num(e$prep), n = num(e$nrep_act),
                 logp = with_logp(num(e$prep))))

## --------------------------------------------------------------------------
## RPP -- one replication per study
## --------------------------------------------------------------------------
r <- read.csv(file.path(RAW_DIR, "RPP", "rpp_data_updates.csv"))
r <- r[!is.na(r$Study.Num), ]
po <- num(r$T_pval_USE..O.); pr <- num(r$T_pval_USE..R.)
build("RPP",
      data.frame(id = r$Study.Num, p = po, n = num(r$T_N..O.),
                 logp = with_logp(po, r$T_Test.Statistic..O., r$T_Test.value..O.,
                                  r$T_df1..O., r$T_df2..O.)),
      data.frame(id = r$Study.Num, p = pr, n = num(r$T_N..R.),
                 logp = with_logp(pr, r$T_Test.Statistic..R., r$T_Test.value..R.,
                                  r$T_df1..R., r$T_df2..R.)))

## --------------------------------------------------------------------------
## RPCB -- one row per replication attempt; original (p, n) repeat within effect
## --------------------------------------------------------------------------
c2 <- read.csv(file.path(RAW_DIR, "RPCB", "rpcb_effect_level.csv"), check.names = TRUE)
id <- paste(c2$Paper.., c2$Experiment.., c2$Effect.., sep = "-")
po <- num(c2$Original.p.value); pr <- num(c2$Replication.p.value)
orig <- data.frame(id = id, p = po, n = num(c2$Original.sample.size),
                   logp = with_logp(po, c2$Original.test.statistic.type,
                                    c2$Original.test.statistic.value,
                                    c2$Original.df1, c2$Original.df2))
orig <- orig[order(is.na(orig$p) | is.na(orig$n)), ]     # prefer a complete row
build("RPCB", orig[!duplicated(orig$id), ],
      data.frame(id = id, p = pr, n = num(c2$Replication.sample.size),
                 logp = with_logp(pr, c2$Replication.test.statistic.type,
                                  c2$Replication.test.statistic.value,
                                  c2$Replication.df1, c2$Replication.df2)))

## --------------------------------------------------------------------------
## SCORE -- claims with a replication report; one row per report.
## Original n = effective sample size where given, else the raw sample size;
## this reproduces score_data.csv for all 287 of its claims.
## --------------------------------------------------------------------------
load(file.path(RAW_DIR, "SCORE", "raw", "analyst data.RData"))
oo <- as.data.frame(orig_outcomes[orig_outcomes$claim_id %in% repli_outcomes$claim_id, ])
ro <- as.data.frame(repli_outcomes)
po <- num(oo$orig_p_value); pr <- num(ro$repli_p_value)
## replication p implied by its own test statistic, to tell a rounded value
## (0.000575 recorded as 0.001) from a bound (2e-10 recorded as 0.001)
pr_stat <- exp(recover_logp(ro$repli_stat_type, ro$repli_stat_value,
                            ro$repli_stat_dof_1, ro$repli_stat_dof_2))
## a t statistic with no df: the normal tail is a lower bound on p for every df
t_no_df <- is.na(pr_stat) & tolower(trimws(ro$repli_stat_type)) %in% "t" &
  !is.na(num(ro$repli_stat_value))
pr_stat[t_no_df] <- 2 * pnorm(-abs(num(ro$repli_stat_value[t_no_df])))
build("SCORE",
      data.frame(id = oo$claim_id, p = po,
                 n = ifelse(!is.na(oo$original_effective_sample_size),
                            oo$original_effective_sample_size, oo$orig_sample_size_value),
                 logp = with_logp(po, oo$orig_stat_type, oo$orig_stat_value,
                                  oo$orig_stat_dof_1, oo$orig_stat_dof_2),
                 censored = !(oo$orig_p_value_type %in% "exact")),
      data.frame(id = ro$claim_id, p = pr, n = num(ro$repli_sample_size_value),
                 logp = with_logp(pr, ro$repli_stat_type, ro$repli_stat_value,
                                  ro$repli_stat_dof_1, ro$repli_stat_dof_2),
                 censored = is_pow10(pr) & (is.na(pr_stat) | pr_stat < pr / 2)))

out <- do.call(rbind, counts)
write.csv(out, file.path(OUT_DIR, "counts.csv"), row.names = FALSE)
