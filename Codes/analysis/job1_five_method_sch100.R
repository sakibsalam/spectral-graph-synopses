## =====================================================================================
## Job 1 (this extension's headline): five-method comparison on nbk_0v1 at V=100
## (Schaefer-100, p = 4,950). Builds on Codes/analysis/multi_resolution_common.R (which
## sources outcome_explorer.R -> abide_classification.R -> abide_pipeline.R ->
## network_bayes_simulation.R, giving access to build_features(), graph_spectral_ard_logit(),
## graph_spectral_plus_sparse_logit() [hybrid], raw_edge_glmnet_logit() [ridge/lasso],
## raw_edge_ard(), spike_slab_ssvs(), ard_fit(), standardize_train_test(), auc_pred()).
##
## Arms and split counts (per the task spec):
##   gs_ard, hybrid, raw_edge_ridge, raw_edge_lasso, raw_edge_ssvs  -- 10 splits each.
##   raw_edge_ard                                                   -- 3 splits ONLY, with
##     max_iter CAPPED (feasibility-checked at 25 first; if a single fit at cap=25 exceeds
##     45 minutes, drops to cap=10 for all 3 splits and records that decision).
##
## raw_edge_ard() in network_bayes_simulation.R has no max_iter argument (hardcoded to
## ard_fit()'s default of 500) -- this script defines a capped wrapper that exposes it and
## additionally reports ard_fit()'s own `iter` field (the achieved iteration count at
## break/exhaustion) so a capped, non-converged fit is never confused with a converged one.
##
## raw_edge_ard / spike_slab_ssvs are Gaussian-likelihood regression fits (see
## network_bayes_simulation.R header) applied here to y coded {0,1} exactly as elsewhere in
## this project (e.g. five_method_multires.R) -- valid for AUC ranking (auc_pred() only needs
## a rank-consistent score) but NOT a calibrated probability; accuracy-at-0.5 is recorded only
## as a secondary, exploratory diagnostic.
##
## Resumable: skips any (outcome, split, method) row already present in out_csv.
## Output: results/five_method_sch100.csv (long, one row per split x method) plus
## results/five_method_sch100_summary.csv (means per method, split counts and timing noted).
## =====================================================================================

suppressMessages(source("Codes/analysis/multi_resolution_common.R"))

## Gate on the capped raw_edge_ard arm (Stage 4 of the 2026-09-11 resume task): this arm alone
## costs ~2h even capped, and the task instructions require deciding whether to run it at all
## AFTER seeing Stages 1-3 results, not automatically as part of finishing the 5 main arms to
## 10 splits (Stage 3). Default OFF; set RUN_ARD_ARM=1 in the environment to include it.
RUN_ARD_ARM <- identical(Sys.getenv("RUN_ARD_ARM", "0"), "1")

OUT_CSV     <- "results/five_method_sch100.csv"
SUMMARY_CSV <- "results/five_method_sch100_summary.csv"
NAME        <- "nbk_0v1"
N_SPLITS_MAIN <- 10
N_SPLITS_ARD  <- 3
SEED0 <- 900
TEST_FRAC <- 0.30
K <- 10
BASIS <- "adjacency"
ARD_MAX_ITER_TRY <- 25
ARD_MAX_ITER_FALLBACK <- 10
ARD_TIME_BUDGET_SEC <- 45 * 60   # 45 minutes

raw_edge_ard_capped <- function(X, y, Xt, max_iter) {
  t0 <- proc.time()[3]
  pp <- standardize_train_test(X, Xt, y)
  fit <- ard_fit(pp$X, pp$y, max_iter = max_iter)
  pred <- as.numeric(pp$Xt %*% fit$coef + pp$ybar)
  list(coef = fit$coef, pred = pred, seconds = proc.time()[3] - t0,
       iter_achieved = fit$iter, iter_cap = max_iter,
       converged_below_cap = fit$iter < max_iter)
}

m <- readRDS("data/osf_schaefer100/master.rds")
cat(sprintf("Loaded Schaefer-100 master: %d connectomes, V=%d, %d subjects\n",
            length(m$conn), m$V, length(unique(m$subjects$SUBJECT))))

set <- make_set_m(m, conds = c("nbk_0back", "nbk_1back"))
conn <- set$conn; meta <- set$meta
y <- as.integer(meta$COND == "nbk_1back")
subj <- meta$SUBJECT; usub <- unique(subj)
cat(sprintf("[%s] N=%d rows, %d subjects\n", NAME, length(conn), length(usub)))

already <- data.frame()
if (file.exists(OUT_CSV)) already <- read.csv(OUT_CSV)
is_done <- function(split, method) {
  nrow(already) > 0 && any(already$outcome == NAME & already$split == split & already$method == method)
}
append_row <- function(row) {
  prev <- if (file.exists(OUT_CSV)) read.csv(OUT_CSV) else NULL
  dir.create(dirname(OUT_CSV), showWarnings = FALSE, recursive = TRUE)
  write.csv(rbind(prev, row), OUT_CSV, row.names = FALSE)
  already <<- if (is.null(prev)) row else rbind(prev, row)
}
mk_row <- function(split, method, auc, seconds, n_nonzero, pred_sd, acc,
                    iter_cap = NA_real_, iter_achieved = NA_real_, converged = NA) {
  data.frame(outcome = NAME, split = split, method = method, auc = auc, acc = acc,
             seconds = seconds, n_nonzero = n_nonzero, pred_sd = pred_sd,
             iter_cap = iter_cap, iter_achieved = iter_achieved, converged = converged,
             p_or_q = NA_integer_, n_train = NA_integer_)
}

## ---- main arms: gs_ard, hybrid, ridge, lasso, ssvs -- 10 splits ---------------------
for (s in seq_len(N_SPLITS_MAIN)) {
  set.seed(SEED0 + s)
  te_s <- sample(usub, max(1, round(TEST_FRAC * length(usub))))
  te <- which(subj %in% te_s); tr <- which(!(subj %in% te_s))
  if (length(unique(y[tr])) < 2 || length(unique(y[te])) < 2) {
    cat(sprintf("[%s] split %d: degenerate label split -- skipping\n", NAME, s)); next
  }
  need_any <- !is_done(s, "gs_ard") || !is_done(s, "hybrid") || !is_done(s, "raw_edge_ridge") ||
              !is_done(s, "raw_edge_lasso") || !is_done(s, "raw_edge_ssvs")
  if (!need_any) { cat(sprintf("[%s] split %d main arms already done -- skipping\n", NAME, s)); next }

  ft <- build_features(conn, tr, K, basis = BASIS)
  Xs <- ft$X_spec; Xe <- ft$X_edge
  y_tr <- y[tr]; y_te <- y[te]
  n_train <- length(tr)

  if (!is_done(s, "gs_ard")) {
    f <- graph_spectral_ard_logit(Xs[tr, ], y_tr, Xs[te, ])
    r <- mk_row(s, "gs_ard", auc_pred(f$prob, y_te), f$seconds, f$nsel, sd(f$prob),
                mean((f$prob > 0.5) == y_te))
    r$p_or_q <- ncol(Xs); r$n_train <- n_train
    append_row(r)
    cat(sprintf("[%s] split %d gs_ard: auc=%.3f sec=%.2f nsel=%d\n", NAME, s, r$auc, r$seconds, r$n_nonzero))
  }
  if (!is_done(s, "hybrid")) {
    f <- graph_spectral_plus_sparse_logit(Xs[tr, ], Xe[tr, ], y_tr, Xs[te, ], Xe[te, ])
    r <- mk_row(s, "hybrid", auc_pred(f$prob, y_te), f$seconds, NA_integer_, sd(f$prob),
                mean((f$prob > 0.5) == y_te))
    r$p_or_q <- ncol(Xe); r$n_train <- n_train
    append_row(r)
    cat(sprintf("[%s] split %d hybrid: auc=%.3f sec=%.2f (n_nonzero not returned by graph_spectral_plus_sparse_logit)\n",
                NAME, s, r$auc, r$seconds))
  }
  if (!is_done(s, "raw_edge_ridge")) {
    f <- raw_edge_glmnet_logit(Xe[tr, ], y_tr, Xe[te, ], alpha = 0)
    r <- mk_row(s, "raw_edge_ridge", auc_pred(f$prob, y_te), f$seconds, f$nsel, sd(f$prob),
                mean((f$prob > 0.5) == y_te))
    r$p_or_q <- ncol(Xe); r$n_train <- n_train
    append_row(r)
    cat(sprintf("[%s] split %d raw_edge_ridge: auc=%.3f sec=%.2f\n", NAME, s, r$auc, r$seconds))
  }
  if (!is_done(s, "raw_edge_lasso")) {
    f <- raw_edge_glmnet_logit(Xe[tr, ], y_tr, Xe[te, ], alpha = 1)
    r <- mk_row(s, "raw_edge_lasso", auc_pred(f$prob, y_te), f$seconds, f$nsel, sd(f$prob),
                mean((f$prob > 0.5) == y_te))
    r$p_or_q <- ncol(Xe); r$n_train <- n_train
    append_row(r)
    cat(sprintf("[%s] split %d raw_edge_lasso: auc=%.3f sec=%.2f nsel=%d\n", NAME, s, r$auc, r$seconds, r$n_nonzero))
  }
  if (!is_done(s, "raw_edge_ssvs")) {
    f <- spike_slab_ssvs(Xe[tr, ], y_tr, Xe[te, ])
    r <- mk_row(s, "raw_edge_ssvs", auc_pred(f$pred, y_te), f$seconds, sum(f$pip > 0.5), sd(f$pred),
                mean((f$pred > 0.5) == y_te))
    r$p_or_q <- ncol(Xe); r$n_train <- n_train
    append_row(r)
    cat(sprintf("[%s] split %d raw_edge_ssvs: auc=%.3f sec=%.2f n_pip>0.5=%d\n", NAME, s, r$auc, r$seconds, r$n_nonzero))
  }
  flush.console()
}

## ---- raw_edge_ard: 3 splits, capped, with an explicit feasibility check first --------
## Gated by RUN_ARD_ARM (see top of file) -- skipped by default.
if (!RUN_ARD_ARM) {
  cat("\n[raw_edge_ard] RUN_ARD_ARM not set -- skipping the capped raw-edge ARD arm (Stage 4).\n")
}
ard_cap <- ARD_MAX_ITER_TRY
first_ard_split <- NA_integer_
if (RUN_ARD_ARM) {
  for (s in seq_len(N_SPLITS_ARD)) if (!is_done(s, "raw_edge_ard")) { first_ard_split <- s; break }
}

if (!is.na(first_ard_split)) {
  s <- first_ard_split
  set.seed(SEED0 + s)
  te_s <- sample(usub, max(1, round(TEST_FRAC * length(usub))))
  te <- which(subj %in% te_s); tr <- which(!(subj %in% te_s))
  ft <- build_features(conn, tr, K, basis = BASIS)
  Xe <- ft$X_edge
  cat(sprintf("\n[%s] FEASIBILITY CHECK: timing raw_edge_ard at max_iter=%d, p=%d, n_train=%d (split %d)...\n",
              NAME, ARD_MAX_ITER_TRY, ncol(Xe), length(tr), s))
  f <- raw_edge_ard_capped(Xe[tr, ], y[tr], Xe[te, ], max_iter = ARD_MAX_ITER_TRY)
  cat(sprintf("  TIMING: raw_edge_ard(max_iter=%d) = %.1fs (iter_achieved=%d)\n",
              ARD_MAX_ITER_TRY, f$seconds, f$iter_achieved))
  if (f$seconds > ARD_TIME_BUDGET_SEC) {
    ard_cap <- ARD_MAX_ITER_FALLBACK
    cat(sprintf("  ** exceeds the %d-minute budget -- DROPPING cap to max_iter=%d for all %d splits **\n",
                ARD_TIME_BUDGET_SEC / 60, ARD_MAX_ITER_FALLBACK, N_SPLITS_ARD))
    ## the split-1 fit above used the wrong cap; redo it at the fallback cap so all 3 splits
    ## are computed under the SAME cap and are directly comparable
    f <- raw_edge_ard_capped(Xe[tr, ], y[tr], Xe[te, ], max_iter = ard_cap)
  } else {
    cat(sprintf("  within budget -- proceeding with max_iter=%d for all %d splits\n", ard_cap, N_SPLITS_ARD))
  }
  y_te <- y[te]
  r <- mk_row(s, "raw_edge_ard", auc_pred(f$pred, y_te), f$seconds, sum(abs(f$coef) > 1e-6), sd(f$pred),
              mean((f$pred > 0.5) == y_te), iter_cap = ard_cap, iter_achieved = f$iter_achieved,
              converged = f$converged_below_cap)
  r$p_or_q <- ncol(Xe); r$n_train <- length(tr)
  append_row(r)
  cat(sprintf("[%s] split %d raw_edge_ard (CAPPED, NOT CONVERGED unless converged=TRUE): auc=%.3f sec=%.1f iter=%d/%d\n",
              NAME, s, r$auc, r$seconds, r$iter_achieved, r$iter_cap))
}

for (s in if (RUN_ARD_ARM) seq_len(N_SPLITS_ARD) else integer(0)) {
  if (is_done(s, "raw_edge_ard") || s == first_ard_split) next
  set.seed(SEED0 + s)
  te_s <- sample(usub, max(1, round(TEST_FRAC * length(usub))))
  te <- which(subj %in% te_s); tr <- which(!(subj %in% te_s))
  ft <- build_features(conn, tr, K, basis = BASIS)
  Xe <- ft$X_edge
  y_tr <- y[tr]; y_te <- y[te]
  f <- raw_edge_ard_capped(Xe[tr, ], y_tr, Xe[te, ], max_iter = ard_cap)
  r <- mk_row(s, "raw_edge_ard", auc_pred(f$pred, y_te), f$seconds, sum(abs(f$coef) > 1e-6), sd(f$pred),
              mean((f$pred > 0.5) == y_te), iter_cap = ard_cap, iter_achieved = f$iter_achieved,
              converged = f$converged_below_cap)
  r$p_or_q <- ncol(Xe); r$n_train <- length(tr)
  append_row(r)
  cat(sprintf("[%s] split %d raw_edge_ard (CAPPED, NOT CONVERGED unless converged=TRUE): auc=%.3f sec=%.1f iter=%d/%d\n",
              NAME, s, r$auc, r$seconds, r$iter_achieved, r$iter_cap))
  flush.console()
}

## ---- summary -------------------------------------------------------------------------
d <- read.csv(OUT_CSV)
agg <- aggregate(cbind(auc, seconds) ~ method, d, mean)
n_splits_tab <- aggregate(split ~ method, d, length); names(n_splits_tab)[2] <- "n_splits"
agg <- merge(agg, n_splits_tab, by = "method")
ard_note <- d[d$method == "raw_edge_ard", c("iter_cap", "iter_achieved", "converged")]
agg$note <- ""
if (nrow(ard_note) > 0) {
  agg$note[agg$method == "raw_edge_ard"] <- sprintf(
    "CAPPED at max_iter=%s (NOT a converged fit); achieved iters: %s",
    paste(unique(ard_note$iter_cap), collapse = "/"), paste(ard_note$iter_achieved, collapse = ","))
}
agg <- agg[order(-agg$auc), ]
write.csv(agg, SUMMARY_CSV, row.names = FALSE)
cat("\n==== Job 1 summary ====\n"); print(agg)
cat("\nWrote", OUT_CSV, "and", SUMMARY_CSV, "\n")
