## =====================================================================================
## Job 4b: the five-method comparison at V=100 (Schaefer-100) -- the new capability this
## atlas unlocks. Raw-edge ARD and raw-edge spike-and-slab SSVS (Codes/core/
## network_bayes_simulation.R::raw_edge_ard / spike_slab_ssvs) are documented throughout this
## project (abide_pipeline.R) as computationally prohibitive at the existing real
## dataset's edge dimensions (p=19,900 at ABIDE cc200, p=26,106 at OSF Power-229) and have
## therefore NEVER been run on real data in this project. At V=100, p = 100*99/2 = 4,950 --
## small enough that both become feasible, so this is the first time the full five-method
## comparison from the simulations (graph-spectral ARD, raw-edge ridge, raw-edge lasso,
## raw-edge ARD, raw-edge spike-slab SSVS) can run on a real dataset.
##
## raw_edge_ard()/spike_slab_ssvs() are Gaussian-likelihood regression functions (see
## network_bayes_simulation.R header) with no logistic counterpart in this codebase. For the
## binary outcome (nbk_0v1) they are run on y coded {0,1} exactly as the regression outcome
## is coded -- i.e. as a LINEAR/shrinkage regression of a 0/1 target, not a calibrated
## probability model. This is scientifically valid for ranking-based AUC (auc_pred() only
## needs the continuous score to rank cases above controls; it never assumes the score is a
## calibrated probability) but the resulting scores are NOT probabilities, so "accuracy" at a
## 0.5 cutoff is reported only as a secondary, exploratory diagnostic -- AUC is the outcome
## that is directly comparable to every other binary result in this project.
##
## Feasibility guard (explicitly required by the task): before committing to the full
## n_splits run, time ONE fit of raw_edge_ard and spike_slab_ssvs at this p on the first
## split. If either exceeds max_seconds_per_fit (default 600s = 10 min), the run proceeds with
## FEWER splits (n_splits_fallback) rather than being abandoned, and the timing is recorded in
## the output CSV's `seconds` column either way so the report can quote it directly.
##
## Resumable: appends one row per (outcome, split, method) to out_csv; already-written splits
## are skipped on rerun.
## =====================================================================================

suppressMessages({
  source("Codes/analysis/multi_resolution_common.R")   # -> outcome_explorer.R -> abide_classification.R
                                                          # -> abide_pipeline.R -> network_bayes_simulation.R
                                                          # brings in: build_features, graph_spectral_ard[_logit],
                                                          # raw_edge_ridge_woodbury, raw_edge_lasso,
                                                          # raw_edge_glmnet_logit, raw_edge_ard, spike_slab_ssvs,
                                                          # auc_pred, rmse, make_set_m
})

build_outcome_5m <- function(m, name) {
  if (name == "nbk_0v1") {
    set <- make_set_m(m, conds = c("nbk_0back", "nbk_1back"))
    y <- as.integer(set$meta$COND == "nbk_1back")
    type <- "binary"
  } else if (name == "age_continuous") {
    set <- make_set_m(m, conds = NULL, aggregate_by_subject = TRUE)
    y <- as.numeric(set$meta$AGE)
    type <- "regression"
  } else stop("unknown outcome for five-method comparison: ", name)
  list(conn = set$conn, meta = set$meta, y = y, type = type)
}

## One split, all five methods. Returns a long data.frame (one row per method).
eval_split_5m <- function(conn, y, type, tr, te, K = 10, basis = "adjacency") {
  ft <- build_features(conn, tr, K, basis = basis)
  Xs <- ft$X_spec; Xe <- ft$X_edge
  y_tr <- y[tr]; y_te <- y[te]
  rows <- list()
  add <- function(method, pred_or_prob, seconds, extra = NA_real_) {
    if (type == "regression") {
      rows[[length(rows) + 1]] <<- data.frame(
        method = method, rmse = rmse(y_te, pred_or_prob),
        spearman = suppressWarnings(cor(y_te, pred_or_prob, method = "spearman")),
        auc = NA_real_, acc = NA_real_, seconds = seconds, n_selected = extra, q_or_p = ncol(if (grepl("edge", method) || method %in% c("raw_edge_ard","raw_edge_ssvs","raw_edge_ridge","raw_edge_lasso")) Xe else Xs)
      )
    } else {
      a <- auc_pred(pred_or_prob, y_te)
      rows[[length(rows) + 1]] <<- data.frame(
        method = method, rmse = NA_real_, spearman = NA_real_,
        auc = a, acc = mean((pred_or_prob > 0.5) == y_te), seconds = seconds, n_selected = extra,
        q_or_p = ncol(if (method %in% c("raw_edge_ard", "raw_edge_ssvs", "raw_edge_ridge", "raw_edge_lasso")) Xe else Xs)
      )
    }
  }

  if (type == "regression") {
    f <- graph_spectral_ard(Xs[tr, ], y_tr, Xs[te, ]);           add("graph_spectral_ard", f$pred, f$seconds, sum(abs(f$coef) > 1e-6))
    f <- raw_edge_ridge_woodbury(Xe[tr, ], y_tr, Xe[te, ]);      add("raw_edge_ridge", f$pred, f$seconds)
    f <- raw_edge_lasso(Xe[tr, ], y_tr, Xe[te, ]);               if (!is.null(f)) add("raw_edge_lasso", f$pred, f$seconds, f$nsel)
    f <- raw_edge_ard(Xe[tr, ], y_tr, Xe[te, ]);                 add("raw_edge_ard", f$pred, f$seconds, sum(abs(f$coef) > 1e-6))
    f <- spike_slab_ssvs(Xe[tr, ], y_tr, Xe[te, ]);              add("raw_edge_ssvs", f$pred, f$seconds, sum(f$pip > 0.5))
  } else {
    f <- graph_spectral_ard_logit(Xs[tr, ], y_tr, Xs[te, ]);     add("graph_spectral_ard", f$prob, f$seconds, f$nsel)
    f <- raw_edge_glmnet_logit(Xe[tr, ], y_tr, Xe[te, ], alpha = 0); add("raw_edge_ridge", f$prob, f$seconds, f$nsel)
    f <- raw_edge_glmnet_logit(Xe[tr, ], y_tr, Xe[te, ], alpha = 1); add("raw_edge_lasso", f$prob, f$seconds, f$nsel)
    ## raw_edge_ard / spike_slab_ssvs are Gaussian regression fits (see file header) applied
    ## to y coded {0,1} -- pred/coef come out on a roughly [0,1]-centered continuous scale
    ## (standardize_train_test() adds back ybar = mean(y) in {0,1}), used here ONLY for AUC
    ## ranking, not as calibrated probabilities.
    f <- raw_edge_ard(Xe[tr, ], y_tr, Xe[te, ]);                 add("raw_edge_ard", f$pred, f$seconds, sum(abs(f$coef) > 1e-6))
    f <- spike_slab_ssvs(Xe[tr, ], y_tr, Xe[te, ]);              add("raw_edge_ssvs", f$pred, f$seconds, sum(f$pip > 0.5))
  }
  df <- do.call(rbind, rows)
  rownames(df) <- NULL
  df
}

run_five_method <- function(m, name, n_splits = 10, seed0 = 900, test_frac = 0.30,
                            K = 10, basis = "adjacency", out_csv,
                            max_seconds_per_fit = 600, n_splits_fallback = 3) {
  oc <- build_outcome_5m(m, name)
  conn <- oc$conn; meta <- oc$meta; y <- oc$y; type <- oc$type
  subj <- meta$SUBJECT; usub <- unique(subj)
  cat(sprintf("[%s] type=%s N=%d rows, %d subjects\n", name, type, length(conn), length(usub)))

  done <- integer(0)
  if (file.exists(out_csv)) { prev <- read.csv(out_csv); done <- unique(prev$split[prev$outcome == name]) }

  ## ---- feasibility check on split 1 (or first not-yet-done split) -------------------
  first_split <- setdiff(seq_len(n_splits), done)[1]
  n_splits_eff <- n_splits
  if (!is.na(first_split)) {
    set.seed(seed0 + first_split)
    te_s <- sample(usub, max(1, round(test_frac * length(usub))))
    te <- which(subj %in% te_s); tr <- which(!(subj %in% te_s))
    ft <- build_features(conn, tr, K, basis = basis)
    Xe <- ft$X_edge
    cat(sprintf("  feasibility check: p=%d, n_train=%d -- timing raw_edge_ard + spike_slab_ssvs on split %d...\n",
                ncol(Xe), length(tr), first_split))
    t_ard0 <- proc.time()[3]
    f_ard <- raw_edge_ard(Xe[tr, ], y[tr], Xe[te, ])
    t_ard <- proc.time()[3] - t_ard0
    t_ssvs0 <- proc.time()[3]
    f_ssvs <- spike_slab_ssvs(Xe[tr, ], y[tr], Xe[te, ])
    t_ssvs <- proc.time()[3] - t_ssvs0
    cat(sprintf("  TIMING: raw_edge_ard = %.1fs, spike_slab_ssvs = %.1fs (per fit, p=%d)\n",
                t_ard, t_ssvs, ncol(Xe)))
    if (max(t_ard, t_ssvs) > max_seconds_per_fit) {
      n_splits_eff <- min(n_splits, n_splits_fallback)
      cat(sprintf("  ** exceeds %ds/fit budget -- reducing to %d splits (was %d) per the feasibility guard **\n",
                  max_seconds_per_fit, n_splits_eff, n_splits))
    } else {
      cat(sprintf("  within %ds/fit budget -- proceeding with the full %d splits\n", max_seconds_per_fit, n_splits))
    }
  }

  for (s in seq_len(n_splits_eff)) {
    if (s %in% done) { cat(sprintf("[%s] split %d already done -- skipping\n", name, s)); next }
    t0 <- proc.time()[3]
    set.seed(seed0 + s)
    te_s <- sample(usub, max(1, round(test_frac * length(usub))))
    te <- which(subj %in% te_s); tr <- which(!(subj %in% te_s))
    if (type != "regression" && (length(unique(y[tr])) < 2 || length(unique(y[te])) < 2)) {
      cat(sprintf("[%s] split %d: degenerate label split -- skipping\n", name, s)); next
    }
    r <- eval_split_5m(conn, y, type, tr, te, K = K, basis = basis)
    r$outcome <- name; r$split <- s; r$type <- type; r$n <- length(conn); r$n_subj <- length(usub)
    prev <- if (file.exists(out_csv)) read.csv(out_csv) else NULL
    dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
    write.csv(rbind(prev, r), out_csv, row.names = FALSE)
    cat(sprintf("[%s] split %d done in %.1fs total -- %s\n", name, s, proc.time()[3] - t0,
                paste(sprintf("%s=%.3f", r$method, if (type == "regression") r$rmse else r$auc), collapse = " ")))
    flush.console()
  }
  invisible(TRUE)
}
