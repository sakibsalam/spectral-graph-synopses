## =====================================================================================
## Fairness-tuned benchmark: old (hardcoded/narrow) vs. new (cross-validated) settings, on
## FOUR outcomes from the OSF Power-229 data, run through IDENTICAL paired splits.
##
## This does NOT tune anything toward a desired answer -- every arm below is tuned on
## TRAINING data only (see comparators_cv.R / features_cv.R headers for exactly how), and the
## script reports whatever numbers come out, including numbers that make our own method look
## worse than before.
##
## Arms (identical splits across all 8, so every comparison below is paired):
##   gs_ard_fixed  -- graph-spectral ARD at the hardcoded adjacency basis, K = 10
##   gs_ard_cv     -- graph-spectral ARD at the basis/K chosen by select_basis_K() on TRAIN
##   hybrid_fixed  -- hybrid (low-rank + sparse) at adjacency, K = 10
##   hybrid_cv     -- hybrid at the SAME basis/K select_basis_K() chose (not re-selected via
##                    the hybrid score -- select_basis_K always scores the pure graph-spectral
##                    fit, per its spec; hybrid_cv just reuses that choice)
##   ridge_old     -- regression: raw_edge_ridge_woodbury(alpha = 50, fixed).
##                    binary: raw_edge_glmnet_logit() (AUC metric, 5-fold, glmnet's own path).
##   ridge_cv      -- regression: raw_edge_ridge_cv() (dual-form CV, scale-free alpha grid).
##                    binary: raw_edge_glmnet_logit_wide() (deviance, 10-fold, widened path).
##   lasso_old / lasso_cv -- the same old-vs-new contrast, for the lasso penalty.
##
## Outcomes (subject-level 70/30 splits, set.seed(900 + s), exactly as run_outcome() in
## Codes/analysis/outcome_explorer.R):
##   nbk_0v1            -- nbk_0back vs nbk_1back, y = 1 for 1-back            (binary)
##   nbk_1v2            -- nbk_1back vs nbk_2back, y = 1 for 2-back            (binary)
##   age_continuous     -- all conditions, averaged per subject, outcome AGE  (regression)
##   age_binary_median  -- same aggregation, outcome AGE_BIN                  (binary)
##
## Outputs (new files; nothing under results/ is overwritten):
##   results/tuned_comparison.csv          -- one row per (outcome, split, arm)
##   results/tuned_comparison_summary.csv  -- one row per (outcome, arm), means across splits
##   results/tuned_comparison_paired_tests.csv -- paired Wilcoxon + t-test for the comparisons
##                                                 Job 3 asks for (cv arms vs ridge_cv, and the
##                                                 old-vs-new contrast for every arm)
##   results/tuned_comparison_permuted.csv -- Job 4's leakage control (age, outcome permuted
##                                             across subjects, 6 splits, full pipeline
##                                             including select_basis_K)
##
## Resumable: both the main loop and the permutation loop check what is already written to
## their target CSV and only compute missing (outcome, split) rows, so a crash loses at most
## the in-progress split.
## =====================================================================================

suppressMessages({
  source("Codes/analysis/outcome_explorer.R")   # make_set, build_features, graph_spectral_ard[_logit],
                                                  # graph_spectral_plus_sparse[_logit], auc_pred, std_feats,
                                                  # raw_edge_ridge_woodbury, raw_edge_glmnet_logit, rmse
  source("Codes/core/comparators_cv.R")          # raw_edge_ridge_cv, raw_edge_glmnet_logit_wide, raw_edge_lasso_wide
  source("Codes/core/features_cv.R")             # build_edge_features, build_spec_features, select_basis_K
})
stopifnot(requireNamespace("glmnet", quietly = TRUE))

## ---- outcome builders + expected-count guard ----------------------------------------

EXPECTED <- list(
  nbk_0v1           = list(n = 272, subj = 144),
  nbk_1v2           = list(n = 276, subj = 144),
  age_continuous    = list(n = 144, subj = 144),
  age_binary_median = list(n = 144, subj = 144)
)

build_outcome <- function(name) {
  if (name == "nbk_0v1") {
    set <- make_set(conds = c("nbk_0back", "nbk_1back"))
    y <- as.integer(set$meta$COND == "nbk_1back")
    type <- "binary"
  } else if (name == "nbk_1v2") {
    set <- make_set(conds = c("nbk_1back", "nbk_2back"))
    y <- as.integer(set$meta$COND == "nbk_2back")
    type <- "binary"
  } else if (name == "age_continuous") {
    set <- make_set(conds = NULL, aggregate_by_subject = TRUE)
    y <- as.numeric(set$meta$AGE)
    type <- "regression"
  } else if (name == "age_binary_median") {
    set <- make_set(conds = NULL, aggregate_by_subject = TRUE)
    y <- as.integer(set$meta$AGE_BIN)
    type <- "binary"
  } else stop("unknown outcome: ", name)
  list(conn = set$conn, meta = set$meta, y = y, type = type)
}

check_counts <- function(name, oc) {
  n <- length(oc$conn); nsub <- length(unique(oc$meta$SUBJECT))
  exp_n <- EXPECTED[[name]]$n; exp_subj <- EXPECTED[[name]]$subj
  cat(sprintf("[%s] N=%d rows, %d subjects (expected %d rows, %d subjects)\n",
              name, n, nsub, exp_n, exp_subj))
  if (n != exp_n || nsub != exp_subj) {
    stop(sprintf(
      "STOP: count mismatch for %s -- got %d rows / %d subjects, expected %d rows / %d subjects.",
      name, n, nsub, exp_n, exp_subj))
  }
  invisible(TRUE)
}

## ---- old-vs-cv raw-edge wrapper functions (diagnostics captured, nothing existing touched) --
## These duplicate the exact fitting recipes of the existing (unmodified) functions purely to
## expose reporting diagnostics (selected alpha/lambda, lambda_frac_of_max, nsel, pred_sd) that
## those functions do not currently return. They do not change raw_edge_ridge_woodbury(),
## raw_edge_glmnet_logit(), or raw_edge_lasso()'s behavior or signatures.

fit_ridge_old_reg <- function(Xe_tr, y_tr, Xe_te) {
  f <- raw_edge_ridge_woodbury(Xe_tr, y_tr, Xe_te, alpha = 50)
  list(pred = f$pred, alpha = 50, lambda_frac_of_max = NA_real_,
       n_nonzero = sum(abs(f$coef) > 1e-10), pred_sd = sd(f$pred), seconds = f$seconds)
}

fit_ridge_cv_reg <- function(Xe_tr, y_tr, Xe_te) {
  f <- raw_edge_ridge_cv(Xe_tr, y_tr, Xe_te)
  list(pred = f$pred, alpha = f$alpha_sel, lambda_frac_of_max = NA_real_,
       n_nonzero = sum(abs(f$coef) > 1e-10), pred_sd = sd(f$pred), seconds = f$seconds)
}

## Same recipe as raw_edge_lasso() in Codes/core/abide_pipeline.R (alpha=1, nfolds=5, glmnet's
## own default lambda path) -- reimplemented here only to expose lambda diagnostics.
fit_lasso_old_reg <- function(Xe_tr, y_tr, Xe_te) {
  t0 <- proc.time()[3]
  pp <- standardize_train_test(Xe_tr, Xe_te, y_tr)
  cv <- glmnet::cv.glmnet(pp$X, pp$y, alpha = 1, nfolds = 5, standardize = FALSE)
  b <- as.numeric(glmnet::coef.glmnet(cv, s = "lambda.min"))[-1]
  pred <- as.numeric(pp$Xt %*% b + pp$ybar)
  lam <- cv$lambda.min; lam_max <- max(cv$lambda)
  list(pred = pred, alpha = lam, lambda_frac_of_max = lam / lam_max,
       n_nonzero = sum(abs(b) > 1e-10), pred_sd = sd(pred), seconds = proc.time()[3] - t0)
}

fit_lasso_cv_reg <- function(Xe_tr, y_tr, Xe_te) {
  f <- raw_edge_lasso_wide(Xe_tr, y_tr, Xe_te, nfolds = 10, lambda_min_ratio = 1e-6)
  list(pred = f$pred, alpha = f$lambda_sel, lambda_frac_of_max = f$lambda_frac_of_max,
       n_nonzero = f$nsel, pred_sd = f$pred_sd, seconds = f$seconds)
}

## Same recipe as raw_edge_glmnet_logit() (AUC metric, 5-fold) -- alpha and lambda_frac_of_max
## are not exposed by that function and are recorded as NA here rather than modifying it.
fit_ridge_old_bin <- function(Xe_tr, y_tr, Xe_te) {
  f <- raw_edge_glmnet_logit(Xe_tr, y_tr, Xe_te, alpha = 0)
  list(prob = f$prob, alpha = NA_real_, lambda_frac_of_max = NA_real_,
       n_nonzero = f$nsel, pred_sd = sd(f$prob), seconds = f$seconds)
}

fit_lasso_old_bin <- function(Xe_tr, y_tr, Xe_te) {
  f <- raw_edge_glmnet_logit(Xe_tr, y_tr, Xe_te, alpha = 1)
  list(prob = f$prob, alpha = NA_real_, lambda_frac_of_max = NA_real_,
       n_nonzero = f$nsel, pred_sd = sd(f$prob), seconds = f$seconds)
}

fit_ridge_cv_bin <- function(Xe_tr, y_tr, Xe_te) {
  f <- raw_edge_glmnet_logit_wide(Xe_tr, y_tr, Xe_te, alpha = 0, type_measure = "deviance",
                                  nfolds = 10, lambda_min_ratio = 1e-6)
  list(prob = f$prob, alpha = f$lambda_sel, lambda_frac_of_max = f$lambda_frac_of_max,
       n_nonzero = f$nsel, pred_sd = f$pred_sd, seconds = f$seconds)
}

fit_lasso_cv_bin <- function(Xe_tr, y_tr, Xe_te) {
  f <- raw_edge_glmnet_logit_wide(Xe_tr, y_tr, Xe_te, alpha = 1, type_measure = "deviance",
                                  nfolds = 10, lambda_min_ratio = 1e-6)
  list(prob = f$prob, alpha = f$lambda_sel, lambda_frac_of_max = f$lambda_frac_of_max,
       n_nonzero = f$nsel, pred_sd = f$pred_sd, seconds = f$seconds)
}

## ---- one split, all 8 arms ------------------------------------------------------------

eval_split_tuned <- function(conn, y, subj, type, tr, te, Xe_all, seed_select = 1) {
  y_tr <- y[tr]; y_te <- y[te]
  Xe_tr <- Xe_all[tr, , drop = FALSE]; Xe_te <- Xe_all[te, , drop = FALSE]

  ## fixed spectral design: the old hardcoded adjacency / K = 10
  sf_fixed <- build_spec_features(conn, tr, K = 10, basis = "adjacency")
  Xsf <- sf_fixed$X_spec
  Xsf_tr <- Xsf[tr, , drop = FALSE]; Xsf_te <- Xsf[te, , drop = FALSE]

  ## cv-selected spectral design (shared by gs_ard_cv AND hybrid_cv, per spec)
  sel <- select_basis_K(conn, tr, y, type = type, groups = subj, seed = seed_select)
  Xsc <- sel$X_spec
  Xsc_tr <- Xsc[tr, , drop = FALSE]; Xsc_te <- Xsc[te, , drop = FALSE]

  rows <- list()
  addrow <- function(arm, pred = NULL, prob = NULL, basis = NA_character_, K = NA_real_,
                     alpha_or_lambda = NA_real_, lambda_frac_of_max = NA_real_,
                     n_nonzero = NA_real_, pred_sd = NA_real_, seconds = NA_real_) {
    if (type == "regression") {
      rmse_v <- rmse(y_te, pred); sp_v <- suppressWarnings(cor(y_te, pred, method = "spearman"))
      auc_v <- NA_real_; acc_v <- NA_real_
    } else {
      rmse_v <- NA_real_; sp_v <- NA_real_
      auc_v <- auc_pred(prob, y_te); acc_v <- mean((prob > 0.5) == y_te)
    }
    rows[[length(rows) + 1]] <<- data.frame(
      arm = arm, rmse = rmse_v, spearman = sp_v, auc = auc_v, acc = acc_v,
      basis = basis, K = K, alpha_or_lambda = alpha_or_lambda,
      lambda_frac_of_max = lambda_frac_of_max, n_nonzero = n_nonzero,
      pred_sd = pred_sd, seconds = seconds, stringsAsFactors = FALSE)
  }

  if (type == "regression") {
    f <- graph_spectral_ard(Xsf_tr, y_tr, Xsf_te)
    addrow("gs_ard_fixed", pred = f$pred, basis = "adjacency", K = 10,
           n_nonzero = sum(abs(f$coef) > 1e-6), pred_sd = sd(f$pred), seconds = f$seconds)

    f <- graph_spectral_ard(Xsc_tr, y_tr, Xsc_te)
    addrow("gs_ard_cv", pred = f$pred, basis = sel$basis, K = sel$K,
           n_nonzero = sum(abs(f$coef) > 1e-6), pred_sd = sd(f$pred), seconds = f$seconds)

    f <- graph_spectral_plus_sparse(Xsf_tr, Xe_tr, y_tr, Xsf_te, Xe_te)
    addrow("hybrid_fixed", pred = f$pred, basis = "adjacency", K = 10,
           n_nonzero = f$n_spec + f$n_edges, pred_sd = sd(f$pred), seconds = f$seconds)

    f <- graph_spectral_plus_sparse(Xsc_tr, Xe_tr, y_tr, Xsc_te, Xe_te)
    addrow("hybrid_cv", pred = f$pred, basis = sel$basis, K = sel$K,
           n_nonzero = f$n_spec + f$n_edges, pred_sd = sd(f$pred), seconds = f$seconds)

    f <- fit_ridge_old_reg(Xe_tr, y_tr, Xe_te)
    addrow("ridge_old", pred = f$pred, alpha_or_lambda = f$alpha,
           lambda_frac_of_max = f$lambda_frac_of_max, n_nonzero = f$n_nonzero,
           pred_sd = f$pred_sd, seconds = f$seconds)

    f <- fit_ridge_cv_reg(Xe_tr, y_tr, Xe_te)
    addrow("ridge_cv", pred = f$pred, alpha_or_lambda = f$alpha,
           lambda_frac_of_max = f$lambda_frac_of_max, n_nonzero = f$n_nonzero,
           pred_sd = f$pred_sd, seconds = f$seconds)

    f <- fit_lasso_old_reg(Xe_tr, y_tr, Xe_te)
    addrow("lasso_old", pred = f$pred, alpha_or_lambda = f$alpha,
           lambda_frac_of_max = f$lambda_frac_of_max, n_nonzero = f$n_nonzero,
           pred_sd = f$pred_sd, seconds = f$seconds)

    f <- fit_lasso_cv_reg(Xe_tr, y_tr, Xe_te)
    addrow("lasso_cv", pred = f$pred, alpha_or_lambda = f$alpha,
           lambda_frac_of_max = f$lambda_frac_of_max, n_nonzero = f$n_nonzero,
           pred_sd = f$pred_sd, seconds = f$seconds)

  } else {   # binary
    f <- graph_spectral_ard_logit(Xsf_tr, y_tr, Xsf_te)
    addrow("gs_ard_fixed", prob = f$prob, basis = "adjacency", K = 10,
           n_nonzero = f$nsel, pred_sd = sd(f$prob), seconds = f$seconds)

    f <- graph_spectral_ard_logit(Xsc_tr, y_tr, Xsc_te)
    addrow("gs_ard_cv", prob = f$prob, basis = sel$basis, K = sel$K,
           n_nonzero = f$nsel, pred_sd = sd(f$prob), seconds = f$seconds)

    f <- graph_spectral_plus_sparse_logit(Xsf_tr, Xe_tr, y_tr, Xsf_te, Xe_te)
    addrow("hybrid_fixed", prob = f$prob, basis = "adjacency", K = 10,
           pred_sd = sd(f$prob), seconds = f$seconds)   # nsel not exposed by this function

    f <- graph_spectral_plus_sparse_logit(Xsc_tr, Xe_tr, y_tr, Xsc_te, Xe_te)
    addrow("hybrid_cv", prob = f$prob, basis = sel$basis, K = sel$K,
           pred_sd = sd(f$prob), seconds = f$seconds)

    f <- fit_ridge_old_bin(Xe_tr, y_tr, Xe_te)
    addrow("ridge_old", prob = f$prob, alpha_or_lambda = f$alpha,
           lambda_frac_of_max = f$lambda_frac_of_max, n_nonzero = f$n_nonzero,
           pred_sd = f$pred_sd, seconds = f$seconds)

    f <- fit_ridge_cv_bin(Xe_tr, y_tr, Xe_te)
    addrow("ridge_cv", prob = f$prob, alpha_or_lambda = f$alpha,
           lambda_frac_of_max = f$lambda_frac_of_max, n_nonzero = f$n_nonzero,
           pred_sd = f$pred_sd, seconds = f$seconds)

    f <- fit_lasso_old_bin(Xe_tr, y_tr, Xe_te)
    addrow("lasso_old", prob = f$prob, alpha_or_lambda = f$alpha,
           lambda_frac_of_max = f$lambda_frac_of_max, n_nonzero = f$n_nonzero,
           pred_sd = f$pred_sd, seconds = f$seconds)

    f <- fit_lasso_cv_bin(Xe_tr, y_tr, Xe_te)
    addrow("lasso_cv", prob = f$prob, alpha_or_lambda = f$alpha,
           lambda_frac_of_max = f$lambda_frac_of_max, n_nonzero = f$n_nonzero,
           pred_sd = f$pred_sd, seconds = f$seconds)
  }

  df <- do.call(rbind, rows)
  rownames(df) <- NULL
  df
}

## ---- resumable outer loop over splits for one outcome --------------------------------

run_outcome_tuned <- function(name, n_splits = 10, test_frac = 0.30, seed0 = 900,
                              out_csv = "results/tuned_comparison.csv",
                              y_override = NULL, permuted = FALSE) {
  oc <- build_outcome(name)
  check_counts(name, oc)
  conn <- oc$conn; meta <- oc$meta; type <- oc$type
  y <- if (!is.null(y_override)) y_override else oc$y
  subj <- meta$SUBJECT; usub <- unique(subj)
  N <- length(conn)

  Xe_all <- build_edge_features(conn)$X_edge   # computed ONCE for this outcome (Job 2a)

  done <- integer(0)
  if (file.exists(out_csv)) {
    prev <- read.csv(out_csv)
    done <- unique(prev$split[prev$outcome == name & prev$permuted == permuted])
  }

  for (s in seq_len(n_splits)) {
    if (s %in% done) { cat(sprintf("[%s] split %d already done -- skipping\n", name, s)); next }
    t0 <- proc.time()[3]
    set.seed(seed0 + s)
    te_s <- sample(usub, max(1, round(test_frac * length(usub))))   # SUBJECT-level split
    te <- which(subj %in% te_s); tr <- which(!(subj %in% te_s))

    if (type != "regression" && (length(unique(y[tr])) < 2 || length(unique(y[te])) < 2)) {
      cat(sprintf("[%s] split %d: degenerate label split -- skipping\n", name, s)); next
    }

    r <- eval_split_tuned(conn, y, subj, type, tr, te, Xe_all, seed_select = seed0 + s)
    r$outcome <- name; r$split <- s; r$type <- type
    r$n <- N; r$n_subj <- length(usub); r$permuted <- permuted

    prev <- if (file.exists(out_csv)) read.csv(out_csv) else NULL
    write.csv(rbind(prev, r), out_csv, row.names = FALSE)

    cat(sprintf("[%s] split %d done in %.1fs -- %s\n", name, s, proc.time()[3] - t0,
                paste(sprintf("%s=%.3f", r$arm,
                              if (type == "regression") r$rmse else r$auc), collapse = " ")))
    flush.console()
  }
  invisible(TRUE)
}

## ---- summary + paired tests (Job 3 report) --------------------------------------------

build_summary <- function(df) {
  key <- interaction(df$outcome, df$arm, drop = TRUE)
  groups <- split(seq_len(nrow(df)), key)
  mode_chr <- function(x) { x <- x[!is.na(x)]; if (!length(x)) return(NA_character_)
                            names(sort(table(x), decreasing = TRUE))[1] }
  rows <- lapply(groups, function(idx) {
    d <- df[idx, ]
    data.frame(
      outcome = d$outcome[1], arm = d$arm[1], type = d$type[1], n_splits = nrow(d),
      mean_rmse = mean(d$rmse, na.rm = TRUE), mean_spearman = mean(d$spearman, na.rm = TRUE),
      mean_auc = mean(d$auc, na.rm = TRUE), mean_acc = mean(d$acc, na.rm = TRUE),
      basis_mode = mode_chr(d$basis),
      K_mode = suppressWarnings(as.numeric(mode_chr(as.character(d$K)))),
      mean_alpha_or_lambda = mean(d$alpha_or_lambda, na.rm = TRUE),
      mean_lambda_frac_of_max = mean(d$lambda_frac_of_max, na.rm = TRUE),
      mean_n_nonzero = mean(d$n_nonzero, na.rm = TRUE),
      mean_pred_sd = mean(d$pred_sd, na.rm = TRUE),
      mean_seconds = mean(d$seconds, na.rm = TRUE),
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  arm_order <- c("gs_ard_fixed", "gs_ard_cv", "hybrid_fixed", "hybrid_cv",
                 "ridge_old", "ridge_cv", "lasso_old", "lasso_cv")
  out$arm <- factor(out$arm, levels = arm_order)
  out[order(out$outcome, out$arm), ]
}

get_metric_vec <- function(df, arm_name, metric_col, splits) {
  sub <- df[df$arm == arm_name, ]
  sub[[metric_col]][match(splits, sub$split)]
}

paired_tests_for_outcome <- function(df, outcome_name, type) {
  sub <- df[df$outcome == outcome_name & df$permuted == FALSE, ]
  splits <- sort(unique(sub$split))
  metric_col <- if (type == "regression") "rmse" else "auc"

  pairs <- list(
    c("gs_ard_cv", "ridge_cv"),
    c("hybrid_cv", "ridge_cv"),
    c("gs_ard_fixed", "gs_ard_cv"),
    c("hybrid_fixed", "hybrid_cv"),
    c("ridge_old", "ridge_cv"),
    c("lasso_old", "lasso_cv")
  )
  out <- lapply(pairs, function(p) {
    a <- get_metric_vec(sub, p[1], metric_col, splits)
    b <- get_metric_vec(sub, p[2], metric_col, splits)
    ok <- is.finite(a) & is.finite(b)
    a <- a[ok]; b <- b[ok]
    wt <- suppressWarnings(try(wilcox.test(a, b, paired = TRUE), silent = TRUE))
    tt <- suppressWarnings(try(t.test(a, b, paired = TRUE), silent = TRUE))
    data.frame(
      outcome = outcome_name, metric = metric_col, arm_a = p[1], arm_b = p[2],
      mean_a = mean(a), mean_b = mean(b), mean_diff = mean(a - b), n = length(a),
      wilcox_p = if (inherits(wt, "try-error")) NA_real_ else wt$p.value,
      t_p = if (inherits(tt, "try-error")) NA_real_ else tt$p.value,
      stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

## ---- Job 4: permutation control (leakage check) ---------------------------------------

run_permutation_control <- function(n_splits = 6,
                                    out_csv = "results/tuned_comparison_permuted.csv") {
  oc <- build_outcome("age_continuous")
  check_counts("age_continuous", oc)
  set.seed(31337)                        # fixed, so the permutation itself is reproducible
  y_perm <- sample(oc$y)                  # permute AGE across subjects, once, before splitting
  run_outcome_tuned("age_continuous", n_splits = n_splits, seed0 = 900, out_csv = out_csv,
                    y_override = y_perm, permuted = TRUE)
  invisible(TRUE)
}

## FIXED (2026-09-10): a flat |signal| > 0.10 threshold is miscalibrated for this test-set
## size. run_permutation_control() splits at test_frac = 0.30 of n_subj = 144 subjects, i.e.
## a test set of ~43 subjects, where the null SD of Spearman's rho is ~= 1/sqrt(43 - 1) =
## 0.154 -- so |rho| > 0.10 fires on pure noise routinely (the previous version did exactly
## that: see the 2026-09-10 workbook entry). A leak can only ever inflate signal upward (never
## push it negative), so the correct test is a one-sample t-test of each arm's per-split
## Spearman values against zero, flagging only when the mean is POSITIVE and significantly so
## (p < ALPHA_LEAK, one-sided). We also report the theoretical null reference SD next to the
## observed spread so the two are directly comparable at a glance.
ALPHA_LEAK <- 0.05

check_permutation_leak <- function(perm_csv = "results/tuned_comparison_permuted.csv") {
  if (!file.exists(perm_csv)) stop("permutation control has not been run yet")
  d <- read.csv(perm_csv)
  d <- d[d$permuted == TRUE, ]

  ## theoretical null SD of Spearman's rho at this test-set size (test_frac * n_subj, same
  ## for every row of a given permutation run)
  n_test_ref <- round(0.30 * d$n_subj[1])
  sd_null_ref <- 1 / sqrt(max(n_test_ref - 1, 1))

  arms <- sort(unique(d$arm))
  rows <- lapply(arms, function(a) {
    x <- d$spearman[d$arm == a]
    x <- x[is.finite(x)]
    n <- length(x)
    if (n < 2) {
      return(data.frame(arm = a, n = n, mean_spearman = if (n == 1) x else NA_real_,
                        sd_spearman = NA_real_, t_p_one_sided = NA_real_,
                        sd_null_ref = sd_null_ref, leaked = FALSE, stringsAsFactors = FALSE))
    }
    tt <- t.test(x, alternative = "greater")   # H1: mean > 0 -- a leak inflates signal upward
    data.frame(arm = a, n = n, mean_spearman = mean(x), sd_spearman = sd(x),
              t_p_one_sided = tt$p.value, sd_null_ref = sd_null_ref,
              leaked = isTRUE(mean(x) > 0 && tt$p.value < ALPHA_LEAK),
              stringsAsFactors = FALSE)
  })
  agg <- do.call(rbind, rows)
  agg <- agg[order(-agg$mean_spearman), ]
  list(table = agg, leaked = any(agg$leaked, na.rm = TRUE), n_splits = length(unique(d$split)),
       n_test_ref = n_test_ref, sd_null_ref = sd_null_ref)
}

## ---- main ------------------------------------------------------------------------------

main <- function() {
  outcomes <- list(
    nbk_0v1           = "binary",
    nbk_1v2           = "binary",
    age_continuous     = "regression",
    age_binary_median  = "binary"
  )

  cat("\n================ JOB 3: tuned comparison, four outcomes ================\n")
  for (nm in names(outcomes)) {
    cat(sprintf("\n---- %s (%s) ----\n", nm, outcomes[[nm]]))
    run_outcome_tuned(nm, n_splits = 10, test_frac = 0.30, seed0 = 900,
                      out_csv = "results/tuned_comparison.csv")
  }

  main_df <- read.csv("results/tuned_comparison.csv")
  main_df <- main_df[main_df$permuted == FALSE, ]

  summ <- build_summary(main_df)
  write.csv(summ, "results/tuned_comparison_summary.csv", row.names = FALSE)
  cat("\nWrote results/tuned_comparison_summary.csv\n")
  print(summ, row.names = FALSE)

  pt <- do.call(rbind, lapply(names(outcomes), function(nm)
    paired_tests_for_outcome(main_df, nm, outcomes[[nm]])))
  write.csv(pt, "results/tuned_comparison_paired_tests.csv", row.names = FALSE)
  cat("\nWrote results/tuned_comparison_paired_tests.csv\n")
  print(pt, row.names = FALSE)

  cat("\n================ JOB 4: permutation (leakage) control ================\n")
  run_permutation_control(n_splits = 6, out_csv = "results/tuned_comparison_permuted.csv")

  leak <- check_permutation_leak()
  cat(sprintf("\n---- permuted age_continuous: one-sample t-test of Spearman vs 0, by arm (%d splits, n_test ~= %d subjects) ----\n",
              leak$n_splits, leak$n_test_ref))
  cat(sprintf("Theoretical null reference SD (1/sqrt(n_test-1)): %.4f -- compare to each arm's sd_spearman below.\n",
              leak$sd_null_ref))
  print(leak$table, row.names = FALSE)
  if (leak$leaked) {
    cat("\n*** LEAK DETECTED: at least one arm's mean permuted Spearman is POSITIVE and",
        sprintf("significantly greater than zero (one-sided p < %.2f). Do NOT treat", ALPHA_LEAK),
        "results/tuned_comparison_summary.csv as trustworthy until this is investigated. ***\n")
  } else {
    cat("\nNo leak detected: no arm's mean permuted Spearman is significantly positive at",
        sprintf("p < %.2f. A leak inflates signal upward, so negative or non-significant means", ALPHA_LEAK),
        "are consistent with clean data, even when their magnitude is comparable to the null SD above.\n")
  }

  invisible(list(summary = summ, paired_tests = pt, leak = leak))
}

if (sys.nframe() == 0) {
  main()
}
