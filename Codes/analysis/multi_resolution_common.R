## =====================================================================================
## Shared helpers for the multi-resolution (Schaefer-100 / Schaefer-200) replication work.
##
## Does NOT modify Codes/analysis/outcome_explorer.R (owned jointly with parallel agents'
## work) -- that file's make_set()/run_outcome() are hardcoded to
## MASTER_RDS = "data/osf_master/master.rds" (the Power-229 master). This file sources
## outcome_explorer.R purely for its STATELESS helpers (eval_split, build_features,
## graph_spectral_ard[_logit], graph_spectral_plus_sparse[_logit], auc_pred, std_feats,
## raw_edge_ridge_woodbury, raw_edge_glmnet_logit, rmse -- none of which touch the global
## MASTER_RDS/.master), plus comparators_cv.R and features_cv.R (also stateless), and adds
## its own master-parameterized replacements for make_set()/run_outcome() so the SAME
## outcome-construction and split logic can run against ANY master.rds path (Power-229,
## Schaefer-100, Schaefer-200) without touching the original global-state functions at all.
## =====================================================================================

suppressMessages({
  source("Codes/analysis/outcome_explorer.R")
  source("Codes/core/comparators_cv.R")
  source("Codes/core/features_cv.R")
})

## ---- master-parameterized replacement for make_set() --------------------------------
## Identical logic to outcome_explorer.R::make_set(), but takes the already-loaded master
## list `m` (from readRDS(<atlas>/master.rds)) as an explicit argument instead of reading the
## global .master/MASTER_RDS.
make_set_m <- function(m, conds = NULL, aggregate_by_subject = FALSE) {
  s <- m$subjects
  keep <- if (is.null(conds)) seq_len(nrow(s)) else which(s$COND %in% conds)
  conn <- m$conn[keep]; meta <- s[keep, , drop = FALSE]
  if (aggregate_by_subject) {
    us <- unique(meta$SUBJECT)
    conn2 <- lapply(us, function(u) Reduce(`+`, conn[meta$SUBJECT == u]) / sum(meta$SUBJECT == u))
    meta2 <- meta[!duplicated(meta$SUBJECT), ]
    meta2 <- meta2[match(us, meta2$SUBJECT), , drop = FALSE]
    conn <- conn2; meta <- meta2
  }
  list(conn = conn, meta = meta)
}

## ---- master-parameterized replacement for run_outcome() ------------------------------
## Same resumable subject-level-split loop and CSV schema as outcome_explorer.R::run_outcome(),
## built on make_set_m() instead of the global make_set(). eval_split() itself (imported via
## the source() above) is unchanged and stateless, so results are directly comparable to the
## Power-229 runs whenever conds/y_col/type/n_splits/K/basis/seed0/methods match.
run_outcome_m <- function(m, name, conds = NULL, aggregate_by_subject = FALSE, y_col = NULL, type,
                          n_splits = 10, test_frac = 0.30, K = 10, basis = "adjacency",
                          methods = c("gs_ard", "hybrid", "ridge", "lasso"), seed0 = 900,
                          out_csv, verbose = TRUE, y_override_fn = NULL) {
  set <- make_set_m(m, conds, aggregate_by_subject)
  meta <- set$meta; conn <- set$conn
  y <- if (!is.null(y_override_fn)) y_override_fn(meta) else meta[[y_col]]
  ok <- which(!is.na(y) & y != "")
  conn <- conn[ok]; meta <- meta[ok, , drop = FALSE]; y <- y[ok]
  if (type != "multiclass") y <- as.numeric(y) else y <- as.character(y)
  subj <- meta$SUBJECT; usub <- unique(subj)
  if (verbose) cat(sprintf("\n[%s] type=%s  N=%d rows, %d subjects, conds=%s\n", name, type,
                           length(conn), length(usub), if (is.null(conds)) "all" else paste(conds, collapse = "+")))

  done <- integer(0)
  if (file.exists(out_csv)) { prev <- read.csv(out_csv)
    done <- unique(prev$split[prev$outcome == name]) }
  acc <- list()
  for (s in seq_len(n_splits)) {
    if (s %in% done) next
    set.seed(seed0 + s)
    te_s <- sample(usub, max(1, round(test_frac * length(usub))))    # SUBJECT-level split
    te <- which(subj %in% te_s); tr <- which(!(subj %in% te_s))
    if (type != "regression" && (length(unique(y[tr])) < 2 || length(unique(y[te])) < 2)) next
    r <- try(eval_split(conn, y, subj, type, tr, te, K = K, basis = basis, methods = methods), silent = TRUE)
    if (inherits(r, "try-error")) { cat("   split", s, "failed:", as.character(r), "\n"); next }
    r$outcome <- name; r$split <- s; r$type <- type; r$n <- length(conn); r$n_subj <- length(usub)
    acc[[length(acc) + 1]] <- r
    dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
    prev <- if (file.exists(out_csv)) read.csv(out_csv) else NULL
    write.csv(rbind(prev, do.call(rbind, acc[length(acc)])), out_csv, row.names = FALSE)
    if (verbose) {
      sig <- r[r$metric == "signal", ]
      cat(sprintf("   split %d: %s\n", s, paste(sprintf("%s=%.2f", sig$method, sig$value), collapse = " ")))
      flush.console()
    }
  }
  invisible(TRUE)
}

## ---- ridge-config helpers (config A = AUC-tuned 5-fold; config W = wide-path deviance) ----
## Mirrors ridge_collapse_diagnostic.R's ARM A / tuned_comparison.R's ridge_old (binary) and
## comparators_cv.R's raw_edge_glmnet_logit_wide (binary) / raw_edge_ridge_cv (regression),
## exposed here as small named wrappers so the nbk_0v1 multi-resolution table can report both
## configurations explicitly, per the task's "the win's configuration-dependence is a live
## issue" instruction.
ridge_configA_bin <- function(Xe_tr, y, Xe_te) raw_edge_glmnet_logit(Xe_tr, y, Xe_te, alpha = 0)
ridge_wide_bin     <- function(Xe_tr, y, Xe_te) raw_edge_glmnet_logit_wide(Xe_tr, y, Xe_te, alpha = 0,
                                                                          type_measure = "deviance",
                                                                          nfolds = 10, lambda_min_ratio = 1e-6)
lasso_configA_bin  <- function(Xe_tr, y, Xe_te) raw_edge_glmnet_logit(Xe_tr, y, Xe_te, alpha = 1)

## ---- Job 4a: nbk_0v1 replication at a given atlas, with BOTH ridge configs exposed ----
## Outcome-explorer's run_outcome_m() already handles gs_ard/hybrid/lasso via the standard
## "ridge"/"lasso" method names (config A everywhere), but Job 4a specifically needs BOTH ridge
## configs side by side (the Power-229 finding was that the nbk_1v2 win's significance flips
## between them) -- so this is a small bespoke loop rather than reusing run_outcome_m() as-is.
run_nbk_0v1_multires <- function(m, n_splits = 10, seed0 = 900, test_frac = 0.30, K = 10,
                                 basis = "adjacency", out_csv) {
  set <- make_set_m(m, conds = c("nbk_0back", "nbk_1back"))
  meta <- set$meta; conn <- set$conn
  y <- as.integer(meta$COND == "nbk_1back")
  subj <- meta$SUBJECT; usub <- unique(subj)
  cat(sprintf("[nbk_0v1] N=%d rows, %d subjects\n", length(conn), length(usub)))

  done <- integer(0)
  if (file.exists(out_csv)) { prev <- read.csv(out_csv); done <- unique(prev$split) }

  for (s in seq_len(n_splits)) {
    if (s %in% done) { cat(sprintf("  split %d already done -- skipping\n", s)); next }
    set.seed(seed0 + s)
    te_s <- sample(usub, max(1, round(test_frac * length(usub))))
    te <- which(subj %in% te_s); tr <- which(!(subj %in% te_s))
    if (length(unique(y[tr])) < 2 || length(unique(y[te])) < 2) { cat("  degenerate split", s, "-- skip\n"); next }

    ft <- build_features(conn, tr, K, basis = basis)
    Xs <- ft$X_spec; Xe <- ft$X_edge
    y_tr <- y[tr]; y_te <- y[te]

    gs <- graph_spectral_ard_logit(Xs[tr, ], y_tr, Xs[te, ])
    hy <- graph_spectral_plus_sparse_logit(Xs[tr, ], Xe[tr, ], y_tr, Xs[te, ], Xe[te, ])
    rgA <- ridge_configA_bin(Xe[tr, ], y_tr, Xe[te, ])
    rgW <- ridge_wide_bin(Xe[tr, ], y_tr, Xe[te, ])
    ls <- lasso_configA_bin(Xe[tr, ], y_tr, Xe[te, ])

    row <- data.frame(
      split = s,
      gs_ard_auc     = auc_pred(gs$prob, y_te),
      hybrid_auc     = auc_pred(hy$prob, y_te),
      ridge_configA_auc = auc_pred(rgA$prob, y_te),
      ridge_wide_auc    = auc_pred(rgW$prob, y_te),
      lasso_auc      = auc_pred(ls$prob, y_te),
      n_train = length(tr), n_test = length(te)
    )
    prev <- if (file.exists(out_csv)) read.csv(out_csv) else NULL
    dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
    write.csv(rbind(prev, row), out_csv, row.names = FALSE)
    cat(sprintf("  split %d: gs=%.3f hybrid=%.3f ridgeA=%.3f ridgeW=%.3f lasso=%.3f\n",
                s, row$gs_ard_auc, row$hybrid_auc, row$ridge_configA_auc, row$ridge_wide_auc, row$lasso_auc))
    flush.console()
  }
  invisible(TRUE)
}
