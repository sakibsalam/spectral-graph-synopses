## =====================================================================================
## Real-data pipeline: graph-spectral Bayesian shrinkage on ABIDE connectomes
##
## Reads the ABIDE PCP connectomes produced by Codes/fetch/abide_fetch.py and runs the
## scalar-on-network regression of the manuscript on a REAL phenotype (default: age).
##
## Follows the manuscript's real-data analysis plan (Sec. "Real data analysis plan"):
##   1. Build subject-level adjacency matrices (done in the Python fetch step).
##   2. Train/test split.
##   3. Build template graph \bar A from TRAINING data only.
##   4. Compute leading K Laplacian eigenvectors -> graph-spectral predictors.
##   5. Fit graph-spectral ARD and compare with scalable raw-edge baselines.
##   6. Report prediction accuracy and runtime.
##
## At V=200 the raw-edge dimension is p = V(V-1)/2 = 19,900. Raw-edge ARD and spike-slab
## SSVS are computationally prohibitive at that p (O(p^3) / dense p x p crossproducts) --
## which is precisely the motivation for the graph-spectral reduction. The scalable
## raw-edge comparators we run are ridge (via the Woodbury identity, an n x n solve) and
## lasso (glmnet), the latter standing in for edge-level selection.
## =====================================================================================

suppressMessages({
  source("Codes/core/network_bayes_simulation.R")   # ard_fit, standardize_train_test, vech_features,
                                                # laplacian_basis, upper_edge_index, rmse,
                                                # graph_spectral_ard, graph_spectral_projection
  has_glmnet <- requireNamespace("glmnet", quietly = TRUE)
})

DATA <- "data/abide"

## ---- load --------------------------------------------------------------------------

## Loads any dataset saved in the connectomes_full.csv / subjects.csv / meta.csv format
## (ABIDE by default; pass dir="data/acpi" etc. for other sources processed the same way).
load_abide <- function(dir = DATA) {
  meta <- read.csv(file.path(dir, "meta.csv"))
  V <- meta$V[1]; N <- meta$N[1]
  flat <- matrix(scan(file.path(dir, "connectomes_full.csv"), sep = ",", quiet = TRUE),
                 ncol = V, byrow = TRUE)              # (N*V) x V
  stopifnot(nrow(flat) == N * V)
  conn <- lapply(seq_len(N), function(i) flat[((i - 1) * V + 1):(i * V), , drop = FALSE])
  subjects <- read.csv(file.path(dir, "subjects.csv"))
  list(conn = conn, subjects = subjects, V = V, N = N, meta = meta)
}

## ---- feature construction ----------------------------------------------------------
## Template graph from TRAINING connectomes only; nonnegative weights via mean |Z| so the
## graph Laplacian is well defined. Spectral covariates use the signed connectomes.

## Build a K-dimensional orthonormal graph basis from a template adjacency Wbar.
##   "laplacian"  : leading (smoothest) eigenvectors of L = D - W       (graph smoothness)
##   "nlaplacian" : leading eigenvectors of the normalized Laplacian    (degree-robust)
##   "adjacency"  : dominant connectivity modes = eigenvectors of W with largest |eigenvalue|
make_basis <- function(Wbar, K, basis = "laplacian") {
  if (basis == "laplacian") {
    lb <- laplacian_basis(Wbar)
    list(Q = lb$vectors[, seq_len(K), drop = FALSE], evals = lb$values[seq_len(K)])
  } else if (basis == "nlaplacian") {
    d <- rowSums(Wbar); dinv <- ifelse(d > 0, 1 / sqrt(d), 0)
    Ln <- diag(length(d)) - (dinv * Wbar) %*% diag(dinv)
    Ln <- (Ln + t(Ln)) / 2
    eg <- eigen(Ln, symmetric = TRUE); ord <- order(eg$values)
    list(Q = eg$vectors[, ord[seq_len(K)], drop = FALSE], evals = eg$values[ord[seq_len(K)]])
  } else if (basis == "adjacency") {
    eg <- eigen(Wbar, symmetric = TRUE); ord <- order(abs(eg$values), decreasing = TRUE)
    list(Q = eg$vectors[, ord[seq_len(K)], drop = FALSE], evals = eg$values[ord[seq_len(K)]])
  } else stop("unknown basis: ", basis)
}

build_features <- function(conn, train_idx, K, basis = "laplacian") {
  V <- nrow(conn[[1]])
  Wbar <- Reduce(`+`, lapply(conn[train_idx], abs)) / length(train_idx)
  diag(Wbar) <- 0
  bs <- make_basis(Wbar, K, basis)
  Q <- bs$Q

  edge_idx <- upper_edge_index(V)
  N <- length(conn)
  q <- K * (K + 1) / 2
  X_spec <- matrix(0, N, q)
  X_edge <- matrix(0, N, nrow(edge_idx))
  for (i in seq_len(N)) {
    A <- conn[[i]]
    X_edge[i, ] <- A[edge_idx]
    C <- t(Q) %*% A %*% Q
    X_spec[i, ] <- vech_features(C)
  }
  list(X_spec = X_spec, X_edge = X_edge, Q = Q, evals = bs$evals, q = q,
       p = nrow(edge_idx))
}

## ---- scalable raw-edge comparators -------------------------------------------------

raw_edge_ridge_woodbury <- function(X, y, Xt, alpha = 50) {
  ## Ridge for n < p via the dual form: beta = X^T (X X^T + alpha I_n)^{-1} y.
  t0 <- proc.time()[3]
  pp <- standardize_train_test(X, Xt, y)
  n <- nrow(pp$X)
  G <- tcrossprod(pp$X) + alpha * diag(n)           # n x n
  w <- solve(G, pp$y)
  beta <- crossprod(pp$X, w)                        # p-vector
  pred <- as.numeric(pp$Xt %*% beta + pp$ybar)
  list(coef = as.numeric(beta), pred = pred, seconds = proc.time()[3] - t0)
}

raw_edge_lasso <- function(X, y, Xt) {
  if (!has_glmnet) return(NULL)
  t0 <- proc.time()[3]
  pp <- standardize_train_test(X, Xt, y)
  cv <- glmnet::cv.glmnet(pp$X, pp$y, alpha = 1, nfolds = 5, standardize = FALSE)
  b <- as.numeric(glmnet::coef.glmnet(cv, s = "lambda.min"))[-1]  # drop intercept
  pred <- as.numeric(pp$Xt %*% b + pp$ybar)
  list(coef = b, pred = pred, seconds = proc.time()[3] - t0,
       nsel = sum(abs(b) > 1e-10))
}

## ---- one train/test evaluation -----------------------------------------------------

run_abide <- function(outcome = "AGE_AT_SCAN",
                      bases = c("laplacian", "nlaplacian", "adjacency"),
                      Ks = c(10, 20, 30),
                      test_frac = 0.30, seed = 1, dir = DATA,
                      out_csv = "results/abide_results.csv") {
  d <- load_abide(dir)
  y_all <- d$subjects[[outcome]]
  keep <- which(is.finite(y_all))
  conn <- d$conn[keep]; y_all <- y_all[keep]
  N <- length(conn)
  cat(sprintf("ABIDE: N=%d usable subjects, V=%d nodes, outcome=%s (range %.1f-%.1f)\n",
              N, d$V, outcome, min(y_all), max(y_all)))

  set.seed(seed)
  test_idx <- sample(seq_len(N), round(test_frac * N))
  train_idx <- setdiff(seq_len(N), test_idx)
  y_tr <- y_all[train_idx]; y_te <- y_all[test_idx]
  baseline_rmse <- rmse(y_te, rep(mean(y_tr), length(y_te)))
  cat(sprintf("Split: %d train / %d test. Baseline (predict train mean) RMSE=%.3f\n",
              length(train_idx), length(test_idx), baseline_rmse))

  rows <- list()
  addrow <- function(method, dim, pred, secs, extra = NA_real_) {
    data.frame(outcome = outcome, method = method, dimension = dim,
               test_rmse = rmse(y_te, pred), seconds = round(secs, 3),
               extra = extra, stringsAsFactors = FALSE)
  }

  ## ---- graph-spectral ARD over (basis x K) ----
  for (basis in bases) {
    for (K in Ks) {
      ft <- build_features(conn, train_idx, K, basis = basis)
      Xs_tr <- ft$X_spec[train_idx, , drop = FALSE]; Xs_te <- ft$X_spec[test_idx, , drop = FALSE]
      f <- graph_spectral_ard(Xs_tr, y_tr, Xs_te)
      rows[[length(rows)+1]] <- addrow(sprintf("graph-spectral ARD [%s, K=%d]", basis, K),
                                       ft$q, f$pred, f$seconds, extra = sum(abs(f$coef) > 1e-6))
    }
  }
  ## graph-spectral projection: one cheap config (it is the slow method; K=10 Laplacian)
  ft <- build_features(conn, train_idx, 10, basis = "laplacian")
  Xs_tr <- ft$X_spec[train_idx, , drop = FALSE]; Xs_te <- ft$X_spec[test_idx, , drop = FALSE]
  f <- graph_spectral_projection(Xs_tr, y_tr, Xs_te)
  rows[[length(rows)+1]] <- addrow("graph-spectral projection [laplacian, K=10]", ft$q, f$pred,
                                   f$seconds, extra = sum(abs(f$coef) > 1e-10))

  ## ---- scalable raw-edge baselines (built once; edge features are K-independent) ----
  ft <- build_features(conn, train_idx, Ks[1])
  Xe_tr <- ft$X_edge[train_idx, , drop = FALSE]; Xe_te <- ft$X_edge[test_idx, , drop = FALSE]

  f <- raw_edge_ridge_woodbury(Xe_tr, y_tr, Xe_te)
  rows[[length(rows)+1]] <- addrow("raw-edge ridge (Woodbury)", ft$p, f$pred, f$seconds)

  f <- raw_edge_lasso(Xe_tr, y_tr, Xe_te)
  if (!is.null(f)) {
    rows[[length(rows)+1]] <- addrow("raw-edge lasso (glmnet)", ft$p, f$pred, f$seconds, extra = f$nsel)
  } else {
    cat("  (glmnet not available -- skipping raw-edge lasso)\n")
  }

  cat(sprintf("\nNOTE: raw-edge ARD and spike-slab SSVS are computationally prohibitive at p=%d\n", ft$p))
  cat("      (dense p x p crossproducts / O(p^3) solves) -- this is the motivation for the\n")
  cat("      graph-spectral reduction. Skipped here by design.\n")

  df <- do.call(rbind, rows)
  df$rmse_vs_baseline <- round(df$test_rmse / baseline_rmse, 3)
  dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
  write.csv(df, out_csv, row.names = FALSE)

  cat("\n================ ABIDE real-data results ================\n")
  print(df[order(df$test_rmse), ], row.names = FALSE)
  cat(sprintf("\n(baseline RMSE = %.3f; 'extra' = # selected/active components; lower test_rmse better)\n",
              baseline_rmse))
  cat("Saved to", out_csv, "\n")
  invisible(df)
}

## ---- repeated-split robustness for the headline comparison -------------------------
## A single 70/30 split is noisy at N=200. Average test RMSE over several random splits
## for the parsimonious graph-spectral model (adjacency, K=10) vs the scalable raw-edge
## baselines, so the conclusion does not hinge on one split.

run_abide_repeated <- function(outcome = "AGE_AT_SCAN", n_splits = 10, K = 10,
                               basis = "adjacency", test_frac = 0.30, dir = DATA,
                               out_csv = "results/abide_repeated.csv") {
  d <- load_abide(dir)
  y_all <- d$subjects[[outcome]]
  keep <- which(is.finite(y_all))
  conn <- d$conn[keep]; y_all <- y_all[keep]
  N <- length(conn)
  cat(sprintf("Repeated-split ABIDE (%s): N=%d, %d splits, graph basis=%s K=%d\n",
              outcome, N, n_splits, basis, K))

  acc <- list()
  for (s in seq_len(n_splits)) {
    set.seed(100 + s)
    test_idx <- sample(seq_len(N), round(test_frac * N))
    train_idx <- setdiff(seq_len(N), test_idx)
    y_tr <- y_all[train_idx]; y_te <- y_all[test_idx]
    base_rmse <- rmse(y_te, rep(mean(y_tr), length(y_te)))

    ft  <- build_features(conn, train_idx, K, basis = basis)
    Xs_tr <- ft$X_spec[train_idx, ]; Xs_te <- ft$X_spec[test_idx, ]
    Xe_tr <- ft$X_edge[train_idx, ]; Xe_te <- ft$X_edge[test_idx, ]

    gs  <- graph_spectral_ard(Xs_tr, y_tr, Xs_te)
    rg  <- raw_edge_ridge_woodbury(Xe_tr, y_tr, Xe_te)
    ls  <- raw_edge_lasso(Xe_tr, y_tr, Xe_te)
    one <- data.frame(
      split = s,
      gs_ard   = rmse(y_te, gs$pred),
      ridge    = rmse(y_te, rg$pred),
      lasso    = if (!is.null(ls)) rmse(y_te, ls$pred) else NA_real_,
      baseline = base_rmse)
    acc[[s]] <- one
  }
  df <- do.call(rbind, acc)
  dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
  write.csv(df, out_csv, row.names = FALSE)

  cat("\n========= Repeated-split mean test RMSE (lower better) =========\n")
  msd <- function(x) sprintf("%.3f +/- %.3f", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
  cat(sprintf("  graph-spectral ARD [%s, K=%d, dim=%d] : %s\n", basis, K, ft$q, msd(df$gs_ard)))
  cat(sprintf("  raw-edge ridge (Woodbury, dim=%d)     : %s\n", ft$p, msd(df$ridge)))
  cat(sprintf("  raw-edge lasso  (glmnet,   dim=%d)    : %s\n", ft$p, msd(df$lasso)))
  cat(sprintf("  baseline (predict train mean)         : %s\n", msd(df$baseline)))
  cat("Saved to", out_csv, "\n")
  invisible(df)
}

if (sys.nframe() == 0) {
  run_abide()
  cat("\n")
  run_abide_repeated()
}
