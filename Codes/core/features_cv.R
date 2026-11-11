## =====================================================================================
## Speed refactor + basis/K selection by cross-validation (fairness fix, mirror of
## comparators_cv.R but for the graph-spectral side of the model).
##
## Background. build_features() (Codes/core/abide_pipeline.R) recomputes the RAW-EDGE design
## matrix (N x p) every time it is called, even though that matrix depends only on `conn` --
## not on train_idx, K, or basis. A basis/K selection loop that tries a dozen (basis, K)
## combinations, each of which currently calls build_features() and therefore rebuilds the
## edge matrix from scratch every time, would waste most of its time on a matrix that never
## changes. This file splits build_features() into its two independent halves:
##
##   build_edge_features(conn)                    -- the N x p raw-edge design; compute ONCE.
##   build_spec_features(conn, train_idx, K, basis) -- the N x q spectral design for one
##                                                    (train_idx, K, basis) combination.
##
## build_features() ITSELF IS NOT TOUCHED and keeps returning exactly what it always has, for
## every script that already depends on it.
##
## select_basis_K() then chooses (basis, K) by cross-validation, using ONLY the training rows,
## mirroring the same untuned-vs-tuned fairness fix as comparators_cv.R for ridge/lasso.
##
## ---------------------------------------------------------------------------------------
## WHY IT IS OKAY (AND CHEAPER) TO BUILD THE SPECTRAL DESIGN ONCE PER (basis, K), RATHER THAN
## RE-DERIVING IT INSIDE EVERY INNER FOLD:
##
## The graph basis Q (and hence the spectral design X_spec = vech(Q' A Q)) is built from
## make_basis() applied to the MEAN of the TRAINING connectomes only (build_spec_features()
## below does exactly what build_features() already does: Wbar <- mean(|conn[train_idx]|),
## then eigen-decompose Wbar). That construction:
##   * never touches a TEST-set connectome (only conn[train_idx] contributes to Wbar), and
##   * never touches the outcome y at all -- it is a purely unsupervised summary of the
##     training design's own graph geometry.
## Because the basis construction is unsupervised and test-blind, re-deriving it separately
## inside each inner CV fold buys NO additional protection against leakage -- the inner test
## fold's connectomes still never enter Wbar either way, since Wbar is built from train_idx
## (the OUTER training rows), a superset of every inner fold's inner-training rows. Rebuilding
## it 5 (or however many) extra times per (basis, K) combination would only spend compute, not
## reduce risk. So select_basis_K() builds the spectral design ONCE per (basis, K) combination
## from the full outer training set, and cross-validates only the SUPERVISED model fit
## (graph_spectral_ard / graph_spectral_ard_logit) across inner folds on top of that fixed
## design. No test-set connectome and no outcome value ever enters basis construction.
## =====================================================================================

suppressMessages(source("Codes/core/abide_classification.R"))  # -> abide_pipeline.R -> network_bayes_simulation.R
                                                                 # brings in: make_basis, vech_features,
                                                                 # upper_edge_index, graph_spectral_ard,
                                                                 # graph_spectral_ard_logit, auc_pred, rmse

## ---- 2a. split build_features() into its two independent halves --------------------

## The N x p raw-edge design. Does not depend on train_idx, K, or basis -- compute ONCE per
## dataset and reuse across every (basis, K) combination and every split.
build_edge_features <- function(conn) {
  V <- nrow(conn[[1]])
  edge_idx <- upper_edge_index(V)
  N <- length(conn)
  X_edge <- matrix(0, N, nrow(edge_idx))
  for (i in seq_len(N)) X_edge[i, ] <- conn[[i]][edge_idx]
  list(X_edge = X_edge, p = nrow(edge_idx), edge_idx = edge_idx)
}

## The N x q spectral design for one (train_idx, K, basis) combination. Basis is built from
## the TRAINING connectomes' mean graph only (see header); returns the design for ALL N rows
## (so it can be sliced by both train_idx and test_idx), plus the basis Q and its eigenvalues.
build_spec_features <- function(conn, train_idx, K, basis = "adjacency") {
  V <- nrow(conn[[1]])
  Wbar <- Reduce(`+`, lapply(conn[train_idx], abs)) / length(train_idx)
  diag(Wbar) <- 0
  bs <- make_basis(Wbar, K, basis)
  Q <- bs$Q

  N <- length(conn)
  q <- K * (K + 1) / 2
  X_spec <- matrix(0, N, q)
  for (i in seq_len(N)) {
    C <- t(Q) %*% conn[[i]] %*% Q
    X_spec[i, ] <- vech_features(C)
  }
  list(X_spec = X_spec, Q = Q, evals = bs$evals, q = q)
}

## ---- 2b. choose (basis, K) by cross-validation on the TRAINING rows only -------------
##
## Scoring: type="regression" -> cross-validated MSE from graph_spectral_ard (lower better);
##          type="binary"     -> cross-validated AUC from graph_spectral_ard_logit (higher
##                                better).
## Inner folds are assigned at the SUBJECT level via `groups` whenever a subject vector is
## supplied (exactly like the outer train/test splits elsewhere in this project), so no
## subject straddles an inner fold boundary either. Any (basis, K) combination with
## q = K(K+1)/2 >= the number of rows in the SMALLEST inner-training fold is skipped (not
## identified there) and recorded as skipped in the returned grid.
##
## `X_edge` is accepted (default NULL) for interface symmetry with the rest of the pipeline
## and reserved for a future edge-aware selection criterion; the scoring rule specified for
## this fairness fix uses only the pure graph-spectral fits, never the hybrid, so X_edge is
## currently unused inside this function.
select_basis_K <- function(conn, train_idx, y, type, X_edge = NULL,
                           bases = c("adjacency", "nlaplacian", "laplacian"),
                           Ks = c(5, 10, 15, 20), inner_folds = 5, seed = 1,
                           groups = NULL) {
  stopifnot(type %in% c("regression", "binary"))
  t0 <- proc.time()[3]
  n_tr <- length(train_idx)

  ## ---- inner fold assignment over the TRAINING rows, at the SUBJECT level if possible ----
  set.seed(seed)
  if (!is.null(groups)) {
    g_tr <- groups[train_idx]
    usub <- unique(g_tr)
    fold_of_subj <- sample(rep(seq_len(inner_folds), length.out = length(usub)))
    fold_id <- fold_of_subj[match(g_tr, usub)]
  } else {
    fold_id <- sample(rep(seq_len(inner_folds), length.out = n_tr))
  }
  inner_train_size <- sapply(seq_len(inner_folds), function(k) sum(fold_id != k))
  min_inner_train <- min(inner_train_size)

  grid <- expand.grid(basis = bases, K = Ks, stringsAsFactors = FALSE)
  grid$q <- grid$K * (grid$K + 1) / 2
  grid$skipped  <- grid$q >= min_inner_train
  grid$cv_score <- NA_real_

  for (g in seq_len(nrow(grid))) {
    if (grid$skipped[g]) next
    basis <- grid$basis[g]; K <- grid$K[g]

    ## Build the spectral design ONCE for this (basis, K), from the full OUTER training set
    ## (unsupervised; see header for why this is legitimate and does not need to be redone
    ## inside every inner fold).
    sf <- build_spec_features(conn, train_idx, K, basis = basis)
    Xs_full <- sf$X_spec

    fold_scores <- rep(NA_real_, inner_folds)
    for (k in seq_len(inner_folds)) {
      in_idx  <- train_idx[fold_id != k]
      out_idx <- train_idx[fold_id == k]
      if (length(out_idx) == 0 || length(in_idx) == 0) next
      Xs_in  <- Xs_full[in_idx, , drop = FALSE]
      Xs_out <- Xs_full[out_idx, , drop = FALSE]
      y_in <- y[in_idx]; y_out <- y[out_idx]

      if (type == "regression") {
        f <- graph_spectral_ard(Xs_in, y_in, Xs_out)
        fold_scores[k] <- mean((y_out - f$pred)^2)              # MSE: lower is better
      } else {
        if (length(unique(y_in)) < 2 || length(unique(y_out)) < 2) next
        f <- graph_spectral_ard_logit(Xs_in, y_in, Xs_out)
        fold_scores[k] <- auc_pred(f$prob, y_out)                # AUC: higher is better
      }
    }
    grid$cv_score[g] <- mean(fold_scores, na.rm = TRUE)
  }

  scored <- which(!is.na(grid$cv_score))
  if (length(scored) == 0) stop("select_basis_K: every (basis, K) combination was skipped -- ",
                                "reduce Ks or inner_folds; smallest inner-training fold has ",
                                min_inner_train, " rows")

  best_g <- if (type == "regression") scored[which.min(grid$cv_score[scored])]
            else scored[which.max(grid$cv_score[scored])]

  ## Rebuild the winning (basis, K) design once more (cheap relative to the whole grid) so the
  ## returned design matches grid$basis[best_g]/grid$K[best_g] exactly.
  sf_best <- build_spec_features(conn, train_idx, grid$K[best_g], basis = grid$basis[best_g])

  list(basis = grid$basis[best_g], K = grid$K[best_g], q = grid$q[best_g],
       grid = grid, cv_score_sel = grid$cv_score[best_g],
       X_spec = sf_best$X_spec, Q = sf_best$Q, evals = sf_best$evals,
       inner_folds = inner_folds, min_inner_train = min_inner_train,
       seconds = proc.time()[3] - t0)
}
