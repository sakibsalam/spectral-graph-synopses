## =====================================================================================
## Properly cross-validated raw-edge comparators (fairness fix).
##
## Background. Codes/core/network_bayes_simulation.R::raw_edge_ridge() and
## Codes/core/abide_pipeline.R::raw_edge_ridge_woodbury() use a HARDCODED penalty (alpha =
## 50, no tuning at all), while the classification comparators (raw_edge_glmnet_logit(), in
## Codes/core/abide_classification.R) at least run 5-fold CV -- but results/
## ridge_collapse_diagnostic.csv shows that CV, with glmnet's default (type.measure="auc",
## narrow default lambda path when n < p) frequently selects a lambda within 88-100% of the
## largest lambda on the path, i.e. an intercept-only / near-constant-probability model. That
## is not a real "raw-edge ridge/lasso can't find the signal" result -- it is a tuning
## artifact from a penalty grid that is too narrow.
##
## This file adds NEW, better-tuned comparators alongside the existing ones. It does NOT
## modify raw_edge_ridge_woodbury(), raw_edge_glmnet_logit(), or build_features() -- those
## keep working exactly as before for any script that still calls them.
##
## Contents:
##   raw_edge_ridge_cv()        -- CV-selected ridge for regression, dual/Gram form, so the
##                                 whole alpha grid costs about one fit per fold.
##   raw_edge_glmnet_logit_wide() -- same recipe as raw_edge_glmnet_logit() but with a
##                                 configurable, much WIDER lambda path (lambda.min.ratio),
##                                 metric, and fold count -- the classification analogue of
##                                 widening the ridge grid.
##   raw_edge_lasso_wide()      -- the same "widen the path" fix applied to the regression
##                                 lasso comparator (raw_edge_lasso() in abide_pipeline.R only
##                                 exercises glmnet's own default path), for parity with
##                                 raw_edge_ridge_cv() when Codes/analysis/tuned_comparison.R
##                                 contrasts "old" vs "cv" settings for BOTH penalty types.
## =====================================================================================

suppressMessages({
  source("Codes/core/abide_classification.R")   # -> abide_pipeline.R -> network_bayes_simulation.R
                                                  # brings in: standardize_train_test, std_feats,
                                                  # has_glmnet, rmse
})

## ---- 1a. CV-selected ridge for regression, dual (Gram) form -------------------------
##
## For n < p (true here: n ~ 70-270, p up to ~26,106) the primal p x p ridge solve is far
## more expensive than working in the dual: with standardized design X (n x p),
##   beta(alpha) = X' (X X' + alpha I_n)^{-1} y  =  X' V diag(1/(dd+alpha)) V' y
## where X X' = V diag(dd) V' is the (cheap, n x n) eigendecomposition of the Gram matrix.
## Because V and dd do not depend on alpha, once they are computed ONE time per CV fold, every
## candidate alpha costs only a diagonal rescale + two matrix-vector products -- so scanning
## n_alpha candidates costs about the same as ONE ridge fit, not n_alpha of them.
##
## The alpha grid is built from the Gram eigenvalues of the FULL standardized training design
## (computed once, outside the fold loop) so that it is scale-free -- i.e. centered at the
## natural scale of X'X for THIS design -- rather than an arbitrary absolute number like the
## old hardcoded 50. The same fixed grid (in absolute alpha units) is then evaluated inside
## each fold using that fold's own (much cheaper, n_fold x n_fold) Gram eigendecomposition.
raw_edge_ridge_cv <- function(X, y, Xt, nfolds = 5, n_alpha = 25, seed = 1) {
  t0 <- proc.time()[3]
  pp <- standardize_train_test(X, Xt, y)     # TRAIN-only centering/scaling; y centered
  Xs <- pp$X; yc <- pp$y
  n <- nrow(Xs)

  ## scale-free alpha grid, anchored on the FULL training design's Gram eigenvalues
  G_full  <- tcrossprod(Xs)                  # n x n
  eg_full <- eigen(G_full, symmetric = TRUE)
  dd_full <- eg_full$values
  alpha_grid <- mean(dd_full) * 10 ^ seq(-6, 6, length.out = n_alpha)

  set.seed(seed)
  fold_id <- sample(rep(seq_len(nfolds), length.out = n))

  cv_err <- matrix(NA_real_, nfolds, n_alpha)
  for (k in seq_len(nfolds)) {
    in_idx  <- which(fold_id != k)
    out_idx <- which(fold_id == k)
    Xf <- Xs[in_idx, , drop = FALSE];  yf <- yc[in_idx]
    Xh <- Xs[out_idx, , drop = FALSE]; yh <- yc[out_idx]

    G  <- tcrossprod(Xf)                     # n_fold x n_fold -- ONE eigendecomp per fold
    eg <- eigen(G, symmetric = TRUE)
    V <- eg$vectors; dd <- eg$values
    Vty <- as.numeric(crossprod(V, yf))

    for (a in seq_along(alpha_grid)) {
      alpha <- alpha_grid[a]
      w    <- V %*% (Vty / (dd + alpha))
      beta <- as.numeric(crossprod(Xf, w))
      pred <- as.numeric(Xh %*% beta)
      cv_err[k, a] <- mean((yh - pred)^2)
    }
  }

  mean_err  <- colMeans(cv_err, na.rm = TRUE)
  a_sel_idx <- which.min(mean_err)
  alpha_sel <- alpha_grid[a_sel_idx]

  ## refit on the FULL training set at the selected alpha (reuse the full-design eigendecomp)
  Vty_full <- as.numeric(crossprod(eg_full$vectors, yc))
  w_full   <- eg_full$vectors %*% (Vty_full / (dd_full + alpha_sel))
  beta     <- as.numeric(crossprod(Xs, w_full))
  pred     <- as.numeric(pp$Xt %*% beta + pp$ybar)

  list(pred = pred, coef = beta,
       alpha_sel = alpha_sel,
       alpha_grid_lo = min(alpha_grid), alpha_grid_hi = max(alpha_grid),
       cv_mse = mean_err[a_sel_idx],
       seconds = proc.time()[3] - t0)
}

## ---- 1b. glmnet logistic with a widened penalty path ---------------------------------
##
## Same fitting recipe as raw_edge_glmnet_logit() (family = "binomial", features standardized
## via std_feats(), predict at s = "lambda.min"), but exposing type.measure, nfolds AND
## lambda.min.ratio. glmnet's default lambda.min.ratio when n < p is 0.01 -- i.e. only 2 orders
## of magnitude of shrinkage are ever tried. Passing a much smaller ratio (default 1e-6 here)
## lets the path reach far less shrunk solutions, which is exactly where a real, non-collapsed
## fit can live in this n << p regime.
##
## lambda_frac_of_max (selected lambda / largest lambda on the path) and pred_sd (SD of the
## predicted test probabilities) are the collapse detectors: lambda_frac_of_max near 1 together
## with pred_sd near 0 means the "selected" model is, numerically, an intercept-only fit.
raw_edge_glmnet_logit_wide <- function(Xe_tr, y, Xe_te, alpha,
                                       type_measure = "deviance", nfolds = 10,
                                       lambda_min_ratio = 1e-6) {
  if (!has_glmnet) return(NULL)
  t0 <- proc.time()[3]
  pp <- std_feats(Xe_tr, Xe_te)
  cv <- glmnet::cv.glmnet(pp$X, y, family = "binomial", alpha = alpha, nfolds = nfolds,
                          standardize = FALSE, type.measure = type_measure,
                          lambda.min.ratio = lambda_min_ratio)
  prob <- as.numeric(predict(cv, newx = pp$Xt, s = "lambda.min", type = "response"))
  b    <- as.numeric(coef(cv, s = "lambda.min"))[-1]
  lam     <- cv$lambda.min
  lam_max <- max(cv$lambda)
  list(prob = prob,
       nsel = sum(abs(b) > 1e-10),
       lambda_sel = lam,
       lambda_frac_of_max = lam / lam_max,
       pred_sd = sd(prob),
       seconds = proc.time()[3] - t0)
}

## ---- 1c. regression-lasso analogue of the same "widen the path" fix -----------------
##
## Not explicitly requested as a named deliverable, but needed for Job 3's "lasso_old vs
## lasso_cv, the same contrast as ridge" on the REGRESSION outcomes (age_continuous):
## raw_edge_lasso() (Codes/core/abide_pipeline.R) only ever exercises glmnet's own default
## path (nfolds = 5, default lambda.min.ratio). This is the widened-path counterpart, mirroring
## raw_edge_glmnet_logit_wide() but for the Gaussian family, so "old" vs "cv" is an apples-to-
## apples contrast on both outcome types. Does not modify raw_edge_lasso() itself.
raw_edge_lasso_wide <- function(X, y, Xt, nfolds = 10, lambda_min_ratio = 1e-6) {
  if (!has_glmnet) return(NULL)
  t0 <- proc.time()[3]
  pp <- standardize_train_test(X, Xt, y)
  cv <- glmnet::cv.glmnet(pp$X, pp$y, alpha = 1, nfolds = nfolds, standardize = FALSE,
                          lambda.min.ratio = lambda_min_ratio)
  b    <- as.numeric(glmnet::coef.glmnet(cv, s = "lambda.min"))[-1]
  pred <- as.numeric(pp$Xt %*% b + pp$ybar)
  lam     <- cv$lambda.min
  lam_max <- max(cv$lambda)
  list(pred = pred, coef = b,
       nsel = sum(abs(b) > 1e-10),
       lambda_sel = lam,
       lambda_frac_of_max = lam / lam_max,
       pred_sd = sd(pred),
       seconds = proc.time()[3] - t0)
}
