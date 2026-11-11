## =====================================================================================
## Hybrid graph-spectral (low-rank) + sparse-edge model.
##
## Motivation: across simulations and real data, graph-spectral ARD wins when the signal is
## low-rank/smooth in a graph basis but loses to dense/sparse raw-edge methods when the signal
## is distributed or edge-localized. A model that carries BOTH parts should be robust across
## regimes:
##
##     B = U_K Theta U_K^T  +  S            (low-rank graph-spectral  +  sparse edges)
##
## Fit by backfitting: (1) graph-spectral ARD on the spectral features for the low-rank part,
## (2) lasso (glmnet) on the raw-edge RESIDUAL for the sparse part; iterate. The lasso scales to
## large p, so the hybrid is feasible even at p = V(V-1)/2 ~ 20k where raw-edge ARD/SSVS are not.
##
## Reuses ard_fit / standardize_train_test from network_bayes_simulation.R; needs glmnet.
## =====================================================================================

graph_spectral_plus_sparse <- function(X_spec, X_edge, y, Xt_spec, Xt_edge,
                                        n_iter = 2L, lasso_nfolds = 5L) {
  stopifnot(requireNamespace("glmnet", quietly = TRUE))
  t0 <- proc.time()[3]
  ps <- standardize_train_test(X_spec, Xt_spec, y)     # spectral block (centers y -> ybar)
  pe <- standardize_train_test(X_edge, Xt_edge, y)     # raw-edge block (same ybar)
  yc <- ps$y

  theta <- rep(0, ncol(ps$X))
  beta  <- rep(0, ncol(pe$X))
  for (it in seq_len(n_iter)) {
    ## (1) low-rank part: ARD on spectral features, against the edge-adjusted residual
    r_spec <- yc - as.numeric(pe$X %*% beta)
    theta  <- ard_fit(ps$X, r_spec)$coef
    ## (2) sparse part: lasso on raw edges, against the spectral-adjusted residual
    r_edge <- yc - as.numeric(ps$X %*% theta)
    cv <- glmnet::cv.glmnet(pe$X, r_edge, alpha = 1, nfolds = lasso_nfolds, standardize = FALSE)
    beta <- as.numeric(coef(cv, s = "lambda.min"))[-1]
  }
  pred <- as.numeric(ps$Xt %*% theta + pe$Xt %*% beta + ps$ybar)
  list(coef_spec = theta, coef_edge = beta, pred = pred,
       seconds = proc.time()[3] - t0,
       n_spec = sum(abs(theta) > 1e-6), n_edges = sum(abs(beta) > 1e-10))
}

## Binary-outcome variant: same backfitting, logistic glmnet for BOTH parts (spectral ridge +
## edge lasso on the linear predictor). Returns predicted probabilities on the test set.
graph_spectral_plus_sparse_logit <- function(X_spec, X_edge, y, Xt_spec, Xt_edge,
                                              n_iter = 2L, nfolds = 5L) {
  stopifnot(requireNamespace("glmnet", quietly = TRUE))
  t0 <- proc.time()[3]
  ctr <- function(X) { m <- colMeans(X); s <- apply(X, 2, sd); s[s == 0 | is.na(s)] <- 1; list(m = m, s = s) }
  cs <- ctr(X_spec); ce <- ctr(X_edge)
  Xs <- scale(X_spec, cs$m, cs$s);  Xse <- scale(X_edge, ce$m, ce$s)
  Ts <- scale(Xt_spec, cs$m, cs$s); Tse <- scale(Xt_edge, ce$m, ce$s)

  off_tr <- rep(0, nrow(Xs)); off_te <- rep(0, nrow(Ts))
  for (it in seq_len(n_iter)) {
    ## spectral part (ridge logistic) with current edge offset
    cv1 <- glmnet::cv.glmnet(Xs, y, family = "binomial", alpha = 0, nfolds = nfolds,
                             standardize = FALSE, offset = off_tr)
    lp_s_tr <- as.numeric(predict(cv1, Xs,  s = "lambda.min", newoffset = off_tr, type = "link"))
    lp_s_te <- as.numeric(predict(cv1, Ts,  s = "lambda.min", newoffset = off_te, type = "link"))
    ## edge part (lasso logistic) with spectral linear predictor as offset
    cv2 <- glmnet::cv.glmnet(Xse, y, family = "binomial", alpha = 1, nfolds = nfolds,
                             standardize = FALSE, offset = lp_s_tr)
    off_tr <- as.numeric(predict(cv2, Xse, s = "lambda.min", newoffset = lp_s_tr, type = "link"))
    off_te <- as.numeric(predict(cv2, Tse, s = "lambda.min", newoffset = lp_s_te, type = "link"))
  }
  prob <- 1 / (1 + exp(-off_te))
  list(prob = prob, seconds = proc.time()[3] - t0)
}
