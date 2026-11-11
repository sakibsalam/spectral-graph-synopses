## =====================================================================================
## Shared helpers for the ICDE experiments (E1-E4). Source from the repo root:
##   source("Codes/icde/icde_common.R")
## Nothing under Codes/core or Codes/analysis is modified; those files are only sourced.
## Outputs go to results/icde/ (created at runtime).
## =====================================================================================

suppressMessages({
  source("Codes/analysis/multi_resolution_common.R")  # -> outcome_explorer, comparators_cv, features_cv,
                                                       #    abide_classification, hybrid_model, core sim
                                                       # gives: make_set_m, make_basis, vech_features,
                                                       #  build_edge_features, auc_pred, std_feats, spearman,
                                                       #  graph_spectral_ard[_logit], ard_fit, ard_logistic_fit,
                                                       #  predict_logit, raw_edge_ridge_cv, raw_edge_glmnet_logit[_wide],
                                                       #  spike_slab_ssvs, standardize_train_test, rmse
  source("Codes/benchmarks/competitor_benchmarks.R")  # kundu_latent_nonlinear
})
stopifnot(requireNamespace("glmnet", quietly = TRUE))

## ---- command line ---------------------------------------------------------------------
has_flag <- function(f) f %in% commandArgs(trailingOnly = TRUE)
arg_val <- function(key, default = NA) {           # --key=value
  a <- grep(paste0("^", key, "="), commandArgs(trailingOnly = TRUE), value = TRUE)
  if (length(a)) sub(paste0("^", key, "="), "", a[1]) else default
}

OUT_DIR <- "results/icde"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
out_path <- function(name, smoke = FALSE) file.path(OUT_DIR, paste0(if (smoke) "smoke_" else "", name))

## ---- data: atlas x outcome -> list(conn, y, subj, type, V) ----------------------------
ATLAS_PATH <- c(power229 = "data/osf_master/master.rds",
                sch100   = "data/osf_schaefer100/master.rds")
.master_cache <- new.env()
load_set <- function(atlas, outcome) {
  if (is.null(.master_cache[[atlas]])) .master_cache[[atlas]] <- readRDS(ATLAS_PATH[[atlas]])
  m <- .master_cache[[atlas]]
  if (outcome == "nbk_0v1") {                      # same construction as tuned_comparison.R
    s <- make_set_m(m, conds = c("nbk_0back", "nbk_1back"))
    y <- as.integer(s$meta$COND == "nbk_1back"); type <- "binary"
  } else if (outcome == "age_continuous") {
    s <- make_set_m(m, conds = NULL, aggregate_by_subject = TRUE)
    y <- as.numeric(s$meta$AGE); type <- "regression"
  } else stop("unknown outcome ", outcome)
  list(conn = s$conn, y = y, subj = s$meta$SUBJECT, type = type, V = nrow(s$conn[[1]]))
}

## ---- splits: subject-level, seed 900 + split (identical to the existing scripts) -------
make_split <- function(subj, s, test_frac = 0.30, seed0 = 900) {
  usub <- unique(subj)
  set.seed(seed0 + s)
  te_s <- sample(usub, max(1, round(test_frac * length(usub))))
  list(tr = which(!(subj %in% te_s)), te = which(subj %in% te_s))
}
## Nested subject subsample of a training set: one fixed permutation per split, take the first
## ceiling(frac * n_subjects) subjects, so 0.25 subset of 0.5 subset of 0.75 subset of 1.0.
subsample_train <- function(tr, subj, frac, s, seed0 = 7000) {
  if (frac >= 1) return(tr)
  set.seed(seed0 + s)
  us <- sample(unique(subj[tr]))
  k <- max(4, ceiling(frac * length(us)))
  tr[subj[tr] %in% us[seq_len(k)]]
}

## ---- synopsis features (train-only template) ------------------------------------------
## Template = mean |A| over TRAINING graphs; basis = K leading adjacency eigenvectors.
template_Q <- function(conn, tr, K, basis = "adjacency") {
  Wbar <- Reduce(`+`, lapply(conn[tr], abs)) / length(tr); diag(Wbar) <- 0
  make_basis(Wbar, K, basis)$Q
}
feat_twosided <- function(conn, Q) {               # vech(Q' A Q), q = K(K+1)/2
  K <- ncol(Q); X <- matrix(0, length(conn), K * (K + 1) / 2)
  for (i in seq_along(conn)) X[i, ] <- vech_features(crossprod(Q, conn[[i]] %*% Q))
  X
}
feat_onesided <- function(conn, Q) {               # A Q flattened, V*R features (vertex identity kept)
  X <- matrix(0, length(conn), nrow(Q) * ncol(Q))
  for (i in seq_along(conn)) X[i, ] <- as.numeric(conn[[i]] %*% Q)
  X
}
feat_rownorm <- function(conn, Q) {                # one scalar per vertex: ||row of A Q||
  X <- matrix(0, length(conn), nrow(Q))
  for (i in seq_along(conn)) X[i, ] <- sqrt(rowSums((conn[[i]] %*% Q)^2))
  X
}

## ---- inner folds by subject, dual ridge with given folds, edge-space PCA -------------------
## Inner cross-validation folds formed over SUBJECTS: all graphs of a subject share a fold. With one
## graph per subject this is ordinary K-fold CV. raw_edge_ridge_cv (Codes/core) forms folds over
## graphs, which splits a subject's two graphs on the binary task and biases penalty selection
## towards heavy shrinkage; the ICDE experiments use the grouped version below on that task.
subject_folds <- function(subj_tr, k = 5, seed = 1) {
  us <- unique(subj_tr); set.seed(seed)
  sample(rep(seq_len(k), length.out = length(us)))[match(subj_tr, us)]
}
## Dual-Gram ridge with a CV penalty; same grid and refit as raw_edge_ridge_cv, folds passed in.
ridge_dual_cv_folds <- function(X, y, Xt, fold_id, n_alpha = 25) {
  pp <- standardize_train_test(X, Xt, y); Xs <- pp$X; yc <- pp$y
  eg0 <- eigen(tcrossprod(Xs), symmetric = TRUE)
  grid <- mean(eg0$values) * 10^seq(-6, 6, length.out = n_alpha)
  err <- matrix(NA_real_, max(fold_id), n_alpha)
  for (k in seq_len(max(fold_id))) {
    i <- which(fold_id != k); o <- which(fold_id == k)
    eg <- eigen(tcrossprod(Xs[i, , drop = FALSE]), symmetric = TRUE); Vty <- as.numeric(crossprod(eg$vectors, yc[i]))
    for (a in seq_along(grid)) {
      b <- as.numeric(crossprod(Xs[i, , drop = FALSE], eg$vectors %*% (Vty / (eg$values + grid[a]))))
      err[k, a] <- mean((yc[o] - as.numeric(Xs[o, , drop = FALSE] %*% b))^2)
    }
  }
  a <- which.min(colMeans(err))
  b <- as.numeric(crossprod(Xs, eg0$vectors %*% (as.numeric(crossprod(eg0$vectors, yc)) / (eg0$values + grid[a]))))
  list(pred = as.numeric(pp$Xt %*% b + pp$ybar), alpha_index = a, alpha_sel = grid[a], n_alpha = n_alpha)
}
## Edge-space principal component scores; loadings from the TRAINING graphs only (n x n Gram matrix).
pca_scores <- function(X, tr, q) {
  mu <- colMeans(X[tr, , drop = FALSE]); Xc <- sweep(X, 2, mu, "-")
  eg <- eigen(tcrossprod(Xc[tr, , drop = FALSE]), symmetric = TRUE); k <- min(q, length(tr) - 1)
  L <- crossprod(Xc[tr, , drop = FALSE], eg$vectors[, seq_len(k), drop = FALSE])
  Xc %*% sweep(L, 2, sqrt(pmax(eg$values[seq_len(k)], 1e-12)), "/")
}

## ---- effective-df fitters (copied from Codes/analysis/model_df_comparison.R, which runs
## its jobs at top level and therefore cannot be sourced). Same estimators, same constants.
gs_ard_logit_with_df <- function(Xs_tr, y_tr, Xs_te) {
  pp <- std_feats(Xs_tr, Xs_te)
  fit <- ard_logistic_fit(pp$X, y_tr)
  prob <- predict_logit(fit$intercept, fit$coef, pp$Xt)
  ## MacKay sum(gamma) recomputed at the reported mode (Hessian of the penalized logistic fit)
  Xd <- cbind(1, pp$X)
  eta <- as.numeric(Xd %*% c(fit$intercept, fit$coef))
  pr <- pmin(pmax(1 / (1 + exp(-eta)), 1e-6), 1 - 1e-6); W <- pr * (1 - pr)
  prec_full <- c(0, fit$prec)
  H <- crossprod(Xd * sqrt(W)) + diag(prec_full, length(prec_full))
  Sig <- tryCatch(solve(H), error = function(e) qr.solve(H))
  list(prob = prob, nominal = ncol(Xs_tr), eff_df = sum(pmax(0, 1 - prec_full * diag(Sig))),
       penalty = NA_real_)
}

## GLM ridge: eff df = 1 + trace(Xw (Xw'Xw + c I)^-1 Xw') = 1 + sum d/(d+c) over eigenvalues d of
## the n x n weighted dual Gram matrix at IRLS convergence; c = n * lambda_CV (glmnet, deviance).
ridge_logit_with_df <- function(Xe_tr, y_tr, Xe_te, nfolds = 5, seed = 1, max_iter = 30) {
  pp <- std_feats(Xe_tr, Xe_te); Xs <- pp$X; Xt <- pp$Xt
  n <- nrow(Xs); p <- ncol(Xs)
  set.seed(seed)
  cv <- glmnet::cv.glmnet(Xs, y_tr, family = "binomial", alpha = 0, nfolds = nfolds,
                          standardize = FALSE, type.measure = "deviance")
  lambda_sel <- cv$lambda.min; c_pen <- n * lambda_sel
  beta0 <- log((mean(y_tr) + 1e-9) / (1 - mean(y_tr) + 1e-9)); beta <- rep(0, p)
  for (it in seq_len(max_iter)) {
    eta <- beta0 + as.numeric(Xs %*% beta)
    pr <- pmin(pmax(1 / (1 + exp(-eta)), 1e-6), 1 - 1e-6); w <- pr * (1 - pr)
    z <- eta + (y_tr - pr) / w
    wbar <- sum(w); xbar <- as.numeric(crossprod(w, Xs)) / wbar; zbar <- sum(w * z) / wbar
    Xw <- sweep(Xs, 2, xbar, "-") * sqrt(w); zw <- (z - zbar) * sqrt(w)
    G <- tcrossprod(Xw)
    sol <- tryCatch(solve(G + c_pen * diag(n), zw), error = function(e) qr.solve(G + c_pen * diag(n), zw))
    beta_new <- as.numeric(crossprod(Xw, sol)); beta0_new <- zbar - sum(xbar * beta_new)
    delta <- sqrt(sum((beta_new - beta)^2) + (beta0_new - beta0)^2)
    beta <- beta_new; beta0 <- beta0_new
    if (delta < 1e-6 * max(1, sqrt(sum(beta^2)))) break
  }
  eta <- beta0 + as.numeric(Xs %*% beta)
  pr <- pmin(pmax(1 / (1 + exp(-eta)), 1e-6), 1 - 1e-6); w <- pr * (1 - pr)
  wbar <- sum(w); xbar <- as.numeric(crossprod(w, Xs)) / wbar
  Xw <- sweep(Xs, 2, xbar, "-") * sqrt(w)
  d <- pmax(eigen(tcrossprod(Xw), symmetric = TRUE, only.values = TRUE)$values, 0)
  prob <- 1 / (1 + exp(-(beta0 + as.numeric(Xt %*% beta))))
  list(prob = prob, nominal = p, eff_df = 1 + sum(d / (d + c_pen)), penalty = c_pen)
}

## Kundu re-implementation (rank R adjacency embedding, RBF kernel ridge); df = trace K(K+lambda I)^-1.
kundu_logit_with_df <- function(conn_tr, y_tr, conn_te, R = 8, ridge = 1e-2) {
  Vv <- nrow(conn_tr[[1]])
  Wbar <- Reduce(`+`, lapply(conn_tr, abs)) / length(conn_tr); diag(Wbar) <- 0
  eg <- eigen(Wbar, symmetric = TRUE); ord <- order(abs(eg$values), decreasing = TRUE)
  U <- eg$vectors[, ord[seq_len(R)], drop = FALSE]
  feat <- function(cl) t(sapply(cl, function(A) as.numeric(A %*% U)))
  Ftr <- feat(conn_tr); Fte <- feat(conn_te)
  cm <- colMeans(Ftr); cs <- apply(Ftr, 2, sd); cs[cs == 0] <- 1
  Ftr <- sweep(sweep(Ftr, 2, cm, "-"), 2, cs, "/"); Fte <- sweep(sweep(Fte, 2, cm, "-"), 2, cs, "/")
  sq <- function(A, B) outer(rowSums(A^2), rowSums(B^2), "+") - 2 * A %*% t(B)
  Dtr <- sq(Ftr, Ftr); sigma2 <- median(Dtr[upper.tri(Dtr)]); sigma2 <- ifelse(sigma2 > 0, sigma2, 1)
  Ktr <- exp(-Dtr / (2 * sigma2)); Kte <- exp(-sq(Fte, Ftr) / (2 * sigma2))
  a <- solve(Ktr + ridge * diag(nrow(Ktr)), y_tr - mean(y_tr))
  prob <- pmin(pmax(as.numeric(Kte %*% a) + mean(y_tr), 1e-6), 1 - 1e-6)
  kap <- pmax(eigen(Ktr, symmetric = TRUE, only.values = TRUE)$values, 0)
  list(prob = prob, nominal = Vv * R, eff_df = sum(kap / (kap + ridge)), penalty = ridge)
}

## Kundu arm for a continuous outcome (no [0,1] clamp); copied from dense_variants.R.
kundu_regress_pred <- function(conn_train, y_train, conn_test, R = 8, ridge = 1e-2) {
  Wbar <- Reduce(`+`, lapply(conn_train, abs)) / length(conn_train); diag(Wbar) <- 0
  eg <- eigen(Wbar, symmetric = TRUE); ord <- order(abs(eg$values), decreasing = TRUE)
  U <- eg$vectors[, ord[seq_len(R)], drop = FALSE]
  feat <- function(cl) t(sapply(cl, function(A) as.numeric(A %*% U)))
  Ftr <- feat(conn_train); Fte <- feat(conn_test)
  cm <- colMeans(Ftr); cs <- apply(Ftr, 2, sd); cs[cs == 0] <- 1
  Ftr <- sweep(sweep(Ftr, 2, cm, "-"), 2, cs, "/"); Fte <- sweep(sweep(Fte, 2, cm, "-"), 2, cs, "/")
  sq <- function(A, B) outer(rowSums(A^2), rowSums(B^2), "+") - 2 * A %*% t(B)
  Dtr <- sq(Ftr, Ftr); sigma2 <- median(Dtr[upper.tri(Dtr)]); sigma2 <- ifelse(sigma2 > 0, sigma2, 1)
  Ktr <- exp(-Dtr / (2 * sigma2)); Kte <- exp(-sq(Fte, Ftr) / (2 * sigma2))
  a <- solve(Ktr + ridge * diag(nrow(Ktr)), y_train - mean(y_train))
  list(pred = as.numeric(Kte %*% a) + mean(y_train))
}

## ---- Gaussian (continuous outcome) effective-df fitters ---------------------------------
## Two-sided ARD regression: ard_fit on train-standardized features, y centered (intercept = 1 df).
## MacKay df = 1 + sum gamma_j, gamma_j = 1 - lambda_j Sigma_jj, Sigma = (a X'X + diag(lambda))^-1
## recomputed at the returned (lambda, alpha_noise) on the active set (pruned features: gamma = 0).
gs_ard_reg_with_df <- function(Xs_tr, y_tr, Xs_te) {
  pp <- standardize_train_test(Xs_tr, Xs_te, y_tr)
  fit <- ard_fit(pp$X, pp$y)
  act <- which(fit$lambda < 1e4)                     # ard_fit's threshold_lambda default
  g <- 0
  if (length(act)) {
    Xa <- pp$X[, act, drop = FALSE]; la <- fit$lambda[act]
    A <- fit$alpha_noise * crossprod(Xa) + diag(la, length(la))
    Sig <- tryCatch(solve(A), error = function(e) qr.solve(A))
    g <- sum(pmax(0, 1 - la * diag(Sig)))
  }
  list(pred = as.numeric(pp$Xt %*% fit$coef + pp$ybar), nominal = ncol(Xs_tr), eff_df = 1 + g,
       penalty = NA_real_)
}
## CV ridge regression (comparators_cv.R::raw_edge_ridge_cv, dual form): df = 1 + sum d/(d + alpha_sel)
## over eigenvalues d of the n x n Gram matrix of the train-standardized design.
ridge_reg_with_df <- function(Xe_tr, y_tr, Xe_te) {
  r <- raw_edge_ridge_cv(Xe_tr, y_tr, Xe_te)
  Xs <- standardize_train_test(Xe_tr, Xe_te, y_tr)$X
  d <- pmax(eigen(tcrossprod(Xs), symmetric = TRUE, only.values = TRUE)$values, 0)
  list(pred = r$pred, nominal = ncol(Xe_tr), eff_df = 1 + sum(d / (d + r$alpha_sel)), penalty = r$alpha_sel)
}
## Kernel ridge (Kundu re-implementation, no clamp): df = trace K (K + lambda I)^-1 (+1 for the mean).
kundu_reg_with_df <- function(conn_tr, y_tr, conn_te, R = 8, ridge = 1e-2) {
  Wbar <- Reduce(`+`, lapply(conn_tr, abs)) / length(conn_tr); diag(Wbar) <- 0
  eg <- eigen(Wbar, symmetric = TRUE); U <- eg$vectors[, order(abs(eg$values), decreasing = TRUE)[seq_len(R)], drop = FALSE]
  feat <- function(cl) t(sapply(cl, function(A) as.numeric(A %*% U)))
  Ftr <- feat(conn_tr); Fte <- feat(conn_te)
  cm <- colMeans(Ftr); cs <- apply(Ftr, 2, sd); cs[cs == 0] <- 1
  Ftr <- sweep(sweep(Ftr, 2, cm, "-"), 2, cs, "/"); Fte <- sweep(sweep(Fte, 2, cm, "-"), 2, cs, "/")
  sq <- function(A, B) outer(rowSums(A^2), rowSums(B^2), "+") - 2 * A %*% t(B)
  Dtr <- sq(Ftr, Ftr); sigma2 <- median(Dtr[upper.tri(Dtr)]); sigma2 <- ifelse(sigma2 > 0, sigma2, 1)
  Ktr <- exp(-Dtr / (2 * sigma2)); Kte <- exp(-sq(Fte, Ftr) / (2 * sigma2))
  a <- solve(Ktr + ridge * diag(nrow(Ktr)), y_tr - mean(y_tr))
  kap <- pmax(eigen(Ktr, symmetric = TRUE, only.values = TRUE)$values, 0)
  list(pred = as.numeric(Kte %*% a) + mean(y_tr), nominal = nrow(conn_tr[[1]]) * R,
       eff_df = 1 + sum(kap / (kap + ridge)), penalty = ridge)
}

## ---- append rows to a CSV (creating it with header on first write) --------------------
append_csv <- function(df, path) {
  write.table(df, path, sep = ",", row.names = FALSE, col.names = !file.exists(path),
              append = file.exists(path), qmethod = "double")
}
