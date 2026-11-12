## =====================================================================================
## Model degrees-of-freedom comparison: nominal parameter count vs EFFECTIVE degrees of
## freedom, per model, on real fits (nbk_0v1 at V=100 / Schaefer-100, p=4,950).
##
## Motivation: "55 vs 26,106 parameters" (the parsimony headline elsewhere in this project)
## counts the length of the coefficient vector each model CARRIES. That is a different
## quantity from how many degrees of freedom the fit actually USES once shrinkage/selection
## is accounted for -- e.g. raw-edge SSVS carries p=4,950 nominal coefficients but its
## effective df (sum of posterior inclusion probabilities) is far smaller, while ridge with
## a very small CV-selected penalty can use almost all p nominal dimensions. This script
## computes both, for every model in the real-data comparison, and reports which of each
## model's dimensions/knobs a user actually controls versus what the fit learns.
##
## Uses ONLY data/osf_schaefer100/master.rds (V=100, p=4,950) -- 5x cheaper than the V=229
## Power-229 master, and the df question is structural, not resolution-specific, per the
## task instructions.
##
## Ownership: this file, results/model_df_comparison.csv, results/figures/fig14_*.png, and
## docs/model_side_by_side.md are owned by this script/session. Nothing under Codes/core,
## Codes/benchmarks, Codes/analysis/outcome_explorer.R or multi_resolution_common.R,
## or any pre-existing results/*.csv is modified --
## effective-df estimators that need internals not RETURNED by an existing fitting function
## (e.g. ARD's posterior gamma) are recomputed HERE from that function's own return values
## (coef, prec) rather than by editing the source function.
## =====================================================================================

suppressMessages({
  source("Codes/analysis/multi_resolution_common.R")   # -> outcome_explorer.R -> abide_classification.R
                                                         # -> abide_pipeline.R -> network_bayes_simulation.R
                                                         # -> comparators_cv.R, features_cv.R
                                                         # brings in: build_features, build_spec_features,
                                                         # ard_logistic_fit, predict_logit, std_feats, auc_pred,
                                                         # graph_spectral_ard_logit, raw_edge_glmnet_logit,
                                                         # spike_slab_ssvs, make_set_m
  source("Codes/benchmarks/competitor_benchmarks.R")    # guha_horseshoe_logit (Guha proxy)
})

OUT_CSV <- "results/model_df_comparison.csv"
FIG_DIR <- "results/figures"
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(OUT_CSV), showWarnings = FALSE, recursive = TRUE)

m <- readRDS("data/osf_schaefer100/master.rds")
cat(sprintf("Loaded Schaefer-100 master: %d connectomes, V=%d\n", length(m$conn), m$V))

set_ <- make_set_m(m, conds = c("nbk_0back", "nbk_1back"))
conn <- set_$conn; meta <- set_$meta
y_all <- as.integer(as.character(meta$COND) == "nbk_1back")
subj <- meta$SUBJECT; usub <- unique(subj)
V <- nrow(conn[[1]])
cat(sprintf("nbk_0v1: %d rows, %d subjects, V=%d, p=%d\n", length(conn), length(usub), V, V*(V-1)/2))

n_splits <- 10
seed0 <- 900
test_frac <- 0.30
K <- 10
basis <- "adjacency"

## ---------------------------------------------------------------------------------------
## Job 1 helpers: one fitting+df function per model. Each returns at least
## list(prob=<test predictions>, nominal=<int>, eff_df=<numeric>).
## ---------------------------------------------------------------------------------------

## ---- ours: graph-spectral ARD-logit. MacKay's sum(gamma) is already what ard_logistic_fit
## computes internally (abide_classification.R:64, gamma <- pmax(0, 1 - prec[feat]*dS[...]))
## but does not return. Recomputed here from the function's own OUTPUT (intercept, coef, prec)
## by reconstructing the IRLS Hessian at the reported mode -- not a re-derivation, an exact
## recomputation of the same quantity the fit already computed and discarded.
gs_ard_logit_with_df <- function(Xs_tr, y_tr, Xs_te) {
  pp <- std_feats(Xs_tr, Xs_te)
  fit <- ard_logistic_fit(pp$X, y_tr)
  prob <- predict_logit(fit$intercept, fit$coef, pp$Xt)

  Xd <- cbind(1, pp$X)
  beta_full <- c(fit$intercept, fit$coef)
  eta <- as.numeric(Xd %*% beta_full)
  pr <- 1 / (1 + exp(-eta)); pr <- pmin(pmax(pr, 1e-6), 1 - 1e-6)
  W <- pr * (1 - pr)
  prec_full <- c(0, fit$prec)                      # intercept unpenalized, as in the original fit
  H <- crossprod(Xd * sqrt(W)) + diag(prec_full, length(prec_full))
  Sig <- tryCatch(solve(H), error = function(e) qr.solve(H))
  dS <- diag(Sig)
  gamma <- pmax(0, 1 - prec_full * dS)
  eff_df <- sum(gamma)                              # includes the always-kept intercept term

  list(prob = prob, nominal = ncol(Xs_tr), eff_df = eff_df,
       nsel = sum(abs(fit$coef) > 1e-6))
}

## ---- raw-edge ridge (logistic): GLM ridge effective df. lambda selected by glmnet's own
## CV ridge logistic (5-fold, deviance); converted to the equivalent IRLS ridge penalty
## c = n*lambda (glmnet minimizes (1/n)*deviance + lambda*(1/2)||beta||^2; the IRLS quadratic
## approximation to (1/n)*deviance is (1/(2n))*sum(w_i*(z_i-x_i'beta)^2); multiplying through
## by 2n gives sum(w_i*(...)^2) + (n*lambda)*||beta||^2, i.e. c = n*lambda in a standard
## weighted-ridge objective). Intercept left unpenalized via the standard weighted-centering
## trick (absorb the intercept by centering X and z at their w-weighted means each iteration).
## Effective df = 1 (unpenalized intercept) + trace(Xw(Xw'Xw+cI)^-1 Xw') at convergence, using
## the DUAL (n x n) form so p=4,950 never needs a p x p solve -- same trick as
## comparators_cv.R::raw_edge_ridge_cv(), generalized to the weighted/logistic case.
ridge_logit_with_df <- function(Xe_tr, y_tr, Xe_te, nfolds = 5, seed = 1, max_iter = 30) {
  pp <- std_feats(Xe_tr, Xe_te)
  Xs <- pp$X; Xt <- pp$Xt
  n <- nrow(Xs); p <- ncol(Xs)

  set.seed(seed)
  cv <- glmnet::cv.glmnet(Xs, y_tr, family = "binomial", alpha = 0, nfolds = nfolds,
                          standardize = FALSE, type.measure = "deviance")
  lambda_sel <- cv$lambda.min
  c_pen <- n * lambda_sel

  beta0 <- log((mean(y_tr) + 1e-9) / (1 - mean(y_tr) + 1e-9))
  beta <- rep(0, p)
  it_used <- 0
  for (it in seq_len(max_iter)) {
    it_used <- it
    eta <- beta0 + as.numeric(Xs %*% beta)
    pr <- 1 / (1 + exp(-eta)); pr <- pmin(pmax(pr, 1e-6), 1 - 1e-6)
    w <- pr * (1 - pr)
    z <- eta + (y_tr - pr) / w

    wbar <- sum(w)
    xbar <- as.numeric(crossprod(w, Xs)) / wbar
    zbar <- sum(w * z) / wbar
    Xc <- sweep(Xs, 2, xbar, "-")
    zc <- z - zbar

    sw <- sqrt(w)
    Xw <- Xc * sw
    zw <- zc * sw
    G <- tcrossprod(Xw)
    sol <- tryCatch(solve(G + c_pen * diag(n), zw),
                    error = function(e) qr.solve(G + c_pen * diag(n), zw))
    beta_new <- as.numeric(crossprod(Xw, sol))
    beta0_new <- zbar - sum(xbar * beta_new)

    delta <- sqrt(sum((beta_new - beta)^2) + (beta0_new - beta0)^2)
    beta <- beta_new; beta0 <- beta0_new
    if (delta < 1e-6 * max(1, sqrt(sum(beta^2)))) break
  }

  ## effective df at the converged fit
  eta <- beta0 + as.numeric(Xs %*% beta)
  pr <- 1 / (1 + exp(-eta)); pr <- pmin(pmax(pr, 1e-6), 1 - 1e-6)
  w <- pr * (1 - pr)
  wbar <- sum(w); xbar <- as.numeric(crossprod(w, Xs)) / wbar
  Xc <- sweep(Xs, 2, xbar, "-")
  Xw <- Xc * sqrt(w)
  G <- tcrossprod(Xw)
  eg <- eigen(G, symmetric = TRUE, only.values = TRUE)
  d <- pmax(eg$values, 0)
  eff_df <- 1 + sum(d / (d + c_pen))

  eta_te <- beta0 + as.numeric(Xt %*% beta)
  prob_te <- 1 / (1 + exp(-eta_te))

  list(prob = prob_te, nominal = p, eff_df = eff_df,
       lambda_sel = lambda_sel, c_pen = c_pen, iters = it_used)
}

## ---- Kundu: kernel ridge on the node-latent embedding. df = trace(K(K+lambda*I)^-1) via the
## eigenvalues of the n x n RBF Gram matrix. Mirrors kundu_latent_nonlinear() in
## Codes/benchmarks/competitor_benchmarks.R (same R, bandwidth heuristic, ridge=1e-2) but also
## returns the ingredients (U, train/all-subject latent scores, kernel eigenvalues) needed for
## Job 2's pair plots, which the original function does not expose.
kundu_logit_with_df <- function(conn_tr, y_tr, conn_te, R = 8, ridge = 1e-2) {
  Vv <- nrow(conn_tr[[1]])
  Wbar <- Reduce(`+`, lapply(conn_tr, abs)) / length(conn_tr); diag(Wbar) <- 0
  eg <- eigen(Wbar, symmetric = TRUE); ord <- order(abs(eg$values), decreasing = TRUE)
  U <- eg$vectors[, ord[seq_len(R)], drop = FALSE]

  feat <- function(cl) t(sapply(cl, function(A) as.numeric(A %*% U)))
  Ftr <- feat(conn_tr); Fte <- feat(conn_te)
  cm <- colMeans(Ftr); cs <- apply(Ftr, 2, sd); cs[cs == 0] <- 1
  Ftr_s <- sweep(sweep(Ftr, 2, cm, "-"), 2, cs, "/")
  Fte_s <- sweep(sweep(Fte, 2, cm, "-"), 2, cs, "/")

  sq <- function(A, B) outer(rowSums(A^2), rowSums(B^2), "+") - 2 * A %*% t(B)
  Dtr <- sq(Ftr_s, Ftr_s); sigma2 <- median(Dtr[upper.tri(Dtr)]); sigma2 <- ifelse(sigma2 > 0, sigma2, 1)
  Ktr <- exp(-Dtr / (2 * sigma2))
  Kte <- exp(-sq(Fte_s, Ftr_s) / (2 * sigma2))
  yc <- y_tr - mean(y_tr)
  alpha_vec <- solve(Ktr + ridge * diag(nrow(Ktr)), yc)
  pred <- as.numeric(Kte %*% alpha_vec) + mean(y_tr)
  prob <- pmin(pmax(pred, 1e-6), 1 - 1e-6)

  eg2 <- eigen(Ktr, symmetric = TRUE, only.values = TRUE)
  kap <- pmax(eg2$values, 0)
  eff_df <- sum(kap / (kap + ridge))

  list(prob = prob, nominal = Vv * R, eff_df = eff_df, R = R, sigma2 = sigma2,
       U = U, Ftr = Ftr, Fte = Fte)
}

## ---------------------------------------------------------------------------------------
## Job 1 main loop
## ---------------------------------------------------------------------------------------

rows <- list()
row_i <- 0

for (s in seq_len(n_splits)) {
  set.seed(seed0 + s)
  te_s <- sample(usub, max(1, round(test_frac * length(usub))))
  te <- which(subj %in% te_s); tr <- which(!(subj %in% te_s))
  if (length(unique(y_all[tr])) < 2 || length(unique(y_all[te])) < 2) {
    cat(sprintf("split %d: degenerate label split -- skipping\n", s)); next
  }
  y_tr <- y_all[tr]; y_te <- y_all[te]
  cat(sprintf("\n=== split %d: n_train=%d n_test=%d ===\n", s, length(tr), length(te)))

  ft <- build_features(conn, tr, K, basis = basis)
  Xs <- ft$X_spec; Xe <- ft$X_edge

  addrow <- function(model, nominal, eff_df, df_estimator, knobs, prob, secs, secondary = NA_real_) {
    row_i <<- row_i + 1
    rows[[row_i]] <<- data.frame(
      model = model, split = s,
      nominal_params = nominal, effective_df = eff_df, secondary_df = secondary,
      df_estimator = df_estimator, controllable_knobs = knobs,
      fitted_auc = auc_pred(prob, y_te), seconds = round(secs, 2),
      stringsAsFactors = FALSE)
  }

  t0 <- proc.time()[3]
  r1 <- gs_ard_logit_with_df(Xs[tr, ], y_tr, Xs[te, ])
  addrow("gs_ard_logit", r1$nominal, r1$eff_df,
        "MacKay sum(gamma), gamma=1-prec*diag(Sigma) at the IRLS mode -- recomputed from ard_logistic_fit()'s own returned coef/prec (abide_classification.R:64), not a new derivation",
        "K (basis size) and basis type (adjacency/laplacian/nlaplacian) -- both FIXED here at K=10/adjacency, not CV-selected in this run; the ARD precisions themselves are LEARNED, not tuned",
        r1$prob, proc.time()[3] - t0)
  cat(sprintf("  gs_ard_logit:  nominal=%d  eff_df=%.2f  auc=%.3f\n", r1$nominal, r1$eff_df, auc_pred(r1$prob, y_te)))

  t0 <- proc.time()[3]
  r2 <- ridge_logit_with_df(Xe[tr, ], y_tr, Xe[te, ], seed = seed0 + s)
  addrow("raw_edge_ridge", r2$nominal, r2$eff_df,
        "GLM ridge effective df: 1 (unpenalized intercept) + trace(Xw(Xw'Xw+cI)^-1 Xw') via eigenvalues of the n x n weighted dual Gram matrix at IRLS convergence; c=n*lambda, lambda from glmnet::cv.glmnet(alpha=0, deviance, 5-fold)",
        "ridge penalty lambda -- CV-selected, not hand-tuned",
        r2$prob, proc.time()[3] - t0)
  cat(sprintf("  raw_edge_ridge: nominal=%d  eff_df=%.2f  auc=%.3f  (lambda=%.3g, c=%.2f, iters=%d)\n",
              r2$nominal, r2$eff_df, auc_pred(r2$prob, y_te), r2$lambda_sel, r2$c_pen, r2$iters))

  t0 <- proc.time()[3]
  l <- raw_edge_glmnet_logit(Xe[tr, ], y_tr, Xe[te, ], alpha = 1)
  addrow("raw_edge_lasso", ncol(Xe), l$nsel,
        "nonzero coefficient count (Zou, Hastie & Tibshirani 2007 lasso df approximation)",
        "lasso penalty lambda -- CV-selected (5-fold, type.measure=auc)",
        l$prob, proc.time()[3] - t0)
  cat(sprintf("  raw_edge_lasso: nominal=%d  eff_df=%d  auc=%.3f\n", ncol(Xe), l$nsel, auc_pred(l$prob, y_te)))

  t0 <- proc.time()[3]
  k <- kundu_logit_with_df(conn[tr], y_tr, conn[te])
  addrow("kundu_kernel_ridge", k$nominal, k$eff_df,
        "trace(K(K+lambda*I)^-1) via eigenvalues of the n x n RBF-kernel Gram matrix; lambda=1e-2 FIXED (as in the original re-implementation)",
        "rank R (fixed=8), RBF bandwidth (fixed, median-heuristic on training latent features), ridge lambda (fixed=1e-2) -- THREE knobs, none CV-tuned in this or the original benchmark",
        k$prob, proc.time()[3] - t0)
  cat(sprintf("  kundu:          nominal=%d  eff_df=%.2f  auc=%.3f  (sigma2=%.3g)\n",
              k$nominal, k$eff_df, auc_pred(k$prob, y_te), k$sigma2))

  t0 <- proc.time()[3]
  ss <- spike_slab_ssvs(Xe[tr, ], y_tr, Xe[te, ])
  pip_sum <- sum(ss$pip); pip_gt05 <- sum(ss$pip > 0.5)
  addrow("raw_edge_ssvs", ncol(Xe), pip_sum,
        "sum of posterior inclusion probabilities (the honest measure per the task instructions) -- NOT count(PIP>0.5), which is reported in secondary_df and can be near-zero even when AUC is competitive",
        "MCMC hyperparameters tau0/tau1/pi_incl/n_iter/burn -- all FIXED defaults (network_bayes_simulation.R::spike_slab_ssvs), not tuned",
        ss$pred, proc.time()[3] - t0, secondary = pip_gt05)
  cat(sprintf("  raw_edge_ssvs:  nominal=%d  eff_df(sum PIP)=%.2f  n(PIP>0.5)=%d  auc=%.3f\n",
              ncol(Xe), pip_sum, pip_gt05, auc_pred(ss$pred, y_te)))

  t0 <- proc.time()[3]
  g <- guha_horseshoe_logit(Xe[tr, ], y_tr, Xe[te, ])
  addrow("guha_proxy_horseshoe", ncol(Xe), g$n_active,
        "APPROXIMATE: active-edge count (|coef|>1e-4*max|coef|) from the final reweighted-ridge fit -- NOT a trace-based effective df (no closed form for the adaptive-reweighting heuristic); reported for completeness, flagged as approximate",
        "reweighting rounds (fixed=3), glmnet ridge lambda EACH round (CV-selected per round)",
        g$prob, proc.time()[3] - t0)
  cat(sprintf("  guha_proxy:     nominal=%d  eff_df(approx, active count)=%d  auc=%.3f\n",
              ncol(Xe), g$n_active, auc_pred(g$prob, y_te)))

  df_now <- do.call(rbind, rows)
  write.csv(df_now, OUT_CSV, row.names = FALSE)
}

cat("\n\n================ Job 1 summary (mean +/- sd across splits) ================\n")
df_all <- do.call(rbind, rows)
agg <- aggregate(cbind(nominal_params, effective_df, fitted_auc) ~ model, data = df_all,
                 FUN = function(x) sprintf("%.2f", mean(x)))
print(agg)
write.csv(df_all, OUT_CSV, row.names = FALSE)
cat("\nWrote", OUT_CSV, "\n")

## =========================================================================================
## Job 2: pair plots -- which latent directions actually carry the nbk_0v1 signal.
## Representative split = 1 (same seed/scheme as the split-1 iteration above).
## =========================================================================================
cat("\n\n================ Job 2: pair plots ================\n")

set.seed(seed0 + 1)
te1s <- sample(usub, max(1, round(test_frac * length(usub))))
te1 <- which(subj %in% te1s); tr1 <- which(!(subj %in% te1s))
y_tr1 <- y_all[tr1]

## ---- Kundu: block score for latent dim r = column mean of A_i U[,r] --------------------
## U built from the TRAIN mean network only (unsupervised); the column-mean score is then
## computed for EVERY subject (train+test) so the pair plot shows the whole nbk_0v1 sample.
R <- 8
Wbar1 <- Reduce(`+`, lapply(conn[tr1], abs)) / length(tr1); diag(Wbar1) <- 0
egW <- eigen(Wbar1, symmetric = TRUE); ordW <- order(abs(egW$values), decreasing = TRUE)
U1 <- egW$vectors[, ordW[seq_len(R)], drop = FALSE]

kundu_scores_all <- sapply(seq_len(R), function(r) sapply(conn, function(A) mean(as.numeric(A %*% U1[, r]))))
colnames(kundu_scores_all) <- paste0("latent_", seq_len(R))

## association of each latent block with y: |t-statistic| on TRAIN subjects only (no test leakage
## into the ranking), then the ranked dims are plotted for the full sample.
tstat_k <- sapply(seq_len(R), function(r) {
  s_tr <- kundu_scores_all[tr1, r]
  g1 <- s_tr[y_tr1 == 1]; g0 <- s_tr[y_tr1 == 0]
  tt <- try(t.test(g1, g0), silent = TRUE)
  if (inherits(tt, "try-error")) 0 else abs(tt$statistic)
})
top5_k <- order(tstat_k, decreasing = TRUE)[1:5]
cat("Kundu top-5 latent dims by |t| (train-only ranking):\n")
print(data.frame(latent_dim = top5_k, abs_t = round(tstat_k[top5_k], 3)))

col_y <- ifelse(y_all == 1, "#d1495b", "#2e86ab")

png(file.path(FIG_DIR, "fig14_kundu_pairplot.png"), width = 1100, height = 1100, res = 130)
pairs(kundu_scores_all[, top5_k, drop = FALSE], col = col_y, pch = 19, cex = 0.6,
     main = "Kundu (one-sided reduction, L=A U): top-5 latent dims by |t| vs nbk_0v1\nscore = column mean of A_i U[,r]  |  red=1back  blue=0back")
dev.off()
cat("Wrote", file.path(FIG_DIR, "fig14_kundu_pairplot.png"), "\n")

## ---- Ours: rank vech(U'AU) spectral features by ARD relevance (smallest learned precision) --
sf1 <- build_spec_features(conn, tr1, K = K, basis = basis)
Xs_all <- sf1$X_spec
pp1 <- std_feats(Xs_all[tr1, , drop = FALSE], Xs_all)   # standardize ALL subjects by TRAIN moments
fit1 <- ard_logistic_fit(pp1$X, y_tr1)
prec_rank <- order(fit1$prec)                           # ascending: smallest precision = most relevant
top5_o <- prec_rank[1:5]
cat("\nOurs top-5 spectral features by ARD relevance (smallest learned precision):\n")
print(data.frame(feature_idx = top5_o, precision = signif(fit1$prec[top5_o], 4)))

ours_scores_all <- pp1$Xt[, top5_o, drop = FALSE]
colnames(ours_scores_all) <- paste0("mode_", top5_o)

png(file.path(FIG_DIR, "fig14_ours_pairplot.png"), width = 1100, height = 1100, res = 130)
pairs(ours_scores_all, col = col_y, pch = 19, cex = 0.6,
     main = "Ours (two-sided reduction, vech(U'AU)): top-5 spectral features by ARD relevance\n(smallest learned precision)  |  red=1back  blue=0back")
dev.off()
cat("Wrote", file.path(FIG_DIR, "fig14_ours_pairplot.png"), "\n")

## ---- side-by-side geometry comparison: top-2 dims of each method -----------------------
png(file.path(FIG_DIR, "fig14_geometry_comparison.png"), width = 1600, height = 800, res = 130)
par(mfrow = c(1, 2))
plot(kundu_scores_all[, top5_k[1]], kundu_scores_all[, top5_k[2]], col = col_y, pch = 19, cex = 0.7,
    xlab = sprintf("Kundu latent_%d", top5_k[1]), ylab = sprintf("Kundu latent_%d", top5_k[2]),
    main = "Kundu: one-sided reduction L=AU\n(V x R, node identity retained)")
legend("topright", legend = c("nbk_1back (y=1)", "nbk_0back (y=0)"), col = c("#d1495b", "#2e86ab"),
      pch = 19, bty = "n", cex = 0.8)
plot(ours_scores_all[, 1], ours_scores_all[, 2], col = col_y, pch = 19, cex = 0.7,
    xlab = colnames(ours_scores_all)[1], ylab = colnames(ours_scores_all)[2],
    main = "Ours: two-sided reduction C=U'AU\n(K x K core, node identity discarded)")
dev.off()
par(mfrow = c(1, 1))
cat("Wrote", file.path(FIG_DIR, "fig14_geometry_comparison.png"), "\n")

cat("\n\nAll done.\n")
