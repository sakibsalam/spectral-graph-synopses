## =====================================================================================
## G6: effective degrees of freedom of logistic edge ridge on the binary task when the inner
## cross-validation folds are formed over subjects instead of graphs.
##   Rscript Codes/icde/G6_df_grouped_t1.R --atlas=sch100 [--smoke]
## Same estimator as ridge_logit_with_df in icde_common.R (glmnet, deviance criterion, default
## path; df = 1 + sum d/(d + c) from the weighted dual Gram matrix at convergence); the only
## change is the foldid argument. Edge lasso (AUC criterion) is refitted the same way; its size is the
## number of nonzero coefficients. The two-sided ARD fit has no inner CV and is refitted here
## on the same splits for a paired comparison.
## Output: results/icde/df_grouped_<atlas>.csv
## =====================================================================================
source("Codes/icde/icde_common.R")
smoke <- has_flag("--smoke"); atlas <- arg_val("--atlas", "sch100"); n_splits <- if (smoke) 1 else 10
out <- out_path(sprintf("df_grouped_%s.csv", atlas), smoke); if (file.exists(out)) file.remove(out)

ridge_logit_df_folds <- function(Xe_tr, y_tr, Xe_te, foldid, max_iter = 30) {
  pp <- std_feats(Xe_tr, Xe_te); Xs <- pp$X; Xt <- pp$Xt; n <- nrow(Xs); p <- ncol(Xs)
  cv <- glmnet::cv.glmnet(Xs, y_tr, family = "binomial", alpha = 0, foldid = foldid, standardize = FALSE, type.measure = "deviance")
  c_pen <- n * cv$lambda.min
  beta0 <- log((mean(y_tr) + 1e-9) / (1 - mean(y_tr) + 1e-9)); beta <- rep(0, p)
  for (it in seq_len(max_iter)) {
    eta <- beta0 + as.numeric(Xs %*% beta); pr <- pmin(pmax(1 / (1 + exp(-eta)), 1e-6), 1 - 1e-6); w <- pr * (1 - pr)
    z <- eta + (y_tr - pr) / w; wbar <- sum(w); xbar <- as.numeric(crossprod(w, Xs)) / wbar; zbar <- sum(w * z) / wbar
    Xw <- sweep(Xs, 2, xbar, "-") * sqrt(w); zw <- (z - zbar) * sqrt(w)
    sol <- solve(tcrossprod(Xw) + c_pen * diag(n), zw)
    bn <- as.numeric(crossprod(Xw, sol)); b0n <- zbar - sum(xbar * bn)
    delta <- sqrt(sum((bn - beta)^2) + (b0n - beta0)^2); beta <- bn; beta0 <- b0n
    if (delta < 1e-6 * max(1, sqrt(sum(beta^2)))) break
  }
  eta <- beta0 + as.numeric(Xs %*% beta); pr <- pmin(pmax(1 / (1 + exp(-eta)), 1e-6), 1 - 1e-6); w <- pr * (1 - pr)
  xbar <- as.numeric(crossprod(w, Xs)) / sum(w); Xw <- sweep(Xs, 2, xbar, "-") * sqrt(w)
  d <- pmax(eigen(tcrossprod(Xw), symmetric = TRUE, only.values = TRUE)$values, 0)
  list(prob = 1 / (1 + exp(-(beta0 + as.numeric(Xt %*% beta)))), eff_df = 1 + sum(d / (d + c_pen)),
       at_path_end = cv$lambda.min == min(cv$lambda), at_path_start = cv$lambda.min == max(cv$lambda))
}

d <- load_set(atlas, "nbk_0v1"); conn <- d$conn; y <- d$y; subj <- d$subj
Xe <- build_edge_features(conn)$X_edge
for (s in seq_len(n_splits)) {
  sp <- make_split(subj, s); tr <- sp$tr; te <- sp$te
  set.seed(1); f_graph <- sample(rep(seq_len(5), length.out = length(tr)))
  f_subj <- subject_folds(subj[tr])
  Xs <- feat_twosided(conn, template_Q(conn, tr, 10)); a <- gs_ard_logit_with_df(Xs[tr, ], y[tr], Xs[te, ])
  rows <- list(data.frame(atlas = atlas, split = s, method = "two_sided_ard", folds = "none", eff_df = a$eff_df,
                          auc = auc_pred(a$prob, y[te]), pred_sd = sd(a$prob), at_path_end = NA, at_path_start = NA))
  for (fold in c("graph", "subject")) {
    r <- ridge_logit_df_folds(Xe[tr, ], y[tr], Xe[te, ], if (fold == "graph") f_graph else f_subj)
    rows[[length(rows) + 1]] <- data.frame(atlas = atlas, split = s, method = "edge_ridge_logit", folds = fold, eff_df = r$eff_df,
      auc = auc_pred(r$prob, y[te]), pred_sd = sd(r$prob), at_path_end = r$at_path_end, at_path_start = r$at_path_start)
  }
  for (fold in c("graph", "subject")) {                      # edge lasso (glmnet, AUC criterion), size = nonzero count
    pp <- std_feats(Xe[tr, ], Xe[te, ])
    cv <- glmnet::cv.glmnet(pp$X, y[tr], family = "binomial", alpha = 1, foldid = if (fold == "graph") f_graph else f_subj,
                            standardize = FALSE, type.measure = "auc")
    pr <- as.numeric(predict(cv, newx = pp$Xt, s = "lambda.min", type = "response"))
    rows[[length(rows) + 1]] <- data.frame(atlas = atlas, split = s, method = "edge_lasso_logit", folds = fold,
      eff_df = sum(abs(as.numeric(coef(cv, s = "lambda.min"))[-1]) > 1e-10), auc = auc_pred(pr, y[te]), pred_sd = sd(pr),
      at_path_end = cv$lambda.min == min(cv$lambda), at_path_start = cv$lambda.min == max(cv$lambda))
  }
  append_csv(do.call(rbind, rows), out); cat("  split", s, "done\n"); flush.console()
}
cat("Wrote", out, "\n")
