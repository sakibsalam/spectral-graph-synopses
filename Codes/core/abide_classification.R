## =====================================================================================
## ABIDE autism-vs-control classification (extension).
##
## Binary-outcome version of the scalar-on-network model. The manuscript notes that for
## binary outcomes the Gaussian likelihood is replaced by a logistic likelihood with the
## "same graph-spectral ARD principle ... via variational logistic regression or a
## Polya-Gamma approximation." Here we implement the Laplace/IRLS evidence-maximization
## form of ARD logistic regression -- i.e. Tipping's Relevance Vector Machine (already
## cited in the manuscript) -- on the graph-spectral features.
##
## Comparators (scalable to the p=19,900 raw edges): glmnet logistic ridge and lasso.
## A full raw-edge ARD/SSVS logistic is again computationally prohibitive at p=19,900.
##
## Metric: test-set AUC (+ accuracy), averaged over repeated random splits.
## =====================================================================================

suppressMessages(source("Codes/core/abide_pipeline.R"))  # load_abide, make_basis, build_features, helpers

auc_pred <- function(prob, label) {
  pos <- prob[label == 1]; neg <- prob[label == 0]
  if (!length(pos) || !length(neg)) return(NA_real_)
  r <- rank(c(pos, neg))
  (sum(r[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) / (length(pos) * length(neg))
}

## ---- ARD logistic (RVM-style: IRLS MAP + MacKay evidence updates) -------------------
## y in {0,1}. beta_j | omega_j ~ N(0, omega_j^{-1}); intercept unpenalized.
ard_logistic_fit <- function(X, y, max_iter = 80, tol = 1e-4,
                             threshold = 1e8, eps = 1e-10) {
  n <- nrow(X); p <- ncol(X)
  Xd <- cbind(1, X)                       # column 1 = intercept
  prec <- c(0, rep(1, p))                 # intercept unpenalized; features start at omega=1
  beta <- c(log((mean(y) + eps) / (1 - mean(y) + eps)), rep(0, p))

  for (iter in seq_len(max_iter)) {
    old <- beta
    keep <- c(1, 1 + which(prec[-1] < threshold))   # always keep intercept
    Xk <- Xd[, keep, drop = FALSE]; pk <- prec[keep]; bk <- beta[keep]

    ## inner IRLS for the MAP estimate given current precisions
    for (inner in seq_len(50)) {
      eta <- as.numeric(Xk %*% bk)
      pr <- 1 / (1 + exp(-eta)); pr <- pmin(pmax(pr, 1e-6), 1 - 1e-6)
      W <- pr * (1 - pr)
      z <- eta + (y - pr) / W
      H <- crossprod(Xk * sqrt(W)) + diag(pk, length(pk))
      bk_new <- tryCatch(solve(H, crossprod(Xk, W * z)),
                         error = function(e) qr.solve(H, crossprod(Xk, W * z)))
      bk_new <- as.numeric(bk_new)
      if (sqrt(sum((bk_new - bk)^2)) < 1e-6 * max(1, sqrt(sum(bk^2)))) { bk <- bk_new; break }
      bk <- bk_new
    }

    ## posterior covariance at the mode -> MacKay precision update
    eta <- as.numeric(Xk %*% bk)
    pr <- 1 / (1 + exp(-eta)); pr <- pmin(pmax(pr, 1e-6), 1 - 1e-6)
    W <- pr * (1 - pr)
    H <- crossprod(Xk * sqrt(W)) + diag(pk, length(pk))
    Sig <- tryCatch(solve(H), error = function(e) qr.solve(H))
    dS <- diag(Sig)

    beta <- rep(0, p + 1); beta[keep] <- bk
    feat <- keep[keep > 1]                          # feature columns currently active
    gamma <- pmax(0, 1 - prec[feat] * dS[match(feat, keep)])
    prec_new <- gamma / pmax(beta[feat]^2, eps)
    prec_new[!is.finite(prec_new)] <- threshold * 10
    prec[feat] <- prec_new
    prec[prec > threshold] <- threshold * 10

    denom <- max(1, sqrt(sum(old^2)))
    if (sqrt(sum((beta - old)^2)) / denom < tol) break
  }
  list(intercept = beta[1], coef = beta[-1], prec = prec[-1])
}

predict_logit <- function(intercept, coef, Xt) {
  eta <- intercept + as.numeric(Xt %*% coef)
  1 / (1 + exp(-eta))
}

## standardize features only (do NOT center the binary outcome)
std_feats <- function(X, Xt) {
  ctr <- colMeans(X); scl <- apply(X, 2, sd); scl[scl == 0 | is.na(scl)] <- 1
  list(X = sweep(sweep(X, 2, ctr, "-"), 2, scl, "/"),
       Xt = sweep(sweep(Xt, 2, ctr, "-"), 2, scl, "/"))
}

graph_spectral_ard_logit <- function(Xs_tr, y, Xs_te) {
  t0 <- proc.time()[3]
  pp <- std_feats(Xs_tr, Xs_te)
  fit <- ard_logistic_fit(pp$X, y)
  prob <- predict_logit(fit$intercept, fit$coef, pp$Xt)
  list(prob = prob, seconds = proc.time()[3] - t0, nsel = sum(abs(fit$coef) > 1e-6))
}

raw_edge_glmnet_logit <- function(Xe_tr, y, Xe_te, alpha) {
  if (!has_glmnet) return(NULL)
  t0 <- proc.time()[3]
  pp <- std_feats(Xe_tr, Xe_te)
  cv <- glmnet::cv.glmnet(pp$X, y, family = "binomial", alpha = alpha, nfolds = 5,
                          standardize = FALSE, type.measure = "auc")
  prob <- as.numeric(predict(cv, newx = pp$Xt, s = "lambda.min", type = "response"))
  b <- as.numeric(coef(cv, s = "lambda.min"))[-1]
  list(prob = prob, seconds = proc.time()[3] - t0, nsel = sum(abs(b) > 1e-10))
}

## ---- repeated-split classification experiment --------------------------------------

run_abide_dx <- function(K = 10, basis = "adjacency", n_splits = 10, test_frac = 0.30,
                         out_csv = "results/abide_dx_results.csv") {
  d <- load_abide()
  ## DX_GROUP: 1 = autism, 2 = control -> y = 1 for autism
  y_all <- as.integer(d$subjects$DX_GROUP == 1)
  conn <- d$conn; N <- length(conn)
  cat(sprintf("ABIDE DX classification: N=%d, autism=%d control=%d; basis=%s K=%d, %d splits\n",
              N, sum(y_all), sum(1 - y_all), basis, K, n_splits))

  acc <- list()
  for (s in seq_len(n_splits)) {
    set.seed(200 + s)
    test_idx <- sample(seq_len(N), round(test_frac * N))
    train_idx <- setdiff(seq_len(N), test_idx)
    y_tr <- y_all[train_idx]; y_te <- y_all[test_idx]

    ft <- build_features(conn, train_idx, K, basis = basis)
    Xs_tr <- ft$X_spec[train_idx, ]; Xs_te <- ft$X_spec[test_idx, ]
    Xe_tr <- ft$X_edge[train_idx, ]; Xe_te <- ft$X_edge[test_idx, ]

    gs <- graph_spectral_ard_logit(Xs_tr, y_tr, Xs_te)
    rg <- raw_edge_glmnet_logit(Xe_tr, y_tr, Xe_te, alpha = 0)   # ridge
    ls <- raw_edge_glmnet_logit(Xe_tr, y_tr, Xe_te, alpha = 1)   # lasso

    acc[[s]] <- data.frame(
      split = s,
      gs_ard_auc = auc_pred(gs$prob, y_te),
      ridge_auc  = if (!is.null(rg)) auc_pred(rg$prob, y_te) else NA_real_,
      lasso_auc  = if (!is.null(ls)) auc_pred(ls$prob, y_te) else NA_real_,
      gs_acc     = mean((gs$prob > 0.5) == y_te),
      ridge_acc  = if (!is.null(rg)) mean((rg$prob > 0.5) == y_te) else NA_real_,
      gs_nsel    = gs$nsel, gs_dim = ft$q, edge_dim = ft$p)
  }
  df <- do.call(rbind, acc)
  dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
  write.csv(df, out_csv, row.names = FALSE)

  msd <- function(x) sprintf("%.3f +/- %.3f", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
  cat("\n=========== ABIDE autism-vs-control, mean test AUC (higher better) ===========\n")
  cat(sprintf("  graph-spectral ARD-logit [%s, K=%d, dim=%d] : AUC %s | acc %s | ~%.0f modes\n",
              basis, K, df$gs_dim[1], msd(df$gs_ard_auc), msd(df$gs_acc), mean(df$gs_nsel)))
  cat(sprintf("  raw-edge logistic ridge  (dim=%d)           : AUC %s | acc %s\n",
              df$edge_dim[1], msd(df$ridge_auc), msd(df$ridge_acc)))
  cat(sprintf("  raw-edge logistic lasso  (dim=%d)           : AUC %s\n",
              df$edge_dim[1], msd(df$lasso_auc)))
  cat("  (chance AUC = 0.500)\n")
  cat("Saved to", out_csv, "\n")
  invisible(df)
}

if (sys.nframe() == 0) {
  run_abide_dx()
}
