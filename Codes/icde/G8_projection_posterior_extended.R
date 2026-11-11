## =====================================================================================
## G8: extended projection-posterior experiment (extends G7) on four collections:
##   sch100, power229 (age_continuous, make_split(subj, s)) and abide, acpi (split as in G2).
##   Rscript Codes/icde/G8_projection_posterior_extended.R [--smoke]
## One fit per split (two-sided core K = 10, q = 55, Gaussian ridge posterior N(m, S)); every
## query reads only the cache (m, S, G, Q, scaling) plus held-out cores:
##   (1) l1 projection to k in {5,10,20} nonzeros  vs  lasso refit on the cores
##   (2) group (mode-level) projection to g in {3,5,7} modes by greedy backward elimination on the
##       cached quadratic  vs  ridge refit on the same modes' core columns
##   (3) 90% predictive intervals from the cache  vs  B-resample bootstrap of dual edge ridge
## Outputs: results/icde/projpost_ext_accuracy.csv, projpost_ext_queries.csv
## =====================================================================================
source("Codes/icde/icde_common.R")
smoke <- has_flag("--smoke"); n_splits <- if (smoke) 1 else 10; B_boot <- 200
o1 <- out_path("projpost_ext_accuracy.csv", smoke); o2 <- out_path("projpost_ext_queries.csv", smoke)
for (f in c(o1, o2)) if (file.exists(f)) file.remove(f)
colls <- if (smoke) c("sch100", "abide") else c("sch100", "power229", "abide", "acpi")
soft <- function(z, t) sign(z) * pmax(abs(z) - t, 0)
proj_l1 <- function(theta, G, lam, b = theta, iters = 200) {
  r <- as.numeric(G %*% theta)
  for (it in seq_len(iters)) { b0 <- b
    for (j in seq_along(b)) b[j] <- soft(r[j] - sum(G[j, -j] * b[-j]), lam / 2) / G[j, j]
    if (max(abs(b - b0)) < 1e-7) break }
  b
}
lam_for_k <- function(theta, G, k) {
  lo <- 0; hi <- 2 * max(abs(G %*% theta))
  for (i in 1:40) { mid <- (lo + hi) / 2; if (sum(proj_l1(theta, G, mid) != 0) > k) lo <- mid else hi <- mid }
  hi
}
idx <- which(upper.tri(diag(10), diag = TRUE), arr.ind = TRUE)      # coefficient j <-> (idx[j,1], idx[j,2])
coefs_in <- function(R) which(idx[, 1] %in% R & idx[, 2] %in% R)
restricted <- function(m, G, Gm, R) { j <- coefs_in(R); b <- numeric(length(m)); b[j] <- solve(G[j, j, drop = FALSE], Gm[j]); b }
quad_loss <- function(m, G, b) { e <- m - b; as.numeric(crossprod(e, G %*% e)) }
backward_modes <- function(m, G, gs = c(7, 5, 3)) {                  # greedy backward elimination of modes
  Gm <- as.numeric(G %*% m); R <- 1:10; out <- list()
  while (length(R) > min(gs)) {
    ls <- sapply(R, function(r) quad_loss(m, G, restricted(m, G, Gm, setdiff(R, r))))
    R <- R[-which.min(ls)]
    if (length(R) %in% gs) out[[as.character(length(R))]] <- list(R = R, b = restricted(m, G, Gm, R))
  }
  out
}
## dual ridge with a FIXED penalty; standardization from the supplied training rows only
dual_fixed <- function(Xtr, ytr, Xte, a) {
  n <- nrow(Xtr); mu <- colMeans(Xtr); sdv <- sqrt(pmax(colMeans(Xtr^2) - mu^2, 0) * n / (n - 1)); sdv[sdv == 0] <- 1
  Xs <- sweep(sweep(Xtr, 2, mu, "-"), 2, sdv, "/"); ybar <- mean(ytr); yc <- ytr - ybar
  al <- solve(tcrossprod(Xs) + a * diag(n), yc); beta <- as.numeric(crossprod(Xs, al))
  Xts <- sweep(sweep(Xte, 2, mu, "-"), 2, sdv, "/")
  list(pred = as.numeric(Xts %*% beta) + ybar, fitted = as.numeric(Xs %*% beta) + ybar, Xs = Xs)
}

for (coll in colls) {
  if (coll %in% c("sch100", "power229")) {
    d <- load_set(coll, "age_continuous"); conn <- d$conn; y <- d$y; subj <- d$subj; V <- d$V; own_subject <- FALSE
  } else {
    d <- load_abide(file.path("data", coll)); y <- d$subjects$AGE_AT_SCAN; keep <- which(is.finite(y))
    conn <- d$conn[keep]; y <- y[keep]; subj <- seq_along(y); V <- d$V; own_subject <- TRUE
  }
  N <- length(conn); Xe <- build_edge_features(conn)$X_edge
  cat(sprintf("[%s] N=%d V=%d p=%d\n", coll, N, V, ncol(Xe))); flush.console()
  for (s in seq_len(n_splits)) {
    if (own_subject) { set.seed(100 + s); te <- sample(seq_len(N), round(0.30 * N)); tr <- setdiff(seq_len(N), te)
    } else { sp <- make_split(subj, s); tr <- sp$tr; te <- sp$te }
    n <- length(tr); fold_tr <- if (own_subject) subject_folds(seq_along(tr)) else subject_folds(subj[tr])
    Q <- template_Q(conn, tr, 10); XT <- feat_twosided(conn, Q)
    ## ---------- fit once (as G7) ----------
    pp <- standardize_train_test(XT[tr, ], XT[te, ], y[tr]); X <- pp$X; yc <- pp$y; q <- ncol(X)
    a <- ridge_dual_cv_folds(XT[tr, ], y[tr], XT[te, ], fold_tr)$alpha_sel
    A <- solve(crossprod(X) + a * diag(q)); m <- as.numeric(A %*% crossprod(X, yc))
    dfr <- sum(diag(X %*% A %*% t(X))); s2 <- sum((yc - X %*% m)^2) / (n - dfr - 1); S <- s2 * A
    G <- crossprod(X) / n
    cache_bytes <- 8 * (length(m) + length(S) + length(G) + length(Q) + 2 * q + 2)
    Xt <- pp$Xt; yte <- y[te]; ybar <- pp$ybar
    acc <- list()
    addacc <- function(method, k, g, b, pred = NULL) {
      if (is.null(pred)) pred <- as.numeric(Xt %*% b) + ybar
      acc[[length(acc) + 1]] <<- data.frame(collection = coll, split = s, method = method, k = k, g = g,
        nnz = sum(b != 0), spearman = spearman(yte, pred), rmse = rmse(yte, pred)) }
    addacc("dense_mean", q, NA, m)
    ## ---------- (1) l1 projection vs lasso refit ----------
    cv <- glmnet::glmnet(X, yc, alpha = 1, standardize = FALSE, intercept = FALSE, nlambda = 200)
    for (k in c(5, 10, 20)) {
      lam <- lam_for_k(m, G, k); b <- proj_l1(m, G, lam); addacc("projection_l1", k, NA, b)
      j <- max(which(cv$df <= k)); addacc("lasso_refit", k, NA, as.numeric(cv$beta[, j]))
    }
    ## ---------- (2) group (mode-level) projection vs ridge refit on the same modes ----------
    bm <- backward_modes(m, G)
    for (g in c(3, 5, 7)) { o <- bm[[as.character(g)]]; cols <- coefs_in(o$R)
      addacc("projection_modes", length(cols), g, o$b)
      rr <- ridge_dual_cv_folds(XT[tr, cols, drop = FALSE], y[tr], XT[te, cols, drop = FALSE], fold_tr)
      b0 <- numeric(q); b0[cols] <- 1   # nnz bookkeeping only
      acc[[length(acc) + 1]] <- data.frame(collection = coll, split = s, method = "ridge_refit_same_modes", k = length(cols), g = g,
        nnz = length(cols), spearman = spearman(yte, rr$pred), rmse = rmse(yte, rr$pred)) }
    append_csv(do.call(rbind, acc), o1)
    ## ---------- (3) predictive intervals: cache vs bootstrap of edge ridge ----------
    z <- qnorm(0.95)
    t_cache <- system.time({ mu <- as.numeric(Xt %*% m) + ybar; sdp <- sqrt(s2 + rowSums((Xt %*% S) * Xt))
      cover <- mean(abs(yte - mu) <= z * sdp); wid <- mean(2 * z * sdp) })[3]
    er <- ridge_dual_cv_folds(Xe[tr, ], y[tr], Xe[te, ], fold_tr); ae <- er$alpha_sel
    edge_sp <- spearman(yte, er$pred)
    fullfit <- dual_fixed(Xe[tr, ], y[tr], Xe[te, ], ae)
    dd <- pmax(eigen(tcrossprod(fullfit$Xs), symmetric = TRUE, only.values = TRUE)$values, 0)
    dfe <- sum(dd / (dd + ae)); rv <- sum((y[tr] - fullfit$fitted)^2) / (n - dfe - 1)
    set.seed(9000 + s)
    t_boot <- system.time({
      P <- matrix(NA_real_, B_boot, length(te))
      for (bb in seq_len(B_boot)) { ii <- sample.int(n, n, replace = TRUE)
        P[bb, ] <- dual_fixed(Xe[tr[ii], , drop = FALSE], y[tr[ii]], Xe[te, , drop = FALSE], ae)$pred }
      mb <- colMeans(P); sb <- sqrt(pmax(apply(P, 2, var) + rv, 0))
      bcover <- mean(abs(yte - mb) <= z * sb); bwid <- mean(2 * z * sb) })[3]
    append_csv(data.frame(collection = coll, split = s, n_train = n, V = V, cache_bytes = cache_bytes,
      edge_store_bytes = 8 * n * ncol(Xe), cache_cover90 = cover, cache_width90 = wid, cache_interval_seconds = t_cache,
      boot_cover90 = bcover, boot_width90 = bwid, boot_seconds = t_boot, boot_B = B_boot, edge_ridge_spearman = edge_sp,
      dense_spearman = spearman(yte, as.numeric(Xt %*% m)), dense_rmse = rmse(yte, as.numeric(Xt %*% m) + ybar)), o2)
    cat(coll, "split", s, "done; boot", round(t_boot, 1), "s\n"); flush.console()
  }
}
cat("Wrote", o1, "and", o2, "\n")
