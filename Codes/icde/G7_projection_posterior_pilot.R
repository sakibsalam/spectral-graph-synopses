## =====================================================================================
## G7 (pilot): fit once, query many times from a cached posterior, without the raw graphs.
##   Rscript Codes/icde/G7_projection_posterior_pilot.R [--smoke]
## Fit once per split: a Gaussian ridge posterior N(m, S) on the two-sided core (K = 10, q = 55),
## penalty by subject-fold CV, noise variance from the residuals. Cached state: m, S, the Gram
## matrix G = X'X / n of the standardized training cores, the V x K basis, and the scaling
## constants. Every query below reads only that cache:
##   Q1 sparse summaries: lasso projection T_lambda(theta) = argmin_b (theta-b)'G(theta-b) + lambda |b|_1
##      of the posterior mean and of 500 draws -> projected predictor with k nonzero mode pairs,
##      and posterior inclusion probabilities per mode pair.
##   Q2 prediction with a 90% predictive interval for held-out graphs (needs only their cores).
##   Q3 vertex-level query: posterior of B = U Theta U' -> node strength with 90% credible band.
## References that need the raw features again: lasso refitted on the two-sided cores at the same
## number of nonzeros, and edge lasso (time only).
## Outputs: results/icde/projpost_accuracy.csv, projpost_queries.csv
## =====================================================================================
source("Codes/icde/icde_common.R")
smoke <- has_flag("--smoke"); n_splits <- if (smoke) 1 else 10; M <- 500
o1 <- out_path("projpost_accuracy.csv", smoke); o2 <- out_path("projpost_queries.csv", smoke)
for (f in c(o1, o2)) if (file.exists(f)) file.remove(f)
soft <- function(z, t) sign(z) * pmax(abs(z) - t, 0)
proj_l1 <- function(theta, G, lam, b = theta, iters = 200) {     # coordinate descent on the cached quadratic
  r <- as.numeric(G %*% theta)
  for (it in seq_len(iters)) { b0 <- b
    for (j in seq_along(b)) b[j] <- soft(r[j] - sum(G[j, -j] * b[-j]), lam / 2) / G[j, j]
    if (max(abs(b - b0)) < 1e-7) break }
  b
}
lam_for_k <- function(theta, G, k) {                             # smallest lambda leaving at most k nonzeros
  lo <- 0; hi <- 2 * max(abs(G %*% theta))
  for (i in 1:40) { mid <- (lo + hi) / 2; if (sum(proj_l1(theta, G, mid) != 0) > k) lo <- mid else hi <- mid }
  hi
}
for (atlas in c("sch100", "power229")) {
  d <- load_set(atlas, "age_continuous"); conn <- d$conn; y <- d$y; subj <- d$subj; V <- d$V
  Xe <- build_edge_features(conn)$X_edge
  for (s in seq_len(n_splits)) {
    sp <- make_split(subj, s); tr <- sp$tr; te <- sp$te; n <- length(tr)
    Q <- template_Q(conn, tr, 10); XT <- feat_twosided(conn, Q)
    ## ---------- fit once ----------
    t_fit <- system.time({
      pp <- standardize_train_test(XT[tr, ], XT[te, ], y[tr]); X <- pp$X; yc <- pp$y; q <- ncol(X)
      a <- ridge_dual_cv_folds(XT[tr, ], y[tr], XT[te, ], subject_folds(subj[tr]))$alpha_sel
      A <- solve(crossprod(X) + a * diag(q)); m <- as.numeric(A %*% crossprod(X, yc))
      dfr <- sum(diag(X %*% A %*% t(X))); s2 <- sum((yc - X %*% m)^2) / (n - dfr - 1); S <- s2 * A
      G <- crossprod(X) / n })[3]
    cache_bytes <- 8 * (length(m) + length(S) + length(G) + length(Q) + 2 * q + 2)
    Xt <- pp$Xt; yte <- y[te]; ybar <- pp$ybar                   # held-out cores only
    set.seed(8000 + s); L <- chol(S); draws <- matrix(rnorm(M * q), M, q) %*% L + rep(m, each = M)
    acc <- list(); addacc <- function(method, k, b, secs) acc[[length(acc) + 1]] <<- data.frame(atlas = atlas, split = s,
      method = method, k = k, nnz = sum(b != 0), spearman = spearman(yte, as.numeric(Xt %*% b)), rmse = rmse(yte, as.numeric(Xt %*% b) + ybar), seconds = secs)
    addacc("dense_posterior_mean", q, m, t_fit)
    ## ---------- Q1: sparse summaries from the cache, and the refit that needs the features ----------
    cv <- glmnet::glmnet(X, yc, alpha = 1, standardize = FALSE, intercept = FALSE, nlambda = 200)
    pip <- list()
    for (k in c(5, 10, 20)) {
      t1 <- system.time({ lam <- lam_for_k(m, G, k); b <- proj_l1(m, G, lam) })[3]
      addacc("projection_of_mean", k, b, t1)
      t2 <- system.time(Bd <- t(apply(draws, 1, proj_l1, G = G, lam = lam, iters = 60)))[3]
      pip[[as.character(k)]] <- colMeans(Bd != 0)
      addacc("projection_draw_mean", k, colMeans(Bd), t2)
      j <- max(which(cv$df <= k)); addacc("lasso_refit_on_cores", k, as.numeric(cv$beta[, j]), NA)
    }
    append_csv(do.call(rbind, acc), o1)
    ## ---------- Q2: predictive intervals; Q3: vertex-level map with credible band ----------
    t3 <- system.time({ mu <- as.numeric(Xt %*% m) + ybar; sdp <- sqrt(s2 + rowSums((Xt %*% S) * Xt))
      cover <- mean(abs(yte - mu) <= qnorm(0.95) * sdp) })[3]
    t4 <- system.time({                                           # coefficients back on the original feature scale
      sc <- apply(XT[tr, ], 2, sd); sc[sc == 0] <- 1
      idx <- which(upper.tri(diag(10), diag = TRUE), arr.ind = TRUE)
      node <- apply(draws[seq_len(200), ], 1, function(th) { Th <- matrix(0, 10, 10); v <- th / sc
        Th[idx] <- v; Th <- Th + t(Th) - diag(diag(Th)); B <- Q %*% Th %*% t(Q); diag(B) <- 0; sqrt(rowSums(B^2)) })
      lo <- apply(node, 1, quantile, 0.05); hi <- apply(node, 1, quantile, 0.95); md <- rowMeans(node)
      top <- order(md, decreasing = TRUE)[seq_len(10)] })[3]
    t_edge <- system.time(glmnet::cv.glmnet(Xe[tr, ], y[tr], alpha = 1, nfolds = 5))[3]
    append_csv(data.frame(atlas = atlas, split = s, n_train = n, V = V, cache_bytes = cache_bytes, edge_store_bytes = 8 * n * ncol(Xe),
      fit_seconds = t_fit, pi90_coverage = cover, pi90_mean_width = mean(2 * qnorm(0.95) * sdp), predict_seconds = t3,
      node_query_seconds = t4, node_band_rel_width = mean((hi - lo) / md), top10_nodes = paste(sort(top), collapse = " "),
      pip10_above_half = sum(pip[["10"]] > 0.5), pip10_max = max(pip[["10"]]), pip10_idx = paste(which(pip[["10"]] > 0.5), collapse = " "),
      edge_lasso_refit_seconds = t_edge), o2)
    cat(atlas, "split", s, "done\n"); flush.console()
  }
}
