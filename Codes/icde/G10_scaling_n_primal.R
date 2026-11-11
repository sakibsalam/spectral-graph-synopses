## =====================================================================================
## G10: primal-form ridge (q x q eigendecomposition per fold) versus n at V = 229, to sit
## beside the dual-form timings of G3. Same data construction and seeds as G3.
##   Rscript Codes/icde/G10_scaling_n_primal.R [--smoke] [--nmax=3200]
## Output: results/icde/scaling_n_primal.csv
## =====================================================================================
source("Codes/icde/icde_common.R")
smoke <- has_flag("--smoke")
nmax <- as.numeric(arg_val("--nmax", 3200))
Ns <- if (smoke) c(100, 200) else c(200, 400, 800, 1600, 3200)
Ns <- Ns[Ns <= nmax]
K <- 10; R <- 10; n_te <- 100; SLOW_S <- 120
out <- out_path("scaling_n_primal.csv", smoke)

## Primal CV ridge: identical grid, folds, selection and refit as raw_edge_ridge_cv.
ridge_primal_cv <- function(X, y, Xt, nfolds = 5, n_alpha = 25, seed = 1) {
  pp <- standardize_train_test(X, Xt, y); Xs <- pp$X; yc <- pp$y; n <- nrow(Xs)
  grid <- (sum(diag(crossprod(Xs))) / n) * 10^seq(-6, 6, length.out = n_alpha)  # = mean eig of XX'
  set.seed(seed); fold_id <- sample(rep(seq_len(nfolds), length.out = n))
  err <- matrix(NA_real_, nfolds, n_alpha)
  for (k in seq_len(nfolds)) {
    i <- which(fold_id != k); o <- which(fold_id == k)
    Xf <- Xs[i, , drop = FALSE]; eg <- eigen(crossprod(Xf), symmetric = TRUE)
    Vm <- eg$vectors; d <- eg$values; Vb <- as.numeric(crossprod(Vm, crossprod(Xf, yc[i])))
    Xh <- Xs[o, , drop = FALSE] %*% Vm                      # held-out scores in eigenbasis
    for (a in seq_len(n_alpha)) {
      bt <- Vb / (d + grid[a])                              # beta in eigenbasis
      err[k, a] <- mean((yc[o] - as.numeric(Xh %*% bt))^2)
    }
  }
  a <- which.min(colMeans(err))
  eg <- eigen(crossprod(Xs), symmetric = TRUE)
  b <- as.numeric(eg$vectors %*% (as.numeric(crossprod(eg$vectors, crossprod(Xs, yc))) / (eg$values + grid[a])))
  list(pred = as.numeric(pp$Xt %*% b + pp$ybar), alpha_sel = grid[a])
}

base <- load_set("power229", "nbk_0v1"); V <- base$V; p <- V * (V - 1) / 2
build_data <- function(n) {
  set.seed(4000 + n); idx <- sample(seq_along(base$conn), n + n_te, replace = TRUE)
  conn <- lapply(idx, function(i) { E <- matrix(rnorm(V * V, sd = 0.05), V, V); E <- (E + t(E)) / 2
                                    A <- base$conn[[i]] + E; diag(A) <- 0; A })
  list(conn = conn, y = base$y[idx])
}

## ---- agreement check (untimed, n = 200, two-sided features) ----
{
  d <- build_data(200); tr <- 1:200; te <- 200 + seq_len(n_te)
  Q <- template_Q(d$conn, tr, K, "adjacency"); Xs <- feat_twosided(d$conn, Q)
  a <- ridge_primal_cv(Xs[tr, ], d$y[tr], Xs[te, ]); b <- raw_edge_ridge_cv(Xs[tr, ], d$y[tr], Xs[te, ])
  bp <- if (!is.null(b$pred)) b$pred else b$prob
  cat(sprintf("AGREEMENT max|pred_primal - pred_dual| = %.3e  (names in dual result: %s)\n",
              max(abs(a$pred - as.numeric(bp))), paste(names(b), collapse = ",")))
  if (!is.null(b$alpha_sel)) cat(sprintf("alpha_sel primal %.6g dual %.6g\n", a$alpha_sel, b$alpha_sel))
  rm(d, Q, Xs); invisible(gc())
}

rows <- list(); wall0 <- proc.time()[3]
timed <- function(n, design, stage, fn, reps, X) {
  for (r in seq_len(reps)) {
    invisible(gc(reset = TRUE)); b0 <- sum(gc()[, 2]); t0 <- proc.time()[3]
    fn(); secs <- proc.time()[3] - t0; peak <- sum(gc()[, 6]) - b0
    rows[[length(rows) + 1]] <<- data.frame(V = V, p = p, n_train = n, design = design, stage = stage,
      repeat_id = r, seconds = secs, peak_mb = peak, n_features = ncol(X),
      train_matrix_mb = 8 * n * ncol(X) / 2^20)
    cat(sprintf("  n=%5d %-10s %-20s rep%d %9.2fs\n", n, design, stage, r, secs)); flush.console()
    if (secs > SLOW_S) break
  }
}

for (n in Ns) {
  cat(sprintf("\n===== n = %d =====\n", n))
  d <- build_data(n); y <- d$y; tr <- seq_len(n); te <- n + seq_len(n_te)
  reps <- if (smoke) 1 else 3
  Q <- template_Q(d$conn, tr, K, "adjacency")
  Xs <- feat_twosided(d$conn, Q); Xo <- feat_onesided(d$conn, Q[, seq_len(R), drop = FALSE])
  rm(d); invisible(gc())
  Xs_tr <- Xs[tr, , drop = FALSE]; Xs_te <- Xs[te, , drop = FALSE]
  Xo_tr <- Xo[tr, , drop = FALSE]; Xo_te <- Xo[te, , drop = FALSE]
  timed(n, "two_sided", "fit_ridge_primal_cv", function() ridge_primal_cv(Xs_tr, y[tr], Xs_te), reps, Xs_tr)
  timed(n, "one_sided", "fit_ridge_primal_cv", function() ridge_primal_cv(Xo_tr, y[tr], Xo_te), reps, Xo_tr)
  rm(Xs, Xo, Xs_tr, Xs_te, Xo_tr, Xo_te); invisible(gc())
  res <- do.call(rbind, rows); names(res)[names(res) == "repeat_id"] <- "repeat"
  write.csv(res, out, row.names = FALSE)
}
cat(sprintf("\nTotal wall time (agreement check excluded from none): %.1f s\n", proc.time()[3] - wall0))
cat("\nWrote", out, "\n")
