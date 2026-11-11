## Graph-adapted Bayesian shrinkage and projection for network-valued predictors
## Base R implementation, designed to mirror the Python prototype.
## It simulates scalar-on-network regression where the true coefficient network
## is sparse in a Laplacian spectral basis.

soft_threshold <- function(z, lambda) {
  sign(z) * pmax(abs(z) - lambda, 0)
}

logspace <- function(a, b, n) {
  10 ^ seq(a, b, length.out = n)
}

rmse <- function(y, yhat) {
  sqrt(mean((y - yhat)^2))
}

make_template_graph <- function(V = 40, n_blocks = 4, within = 0.7, between = 0.12) {
  blocks <- rep(seq_len(n_blocks), each = ceiling(V / n_blocks))[seq_len(V)]
  W <- matrix(0, V, V)
  for (i in seq_len(V - 1)) {
    for (j in (i + 1):V) {
      val <- if (blocks[i] == blocks[j]) within else between
      W[i, j] <- val
      W[j, i] <- val
    }
  }
  diag(W) <- 0
  list(W = W, blocks = blocks)
}

laplacian_basis <- function(W) {
  d <- rowSums(W)
  L <- diag(d) - W
  eg <- eigen(L, symmetric = TRUE)
  ord <- order(eg$values)
  list(values = eg$values[ord], vectors = eg$vectors[, ord, drop = FALSE])
}

upper_edge_index <- function(V) {
  which(upper.tri(matrix(0, V, V)), arr.ind = TRUE)
}

vech_features <- function(C) {
  K <- nrow(C)
  idx <- which(upper.tri(matrix(0, K, K), diag = TRUE), arr.ind = TRUE)
  f <- C[idx]
  off_diag <- idx[, 1] != idx[, 2]
  f[off_diag] <- 2 * f[off_diag]
  as.numeric(f)
}

simulate_spectral_network_regression <- function(n = 120,
                                                 n_test = 300,
                                                 V = 40,
                                                 K = 10,
                                                 n_active = 6,
                                                 sigma = 1.0,
                                                 seed = 1) {
  set.seed(seed)
  tg <- make_template_graph(V = V)
  W0 <- tg$W
  lb <- laplacian_basis(W0)
  Q <- lb$vectors[, seq_len(K), drop = FALSE]

  M_true <- matrix(0, K, K)
  candidate_pairs <- do.call(
    rbind,
    lapply(seq_len(K), function(a) {
      out <- NULL
      for (b in a:K) {
        if (a <= 6 && b <= 6 && abs(a - b) <= 2) {
          out <- rbind(out, c(a, b))
        }
      }
      out
    })
  )

  chosen <- sample(seq_len(nrow(candidate_pairs)), n_active, replace = FALSE)
  active_pairs <- candidate_pairs[chosen, , drop = FALSE]

  for (r in seq_len(nrow(active_pairs))) {
    a <- active_pairs[r, 1]
    b <- active_pairs[r, 2]
    val <- sample(c(-1, 1), size = 1) * runif(1, 0.8, 1.2)
    M_true[a, b] <- val
    M_true[b, a] <- val
  }

  B_true <- Q %*% M_true %*% t(Q)
  diag(B_true) <- 0

  edge_idx <- upper_edge_index(V)
  p <- nrow(edge_idx)
  q <- K * (K + 1) / 2

  draw_networks <- function(N) {
    X_edge <- matrix(0, N, p)
    X_spec <- matrix(0, N, q)

    for (i in seq_len(N)) {
      Z <- matrix(rnorm(K * K, sd = 0.6), K, K)
      Z <- (Z + t(Z)) / 2

      E <- matrix(rnorm(V * V, sd = 0.35), V, V)
      A <- W0 + Q %*% Z %*% t(Q) + E
      A <- (A + t(A)) / 2
      diag(A) <- 0

      X_edge[i, ] <- A[edge_idx]
      C <- t(Q) %*% A %*% Q
      X_spec[i, ] <- vech_features(C)
    }

    list(X_edge = X_edge, X_spec = X_spec)
  }

  train <- draw_networks(n)
  test <- draw_networks(n_test)

  beta_edge <- 2 * B_true[edge_idx]
  eta <- as.numeric(train$X_edge %*% beta_edge)
  scale_signal <- 2.0 / sd(eta)

  B_true <- B_true * scale_signal
  M_true <- M_true * scale_signal
  beta_edge <- 2 * B_true[edge_idx]

  y <- as.numeric(train$X_edge %*% beta_edge + rnorm(n, sd = sigma))
  y_test <- as.numeric(test$X_edge %*% beta_edge + rnorm(n_test, sd = sigma))

  idxK <- which(upper.tri(matrix(0, K, K), diag = TRUE), arr.ind = TRUE)
  theta_true <- M_true[idxK]

  list(
    X_edge = train$X_edge,
    X_spec = train$X_spec,
    y = y,
    Xt_edge = test$X_edge,
    Xt_spec = test$X_spec,
    y_test = y_test,
    beta_edge = beta_edge,
    theta_true = as.numeric(theta_true),
    B_true = B_true,
    Q = Q,
    active_pairs = active_pairs,
    evals = lb$values[seq_len(K)]
  )
}

standardize_train_test <- function(X, Xt, y) {
  x_center <- colMeans(X)
  x_scale <- apply(X, 2, sd)
  x_scale[x_scale == 0 | is.na(x_scale)] <- 1

  Xs <- sweep(sweep(X, 2, x_center, "-"), 2, x_scale, "/")
  Xts <- sweep(sweep(Xt, 2, x_center, "-"), 2, x_scale, "/")

  ybar <- mean(y)
  yc <- y - ybar

  list(X = Xs, Xt = Xts, y = yc, ybar = ybar)
}

ridge_fit <- function(X, y, alpha = 50) {
  p <- ncol(X)
  as.numeric(solve(crossprod(X) + alpha * diag(p), crossprod(X, y)))
}

ard_fit <- function(X,
                    y,
                    max_iter = 500,
                    tol = 1e-4,
                    threshold_lambda = 1e4,
                    eps = 1e-8) {
  ## Empirical-Bayes ARD regression for y = X beta + error.
  ## beta_j | lambda_j ~ N(0, lambda_j^{-1}).
  ## This mirrors the MacKay/Tipping update used by ARD implementations.
  n <- nrow(X)
  p <- ncol(X)

  alpha_noise <- 1 / max(var(y), eps)
  lambda <- rep(1, p)
  coef <- rep(0, p)

  for (iter in seq_len(max_iter)) {
    old_coef <- coef
    active <- which(lambda < threshold_lambda)

    if (length(active) == 0) break

    Xa <- X[, active, drop = FALSE]
    la <- lambda[active]

    A <- alpha_noise * crossprod(Xa) + diag(la, length(la))
    Sigma <- tryCatch(solve(A), error = function(e) qr.solve(A))
    mu <- as.numeric(alpha_noise * Sigma %*% crossprod(Xa, y))

    gamma <- pmax(0, 1 - la * diag(Sigma))
    lambda_new <- gamma / pmax(mu^2, eps)
    lambda_new[!is.finite(lambda_new)] <- threshold_lambda * 10

    resid <- as.numeric(y - Xa %*% mu)
    rss <- sum(resid^2)
    alpha_new <- (n - sum(gamma)) / max(rss, eps)
    alpha_new <- min(max(alpha_new, eps), 1 / eps)

    coef <- rep(0, p)
    coef[active] <- mu
    lambda[active] <- lambda_new
    lambda[lambda > threshold_lambda] <- threshold_lambda * 10
    alpha_noise <- alpha_new

    denom <- max(1, sqrt(sum(old_coef^2)))
    if (sqrt(sum((coef - old_coef)^2)) / denom < tol) break
  }

  list(coef = coef, lambda = lambda, alpha_noise = alpha_noise, iter = iter)
}

lasso_cd <- function(X,
                     y,
                     lambda,
                     max_iter = 5000,
                     tol = 1e-6) {
  ## Coordinate descent for objective:
  ## 0.5 / n * ||y - X beta||^2 + lambda * ||beta||_1.
  n <- nrow(X)
  p <- ncol(X)
  beta <- rep(0, p)
  r <- y
  z <- colSums(X^2) / n
  z[z == 0] <- 1

  for (iter in seq_len(max_iter)) {
    beta_old <- beta
    for (j in seq_len(p)) {
      r <- r + X[, j] * beta[j]
      rho <- sum(X[, j] * r) / n
      beta[j] <- soft_threshold(rho, lambda) / z[j]
      r <- r - X[, j] * beta[j]
    }
    if (sqrt(sum((beta - beta_old)^2)) < tol * max(1, sqrt(sum(beta_old^2)))) break
  }
  beta
}

graph_spectral_ard <- function(X, y, Xt) {
  t0 <- proc.time()[3]
  pp <- standardize_train_test(X, Xt, y)
  fit <- ard_fit(pp$X, pp$y)
  pred <- as.numeric(pp$Xt %*% fit$coef + pp$ybar)
  list(coef = fit$coef, pred = pred, seconds = proc.time()[3] - t0, fit = fit)
}

graph_spectral_projection <- function(X, y, Xt) {
  ## Sparse projection of a dense ridge posterior mean in graph-spectral space.
  ## This is the fast plug-in projection-posterior analogue.
  t0 <- proc.time()[3]
  pp <- standardize_train_test(X, Xt, y)

  dense <- ridge_fit(pp$X, pp$y, alpha = 10)
  target <- as.numeric(pp$X %*% dense)

  best_bic <- Inf
  best_beta <- rep(0, ncol(pp$X))
  best_lambda <- NA_real_

  for (lambda in logspace(-3, 0.5, 70)) {
    b <- lasso_cd(pp$X, target, lambda = lambda)
    rss <- sum((target - pp$X %*% b)^2)
    df <- max(sum(abs(b) > 1e-10), 1)
    bic <- nrow(pp$X) * log(rss / nrow(pp$X) + 1e-12) + df * log(nrow(pp$X))
    if (bic < best_bic) {
      best_bic <- bic
      best_beta <- b
      best_lambda <- lambda
    }
  }

  pred <- as.numeric(pp$Xt %*% best_beta + pp$ybar)
  list(coef = best_beta, pred = pred, seconds = proc.time()[3] - t0,
       lambda = best_lambda, bic = best_bic)
}

raw_edge_ard <- function(X, y, Xt) {
  t0 <- proc.time()[3]
  pp <- standardize_train_test(X, Xt, y)
  fit <- ard_fit(pp$X, pp$y)
  pred <- as.numeric(pp$Xt %*% fit$coef + pp$ybar)
  list(coef = fit$coef, pred = pred, seconds = proc.time()[3] - t0, fit = fit)
}

raw_edge_ridge <- function(X, y, Xt) {
  t0 <- proc.time()[3]
  pp <- standardize_train_test(X, Xt, y)
  b <- ridge_fit(pp$X, pp$y, alpha = 50)
  pred <- as.numeric(pp$Xt %*% b + pp$ybar)
  list(coef = b, pred = pred, seconds = proc.time()[3] - t0)
}

spike_slab_ssvs <- function(X,
                            y,
                            Xt,
                            n_iter = 1200,
                            burn = 400,
                            tau0 = 0.03,
                            tau1 = 1.0,
                            pi_incl = 0.05,
                            seed = 7) {
  set.seed(seed)
  t0 <- proc.time()[3]
  pp <- standardize_train_test(X, Xt, y)
  Xs <- pp$X
  yc <- pp$y
  Xts <- pp$Xt

  n <- nrow(Xs)
  p <- ncol(Xs)
  gamma <- runif(p) < pi_incl
  beta_mean <- rep(0, p)
  pip <- rep(0, p)

  sigma <- 0.5 * sd(yc)
  sig2 <- sigma^2

  for (it in seq_len(n_iter)) {
    d <- ifelse(gamma, tau1^2, tau0^2)
    u <- rnorm(p, sd = sqrt(d))
    delta <- rnorm(n, sd = sigma)
    v <- as.numeric(Xs %*% u + delta)

    A <- (Xs * matrix(d, nrow = n, ncol = p, byrow = TRUE)) %*% t(Xs) + sig2 * diag(n)
    R <- chol(A)
    w <- backsolve(R, forwardsolve(t(R), yc - v))
    beta <- u + d * as.numeric(crossprod(Xs, w))

    log_slab <- log(pi_incl + 1e-16) - log(tau1) - 0.5 * (beta / tau1)^2
    log_spike <- log(1 - pi_incl + 1e-16) - log(tau0) - 0.5 * (beta / tau0)^2
    prob <- 1 / (1 + exp(-(log_slab - log_spike)))
    gamma <- runif(p) < prob

    if (it > burn) {
      beta_mean <- beta_mean + beta
      pip <- pip + as.numeric(gamma)
    }
  }

  beta_mean <- beta_mean / (n_iter - burn)
  pip <- pip / (n_iter - burn)
  pred <- as.numeric(Xts %*% beta_mean + pp$ybar)

  list(coef = beta_mean, pip = pip, pred = pred, seconds = proc.time()[3] - t0)
}

run_one <- function(seed = 1) {
  dat <- simulate_spectral_network_regression(seed = seed)
  y_test <- dat$y_test
  theta_true <- dat$theta_true
  active <- abs(theta_true) > 1e-8
  s <- sum(active)

  rows <- list()

  fit <- graph_spectral_ard(dat$X_spec, dat$y, dat$Xt_spec)
  top <- order(abs(fit$coef), decreasing = TRUE)[seq_len(s)]
  rows[[length(rows) + 1]] <- data.frame(
    seed = seed,
    method = "graph spectral ARD",
    dimension = ncol(dat$X_spec),
    rmse = rmse(y_test, fit$pred),
    seconds = fit$seconds,
    spectral_hits_at_s = sum(active[top])
  )

  fit <- graph_spectral_projection(dat$X_spec, dat$y, dat$Xt_spec)
  top <- order(abs(fit$coef), decreasing = TRUE)[seq_len(s)]
  rows[[length(rows) + 1]] <- data.frame(
    seed = seed,
    method = "graph spectral projection",
    dimension = ncol(dat$X_spec),
    rmse = rmse(y_test, fit$pred),
    seconds = fit$seconds,
    spectral_hits_at_s = sum(active[top])
  )

  fit <- raw_edge_ard(dat$X_edge, dat$y, dat$Xt_edge)
  rows[[length(rows) + 1]] <- data.frame(
    seed = seed,
    method = "raw-edge ARD",
    dimension = ncol(dat$X_edge),
    rmse = rmse(y_test, fit$pred),
    seconds = fit$seconds,
    spectral_hits_at_s = NA
  )

  fit <- raw_edge_ridge(dat$X_edge, dat$y, dat$Xt_edge)
  rows[[length(rows) + 1]] <- data.frame(
    seed = seed,
    method = "raw-edge ridge",
    dimension = ncol(dat$X_edge),
    rmse = rmse(y_test, fit$pred),
    seconds = fit$seconds,
    spectral_hits_at_s = NA
  )

  fit <- spike_slab_ssvs(dat$X_edge, dat$y, dat$Xt_edge, seed = 1000 + seed)
  rows[[length(rows) + 1]] <- data.frame(
    seed = seed,
    method = "raw-edge spike-slab SSVS",
    dimension = ncol(dat$X_edge),
    rmse = rmse(y_test, fit$pred),
    seconds = fit$seconds,
    spectral_hits_at_s = NA
  )

  list(results = do.call(rbind, rows), active_pairs = dat$active_pairs)
}

run_experiment <- function(seeds = c(1, 2, 3)) {
  out <- lapply(seeds, run_one)
  df <- do.call(rbind, lapply(out, `[[`, "results"))

  cat("\nPer-run results:\n")
  print(df)

  cat("\nMean results by method:\n")
  print(aggregate(cbind(dimension, rmse, seconds, spectral_hits_at_s) ~ method,
                  data = df,
                  FUN = function(x) mean(x, na.rm = TRUE)))

  invisible(df)
}

## Uncomment to run:
## results <- run_experiment(c(1, 2, 3))
