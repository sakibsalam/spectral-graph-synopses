## =====================================================================================
## G3: measured cost versus the number of training graphs n at V = 229. Machine IDLE.
##   Rscript Codes/icde/G3_scaling_n.R [--smoke] [--nmax=3200]
## The cohort has 288 graphs for the binary task, so larger collections are simulated by
## resampling those graphs with replacement and adding symmetric Gaussian noise (sd 0.05) to
## every edge; labels follow the resampled graphs. Only cost is read from this experiment,
## never accuracy. Learners and timing conventions are those of E1 (one call = standardize,
## fit with internal penalty selection, score 100 held-out graphs). Cells slower than 120 s
## are measured once; the others three times (report the median).
## Output: results/icde/scaling_n.csv (+ scaling_n_env.txt with the platform)
## =====================================================================================
source("Codes/icde/icde_common.R")
smoke <- has_flag("--smoke")
nmax <- as.numeric(arg_val("--nmax", 3200))
Ns <- if (smoke) c(100, 200) else c(200, 400, 800, 1600, 3200)
Ns <- Ns[Ns <= nmax]
K <- 10; R <- 10; n_te <- 100; SLOW_S <- 120
out <- out_path("scaling_n.csv", smoke)
writeLines(c(capture.output(sessionInfo()),
             paste("cpu:", tryCatch(system("sysctl -n machdep.cpu.brand_string", intern = TRUE), error = function(e) Sys.info()[["machine"]])),
             paste("ram_gb:", tryCatch(as.numeric(system("sysctl -n hw.memsize", intern = TRUE)) / 2^30, error = function(e) NA)),
             paste("threads:", paste(sapply(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"),
                                            function(v) paste0(v, "=", Sys.getenv(v, "unset"))), collapse = ";"))),
           out_path("scaling_n_env.txt", smoke))

base <- load_set("power229", "nbk_0v1"); V <- base$V; p <- V * (V - 1) / 2
invisible(tryCatch({ set.seed(1); Xw <- matrix(rnorm(60 * 30), 60); yw <- rep(0:1, 30)   # untimed warmup
  graph_spectral_ard_logit(Xw[1:40, ], yw[1:40], Xw[41:60, ]); raw_edge_ridge_cv(Xw[1:40, ], yw[1:40], Xw[41:60, ]) },
  error = function(e) NULL))

rows <- list()
timed <- function(n, design, stage, fn, reps) {
  for (r in seq_len(reps)) {
    invisible(gc(reset = TRUE)); b0 <- sum(gc()[, 2]); t0 <- proc.time()[3]
    val <- fn(); secs <- proc.time()[3] - t0; peak <- sum(gc()[, 6]) - b0
    rows[[length(rows) + 1]] <<- data.frame(V = V, p = p, n_train = n, design = design, stage = stage,
      repeat_id = r, seconds = secs, peak_mb = peak, n_features = if (is.matrix(val)) ncol(val) else NA,
      train_matrix_mb = if (is.matrix(val)) 8 * n * ncol(val) / 2^20 else NA)
    cat(sprintf("  n=%5d %-10s %-18s rep%d %9.2fs\n", n, design, stage, r, secs)); flush.console()
    if (secs > SLOW_S) break
  }
  val
}

for (n in Ns) {
  cat(sprintf("\n===== n = %d =====\n", n))
  set.seed(4000 + n); idx <- sample(seq_along(base$conn), n + n_te, replace = TRUE)
  conn <- lapply(idx, function(i) { E <- matrix(rnorm(V * V, sd = 0.05), V, V); E <- (E + t(E)) / 2
                                    A <- base$conn[[i]] + E; diag(A) <- 0; A })
  y <- base$y[idx]; tr <- seq_len(n); te <- n + seq_len(n_te)
  reps <- if (smoke) 1 else 3
  Q  <- timed(n, "template", "mean_plus_eigen", function() template_Q(conn, tr, K, "adjacency"), reps)
  Xs <- timed(n, "two_sided", "build", function() feat_twosided(conn, Q), reps)
  Xo <- timed(n, "one_sided", "build", function() feat_onesided(conn, Q[, seq_len(R), drop = FALSE]), reps)
  Xe <- timed(n, "edge_store", "build", function() build_edge_features(conn)$X_edge, reps)
  rm(conn); invisible(gc())
  timed(n, "two_sided", "fit_ard_logit", function() { graph_spectral_ard_logit(Xs[tr, ], y[tr], Xs[te, ]); Xs[tr, ] }, reps)
  timed(n, "two_sided", "fit_ridge_dual_cv", function() { raw_edge_ridge_cv(Xs[tr, ], y[tr], Xs[te, ]); Xs[tr, ] }, reps)
  timed(n, "one_sided", "fit_ridge_dual_cv", function() { raw_edge_ridge_cv(Xo[tr, ], y[tr], Xo[te, ]); Xo[tr, ] }, reps)
  timed(n, "edge_store", "fit_ridge_dual_cv", function() { raw_edge_ridge_cv(Xe[tr, ], y[tr], Xe[te, ]); Xe[tr, ] }, reps)
  rm(Xs, Xo, Xe); invisible(gc())
  res <- do.call(rbind, rows); names(res)[names(res) == "repeat_id"] <- "repeat"
  write.csv(res, out, row.names = FALSE)
}
cat("\nWrote", out, "\n")
