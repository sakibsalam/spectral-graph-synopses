## =====================================================================================
## E1: measured cost vs number of vertices V.   Run from the repo root, machine IDLE.
##   Rscript Codes/icde/E1_scaling.R [--smoke] [--ram_gb=32]
## Real Power-229 nbk_0v1 graphs on random node subsets (V = 50..229), synthetic graphs for
## V = 300, 400. Per V and design we time: template build + eigendecomposition, synopsis
## build, and each learner's fit. Iterative learners (SSVS, edge ARD, BMGC) are timed for a
## few iterations only; seconds_per_iter is measured, full-fit times are NOT projected here.
## Cost model reminder: edge ARD solves a p x p system per iteration (p^3, p = V(V-1)/2);
## SSVS solves n x n per iteration (n^2 p); two-sided ARD solves q x q (q = K(K+1)/2 = 55).
## =====================================================================================
source("Codes/icde/icde_common.R")
smoke <- has_flag("--smoke")
CAP   <- 600                                        # seconds per cell
K <- 10; R <- 10
Vs_real <- if (smoke) c(30, 60) else c(50, 100, 150, 200, 229)
Vs_syn  <- if (smoke) integer(0) else c(300, 400)
reps    <- if (smoke) 1 else 3
SSVS_IT <- if (smoke) 3 else 5
ARD_IT  <- 2
BMGC_IT <- if (smoke) 2 else 3
SLOW_S  <- 120                                      # after a cell this slow, do not repeat it

## RAM (GB): CLI arg, else sysctl (macOS) or /proc/meminfo (Linux)
ram_gb <- suppressWarnings(as.numeric(arg_val("--ram_gb", NA)))
if (is.na(ram_gb)) ram_gb <- tryCatch({
  if (Sys.info()[["sysname"]] == "Darwin") as.numeric(system("sysctl -n hw.memsize", intern = TRUE)) / 2^30
  else as.numeric(sub(".*:\\s*([0-9]+) kB.*", "\\1", grep("MemTotal", readLines("/proc/meminfo"), value = TRUE))) / 2^20
}, error = function(e) NA_real_)
BLAS_THREADS <- paste(sapply(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS"),
                              function(v) paste0(v, "=", Sys.getenv(v, "unset"))), collapse = ";")
cat(sprintf("RAM = %.1f GB; per-cell cap = %d s; smoke = %s\nthreads: %s\n", ram_gb, CAP, smoke, BLAS_THREADS))
cat("BLAS:", if (!is.null(sessionInfo()$BLAS)) sessionInfo()$BLAS else "unknown", "\n")

## BMGC (authors' sampler); needs gdata, pROC, mvtnorm. L is hardcoded to 2 inside guha_bmgc.R,
## harmless for one layer (both layers receive the same matrix).
bmgc_ok <- tryCatch({ suppressMessages(source("Codes/benchmarks/guha_bmgc_benchmark.R")); TRUE },
                    error = function(e) { cat("BMGC unavailable:", conditionMessage(e), "\n"); FALSE })

## ---- one timed cell: peak R-heap increment (MB), wall seconds, status, return value -----
run_cell <- function(fn, cap = CAP) {
  base <- sum(gc(reset = TRUE)[, 2])
  t0 <- proc.time()[3]; status <- "ok"; msg <- NA_character_
  setTimeLimit(elapsed = cap, transient = FALSE)
  val <- tryCatch(fn(), error = function(e) {
    msg <<- conditionMessage(e)
    status <<- if (grepl("time limit", msg)) "timeout"
               else if (grepl("cannot allocate|memory", msg, ignore.case = TRUE)) "infeasible_memory"
               else "error"
    NULL })
  setTimeLimit(elapsed = Inf, cpu = Inf, transient = FALSE)
  secs <- proc.time()[3] - t0
  list(status = status, seconds = secs, peak_mb = sum(gc()[, 6]) - base, msg = msg, val = val)
}

## ---- data ------------------------------------------------------------------------------
power <- load_set("power229", "nbk_0v1")
gen_synth <- function(V, N = 272, seed = 11) {      # block template + low-rank signal + noise
  set.seed(seed)
  W0 <- make_template_graph(V = V)$W
  Q <- make_basis(W0, K, "adjacency")$Q
  y <- rep(0:1, length.out = N); subj <- rep(seq_len(N / 2), each = 2)
  shift <- matrix(0, K, K); shift[1:3, 1:3] <- 0.25
  conn <- lapply(seq_len(N), function(i) {
    Z <- matrix(rnorm(K * K, sd = 0.6), K, K); Z <- (Z + t(Z)) / 2 + y[i] * shift
    E <- matrix(rnorm(V * V, sd = 0.35), V, V); E <- (E + t(E)) / 2
    A <- W0 + Q %*% Z %*% t(Q) + E; diag(A) <- 0; A })
  list(conn = conn, y = y, subj = subj)
}
get_data <- function(V) {
  if (V <= 229) {
    sel <- if (V == 229) seq_len(229) else { set.seed(1000 + V); sort(sample.int(229, V)) }
    list(conn = lapply(power$conn, function(A) A[sel, sel]), y = power$y, subj = power$subj, source = "power229_subset")
  } else c(gen_synth(V), list(source = "synthetic"))
}

rows <- list(); ard_hist <- NULL                   # ard_hist: last measured (p, sec/iter) for the p^3 guard
slow <- character(0)
add_row <- function(V, p, design, learner, rep, cell, n_iter = NA, bytes = NA, note = "", spi = NA_real_) {
  ok <- identical(cell$status, "ok")
  rows[[length(rows) + 1]] <<- data.frame(
    V = V, p = p, design = design, learner = learner, repeat_id = rep, seconds = if (is.na(cell$seconds)) NA else cell$seconds,
    seconds_per_iter = if (ok && !is.na(spi)) spi else if (ok && !is.na(n_iter)) cell$seconds / n_iter else NA_real_,
    n_iter_timed = n_iter, peak_mb = cell$peak_mb, bytes_per_graph = bytes, status = cell$status,
    measured = ok, blas_threads = BLAS_THREADS, note = if (is.na(cell$msg)) note else paste(note, substr(cell$msg, 1, 120)),
    stringsAsFactors = FALSE)
  cat(sprintf("  V=%3d %-10s %-22s rep%d %-9s %8.2fs\n", V, design, learner, rep, cell$status, cell$seconds)); flush.console()
}
skipped <- function(status, note) list(status = status, seconds = NA_real_, peak_mb = NA_real_, msg = note, val = NULL)

## untimed warmup: first calls of glmnet / ARD / eigen pay one-off loading costs
invisible(tryCatch({ set.seed(1); Xw <- matrix(rnorm(60 * 30), 60); yw <- rep(0:1, 30)
  graph_spectral_ard_logit(Xw[1:40, ], yw[1:40], Xw[41:60, ]); raw_edge_glmnet_logit(Xw[1:40, ], yw[1:40], Xw[41:60, ], alpha = 1)
  raw_edge_ridge_cv(Xw[1:40, ], yw[1:40], Xw[41:60, ]); eigen(crossprod(Xw)) }, error = function(e) NULL))

for (V in c(Vs_real, Vs_syn)) {
  cat(sprintf("\n===== V = %d =====\n", V))
  d <- tryCatch(get_data(V), error = function(e) { cat("data error:", conditionMessage(e), "\n"); NULL })
  if (is.null(d)) next
  conn <- d$conn; y <- d$y; sp <- make_split(d$subj, 1); tr <- sp$tr; te <- sp$te
  p <- V * (V - 1) / 2
  st <- new.env()
  for (r in seq_len(reps)) {
    ## --- template: mean |A| over training graphs + adjacency eigendecomposition (all V modes) ---
    c0 <- run_cell(function() { st$Q <- template_Q(conn, tr, K, "adjacency"); TRUE })
    add_row(V, p, "template", "mean_plus_eigen", r, c0)
    if (is.null(st$Q)) next
    ## --- synopsis builds (all N graphs); bytes per graph = object.size of one stored row ---
    c1 <- run_cell(function() { st$Xs <- feat_twosided(conn, st$Q); TRUE })
    add_row(V, p, "two_sided", "build", r, c1, bytes = 8 * ncol(st$Xs))
    ## one-sided reuses the template eigenvectors (R <= K), so both build times exclude the eigendecomposition
    c2 <- run_cell(function() { st$Xo <- feat_onesided(conn, st$Q[, seq_len(R), drop = FALSE]); TRUE })
    add_row(V, p, "one_sided", "build", r, c2, bytes = 8 * ncol(st$Xo))
    c3 <- run_cell(function() { st$Xe <- build_edge_features(conn)$X_edge; TRUE })
    add_row(V, p, "edge_store", "build", r, c3, bytes = 8 * ncol(st$Xe))
    if (is.null(st$Xs) || is.null(st$Xo) || is.null(st$Xe)) next
    b2 <- 8 * ncol(st$Xs); bo <- 8 * ncol(st$Xo); be <- 8 * ncol(st$Xe)   # float64 payload bytes per graph
    Xs <- st$Xs; Xo <- st$Xo; Xe <- st$Xe

    do <- function(design, learner, bytes, fn, n_iter_fun = function(v) NA, note = "") {
      key <- paste(V, design, learner)
      if (r > 1 && key %in% slow) return(invisible(NULL))
      cell <- run_cell(fn)
      if (cell$seconds > SLOW_S) slow <<- c(slow, key)
      add_row(V, p, design, learner, r, cell, n_iter = n_iter_fun(cell$val), bytes = bytes, note = note)
      invisible(cell)
    }
    ## --- fits ---
    do("two_sided", "ard_logit_K10", b2, function() graph_spectral_ard_logit(Xs[tr, ], y[tr], Xs[te, ]))
    do("one_sided", "ridge_dual_cv_R10", bo, function() raw_edge_ridge_cv(Xo[tr, ], y[tr], Xo[te, ]))
    do("edge_store", "ridge_dual_cv", be, function() raw_edge_ridge_cv(Xe[tr, ], y[tr], Xe[te, ]))
    do("edge_store", "lasso_glmnet", be, function() raw_edge_glmnet_logit(Xe[tr, ], y[tr], Xe[te, ], alpha = 1))
    do("edge_store", "ssvs_iter", be, function() spike_slab_ssvs(Xe[tr, ], y[tr], Xe[te, ], n_iter = SSVS_IT, burn = 1),
       n_iter_fun = function(v) SSVS_IT, note = sprintf("%d iterations incl. standardization", SSVS_IT))

    ## --- edge ARD: first 2 iterations. Guards (recorded, not run): ard_fit holds about six p x p
    ##     double matrices at once (crossprod, diag, sum, identity, inverse, LU copy), so 6*p^2*8 must fit
    ##     in 70% of RAM; and the p^3 projection of ARD_IT iterations from the previous measured V must stay
    ##     under the cap, because setTimeLimit cannot interrupt a running LAPACK solve.
    key <- paste(V, "edge_store", "ard_iter")
    if (!is.na(ram_gb) && 6 * p * p * 8 > 0.7 * ram_gb * 2^30) {
      if (r == 1) add_row(V, p, "edge_store", "ard_iter", r, skipped("infeasible_memory",
                          sprintf("6*p*p*8 = %.1f GB > 70%% of %.0f GB RAM", 6 * p * p * 8 / 2^30, ram_gb)))
    } else if (!is.null(ard_hist) && ard_hist$t * (p / ard_hist$p)^3 * ARD_IT > CAP) {
      if (r == 1) add_row(V, p, "edge_store", "ard_iter", r, skipped("timeout",
                          sprintf("not run (PROJECTED, not measured): p^3 extrapolation from V=%d gives %.0f s/iter; %d iterations exceed the cap", ard_hist$V, ard_hist$t * (p / ard_hist$p)^3, ARD_IT)))
    } else if (!(r > 1 && key %in% slow)) {
      cell <- run_cell(function() {
        pp <- standardize_train_test(Xe[tr, , drop = FALSE], Xe[te, , drop = FALSE], y[tr])
        ard_fit(pp$X, pp$y, max_iter = ARD_IT) })
      nit <- if (is.null(cell$val)) NA else cell$val$iter
      if (cell$seconds > SLOW_S) slow <- c(slow, key)
      add_row(V, p, "edge_store", "ard_iter", r, cell, n_iter = nit, bytes = be, note = "first iterations (all p features active)")
      if (cell$status == "ok" && !is.na(nit)) ard_hist <- list(V = V, p = p, t = cell$seconds / nit)
    }

    ## --- BMGC: authors' Gibbs sampler, per iteration via bmgc_scaling_probe(); L fixed at 2 in source ---
    if (bmgc_ok && !(r > 1 && paste(V, "edge_store", "bmgc_iter") %in% slow)) {
      cell <- run_cell(function() { invisible(capture.output(per <- bmgc_scaling_probe(conn, y, V_sub = NULL, niter = BMGC_IT))); per })
      if (cell$seconds > SLOW_S) slow <- c(slow, paste(V, "edge_store", "bmgc_iter"))
      add_row(V, p, "edge_store", "bmgc_iter", r, cell, n_iter = BMGC_IT, bytes = be,
              spi = if (is.numeric(cell$val)) cell$val else NA_real_,
              note = "seconds = probe incl. setup; seconds_per_iter = probe's own estimate (incl. setup / niter); L hardcoded 2 in guha_bmgc.R")
    }
    rm(list = intersect(c("Xs", "Xo", "Xe"), ls())); st <- new.env(); invisible(gc())
  }
  res <- do.call(rbind, rows); names(res)[names(res) == "repeat_id"] <- "repeat"
  write.csv(res, out_path("scaling.csv", smoke), row.names = FALSE)   # rewritten after every V
}
cat("\nWrote", out_path("scaling.csv", smoke), "\n")
