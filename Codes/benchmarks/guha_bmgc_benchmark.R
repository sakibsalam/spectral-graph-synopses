## =====================================================================================
## Wrapper to run the AUTHORS' REAL Guha code (Bayesian Multiplex Graph Classifier, BMGC) as
## a benchmark, replacing our approximate proxy in Codes/benchmarks/competitor_benchmarks.R.
##
## Source: github.com/sguhastat/Bayesian-Multiplex-Graph-Classifier-for-Functional-Connectivity-
## Across-Cognitive-Control-Tasks (Rodriguez-Acosta, Guha & Dinov, Neuroinformatics 2024). The
## exact function is in Codes/benchmarks/guha_bmgc.R (a copy of their `binary_multi_network_function_with_
## results.R` with TWO minimal patches, both marked `# PATCH`: accept an external test split so
## the comparison is subject-stratified & comparable; and return the test predictions).
##
## Requires: BayesLogit, GIGrvg, pscl, MCMCpack, LaplacesDemon, mvtnorm, gdata, MASS, pROC
## (installed). BMGC is multiplex (L=2 layers); for a single connectome we pass it as both layers.
##
## !!! SCALABILITY (measured on this data) !!!
##   BMGC's per-node Gibbs step does an O(N^3) dmvnorm inside a loop over all V nodes, so cost is
##   ~ niter * V * N^3. On the OSF data (V=229, p=26,106, N~276) this is **~20 s PER ITERATION**
##   => a converged run (niter=2000-10000) is ~11-57 HOURS for a SINGLE fit; a 10-split benchmark
##   is weeks of compute. Reducing V to make it tractable destroys the (distributed) signal
##   (all methods fall to ~chance at V=40). This is precisely the edge-dimension bottleneck our
##   graph-spectral reduction removes (our fit: milliseconds in the K=10 spectral basis).
##   => In the headline V=229 benchmark table, "Guha" is the fast global-local-shrinkage proxy
##      (Codes/benchmarks/competitor_benchmarks.R); this exact-code path is provided for cluster runs / smaller V.
## =====================================================================================

suppressMessages({
  library(gdata)
  source("Codes/benchmarks/guha_bmgc.R")   # the authors' patched Gibbs sampler
})

## Run the real BMGC for a single train/test split of single-layer connectomes.
##   conn: list of V x V matrices (all subjects);  y: 0/1 outcome;  test_idx: indices of test rows.
guha_bmgc_predict <- function(conn, y, test_idx, R = 6, niter = 2000, n_burnin = 1000, seed = 11) {
  Xmat <- t(sapply(conn, function(A) upperTriangle(A, diag = FALSE, byrow = TRUE)))
  t0 <- proc.time()[3]
  m <- binary_multilayer_network_model(y = y, Xmat1 = Xmat, Xmat2 = Xmat, test_idxs = test_idx,
                                       b1.tr = NULL, b2.tr = NULL, R = R,
                                       niter = niter, n_burnin = n_burnin, seed = seed)
  list(prob = m$test_pred, y_test = m$test_y, auc = as.numeric(m$auc),
       node_probabilities = m$node_probabilities, seconds = proc.time()[3] - t0)
}

## Empirical per-iteration cost (to document the scaling barrier without a full run).
bmgc_scaling_probe <- function(conn, y, V_sub = NULL, niter = 3) {
  if (!is.null(V_sub)) {
    av <- apply(t(sapply(conn, function(A) rowSums(abs(A)))), 2, var)
    sel <- order(av, decreasing = TRUE)[seq_len(V_sub)]
    conn <- lapply(conn, function(A) A[sel, sel])
  }
  V <- nrow(conn[[1]]); N <- length(conn)
  te <- sample(N, round(0.3 * N))
  t0 <- proc.time()[3]
  invisible(guha_bmgc_predict(conn, y, te, niter = niter, n_burnin = 0))
  per <- (proc.time()[3] - t0) / niter
  cat(sprintf("BMGC @ V=%d N=%d: %.1f s/iter -> niter=2000 ~ %.1f h, niter=10000 ~ %.1f h\n",
              V, N, per, per * 2000 / 3600, per * 10000 / 3600))
  invisible(per)
}
