## =====================================================================================
## Graph-spectral compressibility of real connectomes.
##
## Label-free empirical test of the method's CORE PREMISE and of the paper's "Controlled
## truncation" assumption (E[<A_i, R_{0,K}>^2] <= r_K^2). The graph-spectral model only
## "sees" each connectome A_i through the top-K x K block C_i = Q_K^T A_i Q_K. So a direct
## measure of how much information truncation discards is the fraction of each connectome's
## Frobenius energy captured by that block:
##
##     energy(K) = mean_i  || Q_K^T A_i Q_K ||_F^2  /  || A_i ||_F^2
##
## (Q = Laplacian eigenbasis of the mean template graph; energy in [0,1], rising in K.)
## A higher curve at small K => more spectrally compressible => the low-rank graph-spectral
## representation is better justified. Hypothesis: STRUCTURAL connectomes (OASIS-3, braingraph)
## are more compressible than FUNCTIONAL ones (ABIDE, ACPI).
##
## Needs NO outcome labels -- usable on the braingraph OASIS-3 data even without demographics.
## Outputs results/spectral_compressibility.csv + results/figures/fig5_compressibility.png
## =====================================================================================

suppressMessages({
  source("Codes/pipelines/braingraph_oasis_pipeline.R")   # also sources abide_pipeline.R (load_abide, laplacian_basis)
})

## energy curve for a list of V x V connectomes over a grid of K (number of modes)
energy_curve <- function(conn, Ks) {
  V <- nrow(conn[[1]])
  Wbar <- Reduce(`+`, lapply(conn, abs)) / length(conn); diag(Wbar) <- 0
  Q <- laplacian_basis(Wbar)$vectors                 # V x V orthonormal Laplacian basis
  ## per-connectome cumulative block energy
  num <- numeric(length(Ks)); den <- 0
  for (A in conn) {
    M <- crossprod(Q, A %*% Q)                       # Q^T A Q  (V x V)
    s2 <- M^2
    den <- den + sum(s2)
    ## cumulative sum of squares of leading k x k block, for each k in Ks
    csum <- function(k) sum(s2[1:k, 1:k])
    num <- num + vapply(Ks, csum, numeric(1))
  }
  num / den
}

## threshold-crossing: smallest K (and K/V) capturing >= p of total energy
k_for <- function(Ks, energy, V, p) {
  i <- which(energy >= p)[1]
  if (is.na(i)) return(c(K = NA, frac = NA))
  c(K = Ks[i], frac = round(Ks[i] / V, 3))
}

run_compressibility <- function() {
  datasets <- list()
  ## functional (already built, connectomes_full format)
  if (file.exists("data/abide/meta.csv"))
    datasets[["ABIDE (func, V=200)"]] <- load_abide("data/abide")$conn
  if (file.exists("data/acpi/meta.csv"))
    datasets[["ACPI (func, V=200)"]]  <- load_abide("data/acpi")$conn
  ## structural (braingraph OASIS-3), one session per subject
  for (s in 1:3) {
    d <- sprintf("data/braingraph/oasis3_graphmls_scale%d", s)
    if (dir.exists(d) && length(list.files(d, pattern = "[.]graphml$")))
      datasets[[sprintf("OASIS-3 (struct, scale%d)", s)]] <- load_braingraph(d)$conn
  }
  stopifnot(length(datasets) >= 1)

  rows <- list(); curves <- list()
  for (nm in names(datasets)) {
    conn <- datasets[[nm]]; V <- nrow(conn[[1]])
    Ks <- sort(unique(pmin(c(2, 5, 10, 15, 20, 25, 30, 40, 50, 60, 80, 100), V - 1)))
    e <- energy_curve(conn, Ks)
    curves[[nm]] <- list(Ks = Ks, e = e, V = V)
    for (p in c(0.5, 0.8, 0.9)) {
      kf <- k_for(Ks, e, V, p)
      rows[[length(rows) + 1]] <- data.frame(dataset = nm, V = V, N = length(conn),
        threshold = p, K = unname(kf["K"]), K_over_V = unname(kf["frac"]))
    }
  }
  df <- do.call(rbind, rows)
  dir.create("results", showWarnings = FALSE)
  write.csv(df, "results/spectral_compressibility.csv", row.names = FALSE)

  ## figure: energy captured vs fraction of spectrum used (K/V)
  dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
  png("results/figures/fig5_compressibility.png", width = 1300, height = 950, res = 150)
  par(mar = c(5, 5, 4, 2))
  cols <- c("#1b9e77", "#66c2a5", "#7570b3", "#d95f02", "#e7298a")
  plot(NA, xlim = c(0, 0.35), ylim = c(0, 1), xlab = "fraction of graph spectrum used (K / V)",
       ylab = "fraction of connectome energy captured",
       main = "Graph-spectral compressibility of real connectomes")
  abline(h = c(0.5, 0.8, 0.9), col = "grey90", lty = 3)
  nm_all <- names(curves)
  for (i in seq_along(nm_all)) {
    cv <- curves[[nm_all[i]]]
    lines(cv$Ks / cv$V, cv$e, col = cols[i], lwd = 2, type = "b", pch = 19, cex = 0.7)
  }
  legend("right", legend = nm_all, col = cols[seq_along(nm_all)], lwd = 2, pch = 19, bty = "n", cex = 0.8)
  text(0.30, 0.70, "ACPI spike = retained\nglobal signal (DC mode),\nnot predictive low-rank", cex = 0.62, col = "#1f9e89")
  mtext("Diagnostic: at K/V~0.05-0.15 (the method's K=10-30), real connectomes keep little low-frequency energy -> low-rank truncation is lossy here.",
        side = 1, line = 3.6, cex = 0.7)
  dev.off()

  ## energy captured at the K values the method actually uses (more interpretable than
  ## threshold-crossing, which most datasets never reach in the low-frequency basis)
  cat("\n===== Energy captured at the method's typical K (Laplacian low-freq basis) =====\n")
  atK <- do.call(rbind, lapply(names(curves), function(nm) {
    cv <- curves[[nm]]
    g <- function(k) { i <- which(cv$Ks == k); if (length(i)) round(cv$e[i], 3) else NA }
    data.frame(dataset = nm, V = cv$V, `K=10` = g(10), `K=20` = g(20), `K=30` = g(30),
               check.names = FALSE)
  }))
  print(atK, row.names = FALSE)

  cat("\n===== K (and K/V) to capture energy thresholds =====\n")
  print(df, row.names = FALSE)
  cat("\nCAVEAT: ACPI retains the global signal (not global-signal-regressed like ABIDE),\n")
  cat("  so its low-K spike is a single dominant DC mode, not low-rank predictive structure.\n")
  cat("  Net: at K=10-30, real connectomes (structural included) keep only a modest fraction\n")
  cat("  of energy in low-frequency Laplacian modes -- consistent with the supervised finding\n")
  cat("  that low-rank graph-spectral truncation discards predictive signal for these phenotypes.\n")
  cat("Saved results/spectral_compressibility.csv and results/figures/fig5_compressibility.png\n")
  invisible(list(df = df, curves = curves, atK = atK))
}

if (sys.nframe() == 0) run_compressibility()
