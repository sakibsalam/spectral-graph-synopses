## =====================================================================================
## G2: the design comparison on two further functional collections (ABIDE, ACPI; V = 200,
## one graph per subject, outcome = age at scan).
##   Rscript Codes/icde/G2_external_collections.R [--smoke]
## Splits: 10 random 70/30 splits, seed 100 + split (the scheme of results/abide_repeated.csv,
## so the two-sided ARD and fixed-penalty ridge rows reproduce that file).
## Arms: two-sided K in 5, 10 (ARD) and K = 10 (dual ridge, CV); one-sided R in 5, 10 (dual
## ridge, CV); edge dual ridge with CV penalty and with the fixed penalty alpha = 50;
## edge-space PCA with 55 components (ARD, ridge); Gaussian random projection to 55 and 1000
## features (ridge, 5 draws); training-mean baseline.
## Output: results/icde/external_collections.csv
## =====================================================================================
source("Codes/icde/icde_common.R")
smoke <- has_flag("--smoke")
n_splits <- if (smoke) 1 else 10
n_draws <- if (smoke) 1 else 5
out <- out_path("external_collections.csv", smoke)
if (file.exists(out)) file.remove(out)

for (coll in c("abide", "acpi")) {
  d <- load_abide(file.path("data", coll))
  y <- d$subjects$AGE_AT_SCAN; keep <- which(is.finite(y)); conn <- d$conn[keep]; y <- y[keep]
  N <- length(conn); Xe <- build_edge_features(conn)$X_edge; p <- ncol(Xe)
  cat(sprintf("[%s] N=%d V=%d p=%d\n", coll, N, d$V, p))
  for (s in seq_len(n_splits)) {
    set.seed(100 + s); te <- sample(seq_len(N), round(0.30 * N)); tr <- setdiff(seq_len(N), te)
    rows <- list()
    add <- function(design, size, learner, draw, pred, nf)
      rows[[length(rows) + 1]] <<- data.frame(collection = coll, V = d$V, split = s, n_train = length(tr),
        design = design, size = size, learner = learner, draw = draw, n_features = nf,
        spearman = spearman(y[te], pred), rmse = rmse(y[te], pred))
    ard <- function(X) graph_spectral_ard(X[tr, , drop = FALSE], y[tr], X[te, , drop = FALSE])$pred
    rdg <- function(X) raw_edge_ridge_cv(X[tr, , drop = FALSE], y[tr], X[te, , drop = FALSE])$pred
    add("baseline", NA, "train_mean", 0, rep(mean(y[tr]), length(te)) + rnorm(length(te), 0, 1e-9), 0)
    for (K in c(5, 10)) {
      Q <- template_Q(conn, tr, K); XT <- feat_twosided(conn, Q); XO <- feat_onesided(conn, Q)
      add("two_sided", K, "ard", 0, ard(XT), ncol(XT))
      if (K == 10) add("two_sided", K, "ridge", 0, rdg(XT), ncol(XT))
      add("one_sided", K, "ridge", 0, rdg(XO), ncol(XO))
    }
    add("edge", NA, "ridge", 0, rdg(Xe), p)
    add("edge", NA, "ridge_fixed50", 0, raw_edge_ridge_woodbury(Xe[tr, ], y[tr], Xe[te, ])$pred, p)
    XP <- pca_scores(Xe, tr, 55)
    add("pca", 55, "ard", 0, ard(XP), ncol(XP)); add("pca", 55, "ridge", 0, rdg(XP), ncol(XP))
    for (dr in seq_len(n_draws)) for (q in c(55, 1000)) {
      set.seed(5000 + 100 * s + dr); XR <- Xe %*% (matrix(rnorm(p * q), p, q) / sqrt(q))
      add("rp", q, "ridge", dr, rdg(XR), q)
    }
    append_csv(do.call(rbind, rows), out)
    cat(sprintf("  split %d done\n", s)); flush.console()
  }
}
cat("Wrote", out, "\n")
