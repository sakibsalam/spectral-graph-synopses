## =====================================================================================
## G1: generic reductions at the size of the two-sided synopsis, and a learner-matched
## comparison of the three designs.
##   Rscript Codes/icde/G1_generic_synopses.R --atlas=sch100 --outcome=age_continuous [--smoke]
## Arms (all on the E3 splits, seed 900 + split, nested training fractions):
##   two_sided  K = 10 (q = 55)      learners: ard, ridge   (ridge = dual-Gram CV, inner folds by subject)
##   one_sided  R = 10 (V R)         learner : ridge
##   edge       p                    learner : ridge
##   pca        q = 55               learners: ard, ridge   (edge-space PCA, training graphs only)
##   rp         q in 55, 210, 1000   learners: ard (55, 210), ridge (all); 5 Gaussian draws
## PCA keeps min(q, n_train - 1) components. The random projection is data independent:
## Z = X G / sqrt(q), G a p x q standard Gaussian matrix with seed 5000 + 100 split + draw.
## One output file per (atlas, outcome) so the four jobs can run side by side:
##   results/icde/generic_<atlas>_<outcome>.csv
## =====================================================================================
source("Codes/icde/icde_common.R")
smoke   <- has_flag("--smoke")
atlas   <- arg_val("--atlas", "sch100")
oc      <- arg_val("--outcome", "age_continuous")
fracs   <- if (smoke) 1 else c(0.25, 0.5, 0.75, 1.0)
n_splits <- if (smoke) 1 else 10
n_draws <- if (smoke) 1 else 5
rp_sizes <- c(55, 210, 1000)
out <- out_path(sprintf("generic_%s_%s.csv", atlas, oc), smoke)
done <- if (file.exists(out)) unique(read.csv(out)$split) else integer(0)

d <- load_set(atlas, oc); conn <- d$conn; y <- d$y; subj <- d$subj
binary <- d$type == "binary"
metric_name <- if (binary) "auc" else "spearman"
score <- function(pred, yy) if (binary) auc_pred(pred, yy) else spearman(yy, pred)
Xe <- build_edge_features(conn)$X_edge; p <- ncol(Xe)
cat(sprintf("[%s %s] %d graphs, V=%d, p=%d\n", atlas, oc, length(conn), d$V, p))

fit_ard <- function(X, tr, te) {
  if (binary) graph_spectral_ard_logit(X[tr, , drop = FALSE], y[tr], X[te, , drop = FALSE])$prob
  else graph_spectral_ard(X[tr, , drop = FALSE], y[tr], X[te, , drop = FALSE])$pred
}
## dual ridge, CV penalty; inner folds over subjects (one graph per subject on age, two on the binary task)
fit_ridge <- function(X, tr, te) ridge_dual_cv_folds(X[tr, , drop = FALSE], y[tr], X[te, , drop = FALSE], subject_folds(subj[tr]))$pred

for (s in seq_len(n_splits)) {
  if (s %in% done) next
  sp <- make_split(subj, s); te <- sp$te
  if (binary && length(unique(y[te])) < 2) next
  rows <- list()
  run <- function(design, size, learner, draw, fr, tr, X, build_s) {
    t0 <- proc.time()[3]
    pr <- tryCatch(if (learner == "ard") fit_ard(X, tr, te) else fit_ridge(X, tr, te),
                   error = function(e) { cat("   failed:", design, size, learner, conditionMessage(e), "\n"); NULL })
    rows[[length(rows) + 1]] <<- data.frame(atlas = atlas, outcome = oc, design = design, size = size,
      learner = learner, draw = draw, train_frac = fr, split = s, n_train = length(tr), n_features = ncol(X),
      metric_name = metric_name, metric = if (is.null(pr)) NA else score(pr, y[te]),
      rmse = if (is.null(pr) || binary) NA else rmse(y[te], pr),
      build_seconds = build_s, fit_seconds = proc.time()[3] - t0)
  }
  for (fr in fracs) {
    tr <- subsample_train(sp$tr, subj, fr, s)
    if (binary && length(unique(y[tr])) < 2) next
    t0 <- proc.time()[3]; Q <- template_Q(conn, tr, 10); XT <- feat_twosided(conn, Q); bT <- proc.time()[3] - t0
    run("two_sided", 10, "ard", 0, fr, tr, XT, bT); run("two_sided", 10, "ridge", 0, fr, tr, XT, bT)
    t0 <- proc.time()[3]; XO <- feat_onesided(conn, Q); bO <- proc.time()[3] - t0
    run("one_sided", 10, "ridge", 0, fr, tr, XO, bO)
    run("edge", NA, "ridge", 0, fr, tr, Xe, 0)
    t0 <- proc.time()[3]; XP <- pca_scores(Xe, tr, 55); bP <- proc.time()[3] - t0
    run("pca", 55, "ard", 0, fr, tr, XP, bP); run("pca", 55, "ridge", 0, fr, tr, XP, bP)
    for (dr in seq_len(n_draws)) for (q in rp_sizes) {
      t0 <- proc.time()[3]
      set.seed(5000 + 100 * s + dr); G <- matrix(rnorm(p * q), p, q) / sqrt(q)
      XR <- Xe %*% G; bR <- proc.time()[3] - t0
      if (q <= 210) run("rp", q, "ard", dr, fr, tr, XR, bR)
      run("rp", q, "ridge", dr, fr, tr, XR, bR)
    }
  }
  append_csv(do.call(rbind, rows), out)
  cat(sprintf("  split %d done (%d rows)\n", s, length(rows))); flush.console()
}
cat("Wrote", out, "\n")
