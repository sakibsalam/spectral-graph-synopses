## G9: template rank K sweep for the two-sided synopsis with dual ridge (CV penalty).
##   Rscript Codes/icde/G9_ksweep_ridge.R [--smoke] [--atlas=sch100]
## Same splits / nested training fractions / scoring as G1. Output: results/icde/ksweep_ridge.csv
## (with --atlas=<a>, results/icde/ksweep_ridge_<a>.csv so atlases can run in parallel).
source("Codes/icde/icde_common.R")
smoke   <- has_flag("--smoke")
atlases <- if (is.na(arg_val("--atlas"))) c("sch100", "power229") else arg_val("--atlas")
outcomes <- c("age_continuous", "nbk_0v1")
fracs   <- if (smoke) 1 else c(0.25, 0.5, 0.75, 1.0)
n_splits <- if (smoke) 1 else 10
Ks <- c(5, 10, 15, 20, 30)
out <- out_path(if (length(atlases) == 1) sprintf("ksweep_ridge_%s.csv", atlases) else "ksweep_ridge.csv", smoke)
if (file.exists(out)) file.remove(out)

for (atlas in atlases) for (oc in outcomes) {
  d <- load_set(atlas, oc); conn <- d$conn; y <- d$y; subj <- d$subj
  binary <- d$type == "binary"
  metric_name <- if (binary) "auc" else "spearman"
  cat(sprintf("[%s %s] %d graphs\n", atlas, oc, length(conn)))
  for (s in seq_len(n_splits)) {
    sp <- make_split(subj, s); te <- sp$te
    if (binary && length(unique(y[te])) < 2) next
    rows <- list()
    for (fr in fracs) {
      tr <- subsample_train(sp$tr, subj, fr, s)
      if (binary && length(unique(y[tr])) < 2) next
      for (K in Ks) {
        Q <- template_Q(conn, tr, K); X <- feat_twosided(conn, Q)
        pr <- ridge_dual_cv_folds(X[tr, , drop = FALSE], y[tr], X[te, , drop = FALSE], subject_folds(subj[tr]))$pred
        rows[[length(rows) + 1]] <- data.frame(atlas = atlas, outcome = oc, K = K, n_features = ncol(X),
          train_frac = fr, split = s, n_train = length(tr), metric_name = metric_name,
          metric = if (binary) auc_pred(pr, y[te]) else spearman(y[te], pr),
          rmse = if (binary) NA else rmse(y[te], pr))
      }
    }
    if (length(rows)) append_csv(do.call(rbind, rows), out)
    cat(sprintf("  split %d done\n", s)); flush.console()
  }
}
cat("Wrote", out, "\n")
