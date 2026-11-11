## =====================================================================================
## G5: does forming the inner cross-validation folds by subject change the ridge arms on the
## binary task? Each subject contributes two graphs with opposite labels, and the default
## dual-Gram CV ridge forms its folds over graphs.
##   Rscript Codes/icde/G5_grouped_cv_t1.R --atlas=sch100 [--smoke]
## Arms: edge, one-sided (R = 5, 10), two-sided (K = 10) and edge-space PCA (55), each with
## dual ridge and folds over graphs ("graph") or over subjects ("subject"). Same splits and
## nested training fractions as E3. Output: results/icde/grouped_cv_<atlas>.csv
## =====================================================================================
source("Codes/icde/icde_common.R")
smoke <- has_flag("--smoke"); atlas <- arg_val("--atlas", "sch100")
fracs <- if (smoke) 1 else c(0.25, 0.5, 0.75, 1.0); n_splits <- if (smoke) 1 else 10
out <- out_path(sprintf("grouped_cv_%s.csv", atlas), smoke)
if (file.exists(out)) file.remove(out)


d <- load_set(atlas, "nbk_0v1"); conn <- d$conn; y <- d$y; subj <- d$subj
Xe <- build_edge_features(conn)$X_edge
for (s in seq_len(n_splits)) {
  sp <- make_split(subj, s); te <- sp$te; rows <- list()
  for (fr in fracs) {
    tr <- subsample_train(sp$tr, subj, fr, s)
    set.seed(1); f_graph <- sample(rep(seq_len(5), length.out = length(tr)))          # as raw_edge_ridge_cv
    f_subj <- subject_folds(subj[tr])
    Q10 <- template_Q(conn, tr, 10)
    feats <- list(edge = Xe, one_sided_5 = feat_onesided(conn, Q10[, 1:5]), one_sided_10 = feat_onesided(conn, Q10),
                  two_sided_10 = feat_twosided(conn, Q10), pca_55 = pca_scores(Xe, tr, 55))
    for (nm in names(feats)) for (fold in c("graph", "subject")) {
      r <- ridge_dual_cv_folds(feats[[nm]][tr, , drop = FALSE], y[tr], feats[[nm]][te, , drop = FALSE],
                               if (fold == "graph") f_graph else f_subj)
      rows[[length(rows) + 1]] <- data.frame(atlas = atlas, arm = nm, folds = fold, train_frac = fr, split = s,
        n_train = length(tr), auc = auc_pred(r$pred, y[te]), alpha_index = r$alpha_index, pred_sd = sd(r$pred))
    }
  }
  append_csv(do.call(rbind, rows), out); cat("  split", s, "done\n"); flush.console()
}
cat("Wrote", out, "\n")
