## =====================================================================================
## E3: synopsis size x training-set size grid (both atlases, two outcomes).
##   Rscript Codes/icde/E3_k_by_n.R [--smoke] [--no-glmnet-ref]
## Designs: two-sided ARD (K in 5,10,15,20,30), one-sided ridge on A U_R (R in 5,10; dual CV
## ridge), edge ridge CV reference (dual form), and for the binary outcome the wide-path logistic
## ridge (the configuration that collapses to a constant at Power-229; kept as a second reference).
## Training subjects are subsampled within each split's training set (nested, fixed seed);
## template, standardization and ridge penalty use the subsampled training graphs only.
## Resumable: a finished (atlas, outcome, split) block is skipped on rerun.
## =====================================================================================
source("Codes/icde/icde_common.R")
smoke <- has_flag("--smoke"); glm_ref <- !has_flag("--no-glmnet-ref") && !smoke
atlases  <- if (smoke) "sch100" else c("power229", "sch100")
outcomes <- c("nbk_0v1", "age_continuous")
Ks   <- if (smoke) c(5, 10) else c(5, 10, 15, 20, 30)
Rs   <- if (smoke) 5 else c(5, 10)
fracs <- if (smoke) c(0.5, 1) else c(0.25, 0.5, 0.75, 1.0)
n_splits <- if (smoke) 1 else 10
out <- out_path("k_by_n.csv", smoke)
done <- if (file.exists(out)) unique(with(read.csv(out), paste(atlas, outcome, split))) else character(0)

for (atlas in atlases) for (oc in outcomes) {
  d <- load_set(atlas, oc); conn <- d$conn; y <- d$y; subj <- d$subj
  metric_name <- if (d$type == "binary") "auc" else "spearman"
  score <- function(pred, yy) if (d$type == "binary") auc_pred(pred, yy) else spearman(yy, pred)
  Xe <- build_edge_features(conn)$X_edge
  cat(sprintf("\n[%s %s] %d graphs, %d subjects, V=%d\n", atlas, oc, length(conn), length(unique(subj)), d$V))
  for (s in seq_len(n_splits)) {
    if (paste(atlas, oc, s) %in% done) next
    sp <- make_split(subj, s); te <- sp$te
    if (d$type == "binary" && length(unique(y[te])) < 2) next
    rows <- list()
    add <- function(design, kr, frac, ntr, val, secs)
      rows[[length(rows) + 1]] <<- data.frame(atlas = atlas, outcome = oc, design = design, K_or_R = kr,
        train_frac = frac, split = s, n_train = ntr, metric_name = metric_name, metric = val, seconds = secs)
    fit_try <- function(fn) tryCatch(fn(), error = function(e) { cat("   fit failed:", conditionMessage(e), "\n"); NULL })
    for (fr in fracs) {
      tr <- subsample_train(sp$tr, subj, fr, s)
      if (d$type == "binary" && length(unique(y[tr])) < 2) next
      for (K in Ks) {                                          # two-sided ARD
        t0 <- proc.time()[3]
        r <- fit_try(function() {
          X <- feat_twosided(conn, template_Q(conn, tr, K))
          if (d$type == "binary") graph_spectral_ard_logit(X[tr, ], y[tr], X[te, ])$prob
          else graph_spectral_ard(X[tr, ], y[tr], X[te, ])$pred })
        add("two_sided_ard", K, fr, length(tr), if (is.null(r)) NA else score(r, y[te]), proc.time()[3] - t0)
      }
      for (R in Rs) {                                          # one-sided ridge (dual, CV penalty)
        t0 <- proc.time()[3]
        r <- fit_try(function() {
          X <- feat_onesided(conn, template_Q(conn, tr, R))
          raw_edge_ridge_cv(X[tr, ], y[tr], X[te, ])$pred })
        add("one_sided_ridge", R, fr, length(tr), if (is.null(r)) NA else score(r, y[te]), proc.time()[3] - t0)
      }
      t0 <- proc.time()[3]                                     # edge ridge reference (dual form, CV)
      r <- fit_try(function() raw_edge_ridge_cv(Xe[tr, ], y[tr], Xe[te, ])$pred)
      add("edge_ridge_cv", NA, fr, length(tr), if (is.null(r)) NA else score(r, y[te]), proc.time()[3] - t0)
      if (glm_ref && d$type == "binary") {                     # logistic wide-path ridge (collapse-prone)
        t0 <- proc.time()[3]
        r <- fit_try(function() raw_edge_glmnet_logit_wide(Xe[tr, ], y[tr], Xe[te, ], alpha = 0)$prob)
        add("edge_ridge_glmnet_wide", NA, fr, length(tr), if (is.null(r)) NA else score(r, y[te]), proc.time()[3] - t0)
      }
    }
    append_csv(do.call(rbind, rows), out)
    cat(sprintf("  split %d done (%d rows)\n", s, length(rows))); flush.console()
  }
}
cat("Wrote", out, "\n")
