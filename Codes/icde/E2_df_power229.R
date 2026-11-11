## =====================================================================================
## E2 / E2b: effective degrees of freedom vs nominal dimension.
##   Rscript Codes/icde/E2_df_power229.R [--smoke]                                   (E2:  Power-229, nbk_0v1 + age_continuous)
##   Rscript Codes/icde/E2_df_power229.R --atlas=sch100 --outcomes=age_continuous     (E2b: Schaefer-100 age; nbk_0v1 there is
##                                                                                      already in results/model_df_comparison.csv)
## Same 10 subject-level splits (seed 900 + split, test share 0.30) and K = 10 adjacency basis as
## Codes/analysis/model_df_comparison.R. Estimators:
##   binary (nbk_0v1), metric auc:
##     two-sided ARD logistic : MacKay sum(gamma_j) at the IRLS mode (intercept counts 1)
##     CV ridge (logistic)    : 1 + trace of the weighted ridge hat matrix, n x n dual Gram eigenvalues
##     kernel ridge (Kundu)   : trace K (K + lambda I)^-1 from the RBF Gram eigenvalues
##   continuous (age_continuous, subject-mean connectomes), metric spearman:
##     two-sided ARD (ard_fit): 1 + MacKay sum(gamma_j) at the returned precisions and noise precision
##     CV ridge (dual form)   : 1 + sum d/(d + alpha_CV) over train Gram eigenvalues (raw_edge_ridge_cv)
##     kernel ridge (Kundu)   : 1 + trace K (K + lambda I)^-1
## Method names match results/model_df_comparison.csv. Ridge collapses to a constant on nbk_0v1 at V = 229
## (df near 1, pred_sd near 0); that is reported as is, which is why age_continuous is added.
## Output: results/icde/df_<atlas>.csv (df_power229.csv for E2, df_sch100.csv for E2b)
##   columns: atlas,outcome,split,method,nominal_dim,eff_df,metric_name,metric,seconds,pred_sd,penalty
## =====================================================================================
source("Codes/icde/icde_common.R")
smoke <- has_flag("--smoke")
atlas <- arg_val("--atlas", "power229")
outcomes <- strsplit(arg_val("--outcomes", "nbk_0v1,age_continuous"), ",")[[1]]
n_splits <- if (smoke) 1 else 10
K <- 10
out <- out_path(paste0("df_", atlas, ".csv"), smoke); if (file.exists(out)) file.remove(out)

for (oc in outcomes) {
  d <- load_set(atlas, oc)
  conn <- d$conn; y <- d$y; subj <- d$subj; bin <- d$type == "binary"
  cat(sprintf("%s %s: %d graphs, %d subjects, V=%d\n", atlas, oc, length(conn), length(unique(subj)), d$V))
  Xe <- build_edge_features(conn)$X_edge             # N x p, independent of the split
  for (s in seq_len(n_splits)) {
    sp <- make_split(subj, s); tr <- sp$tr; te <- sp$te
    if (bin && (length(unique(y[tr])) < 2 || length(unique(y[te])) < 2)) { cat("split", s, "degenerate, skipped\n"); next }
    Xs <- build_spec_features(conn, tr, K, "adjacency")$X_spec
    one <- function(method, fn) {
      t0 <- proc.time()[3]
      r <- tryCatch(fn(), error = function(e) { cat("  ", method, "failed:", conditionMessage(e), "\n"); NULL })
      if (is.null(r)) return(invisible(NULL))
      sc <- if (bin) r$prob else r$pred
      row <- data.frame(atlas = atlas, outcome = oc, split = s, method = method, nominal_dim = r$nominal,
                        eff_df = r$eff_df, metric_name = if (bin) "auc" else "spearman",
                        metric = if (bin) auc_pred(sc, y[te]) else spearman(y[te], sc),
                        seconds = proc.time()[3] - t0, pred_sd = sd(sc), penalty = r$penalty)
      append_csv(row, out)
      cat(sprintf("%s split %2d %-20s nominal=%6d eff_df=%7.2f %s=%.3f pred_sd=%.3f\n", oc, s, method,
                  row$nominal_dim, row$eff_df, row$metric_name, row$metric, row$pred_sd))
    }
    if (bin) {
      one("gs_ard_logit",       function() gs_ard_logit_with_df(Xs[tr, ], y[tr], Xs[te, ]))
      one("raw_edge_ridge",     function() ridge_logit_with_df(Xe[tr, ], y[tr], Xe[te, ], seed = 900 + s))
      one("kundu_kernel_ridge", function() kundu_logit_with_df(conn[tr], y[tr], conn[te]))
    } else {
      one("gs_ard",             function() gs_ard_reg_with_df(Xs[tr, ], y[tr], Xs[te, ]))
      one("raw_edge_ridge",     function() ridge_reg_with_df(Xe[tr, ], y[tr], Xe[te, ]))
      one("kundu_kernel_ridge", function() kundu_reg_with_df(conn[tr], y[tr], conn[te]))
    }
  }
}
cat("Wrote", out, "\n")
