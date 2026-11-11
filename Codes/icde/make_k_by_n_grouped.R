## Training-size grid with the ridge arms of the binary task refitted with inner folds by subject.
##   Rscript Codes/icde/make_k_by_n_grouped.R
## Copies results/icde/k_by_n.csv and replaces the metric of the nbk_0v1 one-sided ridge (R = 5, 10)
## and edge ridge rows by the subject-fold values of G5 (results/icde/grouped_cv_<atlas>.csv). Two-sided
## ARD rows (no inner CV) and all age rows (one graph per subject) are unchanged. The logistic
## wide-path reference rows are dropped. Output: results/icde/k_by_n_grouped.csv
kb <- read.csv("results/icde/k_by_n.csv")
q <- do.call(rbind, lapply(list.files("results/icde", "^grouped_cv_.*\\.csv$", full.names = TRUE), read.csv))
q <- q[q$folds == "subject", ]
kb <- kb[kb$design != "edge_ridge_glmnet_wide", ]
arm <- ifelse(kb$design == "one_sided_ridge", paste0("one_sided_", kb$K_or_R), ifelse(kb$design == "edge_ridge_cv", "edge", NA))
i <- which(kb$outcome == "nbk_0v1" & !is.na(arm))
j <- match(paste(kb$atlas[i], arm[i], kb$train_frac[i], kb$split[i]), paste(q$atlas, q$arm, q$train_frac, q$split))
stopifnot(!anyNA(j), all(kb$n_train[i] == q$n_train[j]))
kb$metric[i] <- q$auc[j]; kb$seconds[i] <- NA
write.csv(kb, "results/icde/k_by_n_grouped.csv", row.names = FALSE)
cat("replaced", length(i), "rows; wrote results/icde/k_by_n_grouped.csv\n")
