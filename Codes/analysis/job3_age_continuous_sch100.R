## =====================================================================================
## Job 3: continuous age at V=100 (Schaefer-100), 10 splits, arms gs_ard/hybrid/ridge/lasso,
## Spearman + RMSE. Same subject-level split scheme (seed0=900, 30% test) as Jobs 1-2.
##
## Reuses run_outcome_m() from Codes/analysis/multi_resolution_common.R (the master-
## parameterized twin of outcome_explorer.R's run_outcome(), built on the SAME eval_split()
## so results are directly comparable to the Power-229 age_continuous numbers). One
## connectome per subject (aggregate_by_subject = TRUE), y = AGE (years).
## =====================================================================================

suppressMessages(source("Codes/analysis/multi_resolution_common.R"))

OUT_CSV <- "results/age_continuous_sch100.csv"
m <- readRDS("data/osf_schaefer100/master.rds")
cat(sprintf("Loaded Schaefer-100 master: %d connectomes, V=%d, %d subjects\n",
            length(m$conn), m$V, length(unique(m$subjects$SUBJECT))))

run_outcome_m(m, name = "age_continuous_sch100", conds = NULL, aggregate_by_subject = TRUE,
             y_col = "AGE", type = "regression", n_splits = 10, test_frac = 0.30, K = 10,
             basis = "adjacency", methods = c("gs_ard", "hybrid", "ridge", "lasso"),
             seed0 = 900, out_csv = OUT_CSV)

d <- read.csv(OUT_CSV)
cat("\n==== Job 3 summary (mean per method, metric) ====\n")
if (nrow(d) > 0) {
  agg <- aggregate(value ~ method + metric, d, mean)
  print(agg[agg$metric %in% c("rmse", "spearman"), ])
}
cat("\nDone. Wrote", OUT_CSV, "\n")
