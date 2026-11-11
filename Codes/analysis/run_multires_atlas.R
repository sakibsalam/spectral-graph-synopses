## =====================================================================================
## Driver for Job 4a (nbk_0v1 multi-resolution replication) and Job 4c (continuous age at
## this resolution), for ONE atlas at a time (Schaefer-100 or Schaefer-200). Called with
## Rscript --vanilla, passing the atlas's master.rds path and output suffix as CLI args so the
## exact same code runs both atlases without duplication.
##
## Usage:
##   Rscript Codes/analysis/run_multires_atlas.R data/osf_schaefer100/master.rds sch100
##   Rscript Codes/analysis/run_multires_atlas.R data/osf_schaefer200/master.rds sch200
## =====================================================================================

suppressMessages(source("Codes/analysis/multi_resolution_common.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript run_multires_atlas.R <master.rds path> <suffix>")
MASTER_PATH <- args[1]
SUFFIX <- args[2]

cat(sprintf("\n================ Multi-resolution jobs for %s (suffix=%s) ================\n",
            MASTER_PATH, SUFFIX))
m <- readRDS(MASTER_PATH)
cat(sprintf("Loaded master: %d connectomes, V=%d, %d subjects\n",
            length(m$conn), m$V, length(unique(m$subjects$SUBJECT))))

## ---- Job 4a: nbk_0v1, both ridge configs, 10 splits ----------------------------------
cat("\n---- Job 4a: nbk_0v1 ----\n")
run_nbk_0v1_multires(m, n_splits = 10, seed0 = 900, test_frac = 0.30, K = 10, basis = "adjacency",
                     out_csv = sprintf("results/nbk_0v1_%s.csv", SUFFIX))

## ---- Job 4c: age_continuous, standard 4-method set (gs_ard/hybrid/ridge/lasso), 10 splits ----
cat("\n---- Job 4c: age_continuous ----\n")
run_outcome_m(m, name = "age_continuous", conds = NULL, aggregate_by_subject = TRUE,
             y_col = "AGE", type = "regression", n_splits = 10, seed0 = 900,
             K = 10, basis = "adjacency", out_csv = sprintf("results/age_continuous_%s.csv", SUFFIX))

cat(sprintf("\nDone with %s.\n", SUFFIX))
