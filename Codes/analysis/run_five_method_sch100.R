## =====================================================================================
## Driver for Job 4b: the five-method comparison at V=100 (Schaefer-100), on nbk_0v1 and
## age_continuous. See Codes/analysis/five_method_multires.R for the method implementations
## and the feasibility guard. Output: results/five_method_sch100.csv (both outcomes, long
## format, resumable).
## =====================================================================================

suppressMessages(source("Codes/analysis/five_method_multires.R"))

m <- readRDS("data/osf_schaefer100/master.rds")
cat(sprintf("Loaded Schaefer-100 master: %d connectomes, V=%d, %d subjects\n",
            length(m$conn), m$V, length(unique(m$subjects$SUBJECT))))

OUT <- "results/five_method_sch100.csv"

cat("\n================ Job 4b: five-method comparison, nbk_0v1 (V=100) ================\n")
run_five_method(m, "nbk_0v1", n_splits = 10, seed0 = 900, test_frac = 0.30, K = 10,
                basis = "adjacency", out_csv = OUT, max_seconds_per_fit = 600, n_splits_fallback = 3)

cat("\n================ Job 4b: five-method comparison, age_continuous (V=100) ================\n")
run_five_method(m, "age_continuous", n_splits = 10, seed0 = 900, test_frac = 0.30, K = 10,
                basis = "adjacency", out_csv = OUT, max_seconds_per_fit = 600, n_splits_fallback = 3)

cat("\nDone. Wrote", OUT, "\n")
