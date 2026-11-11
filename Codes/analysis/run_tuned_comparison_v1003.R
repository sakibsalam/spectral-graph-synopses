## =====================================================================================
## v1003 rerun of Codes/analysis/tuned_comparison.R on the completed OSF cache (1003 files,
## 144 subjects -- up from 958/139). Does NOT modify tuned_comparison.R; sources it (which,
## per the sys.nframe()==0 guard at its bottom, does not auto-run when sourced) and re-invokes
## its own exported functions with new output paths, so nothing under results/ is overwritten
## and no statistical setting, seed, or split scheme changes -- only the input sample size.
##
## The one thing that has to change here: tuned_comparison.R's EXPECTED row/subject-count
## guardrail is a data-integrity check, not a statistical setting, and it is hardcoded to the
## OLD counts (nbk_0v1=272/144, nbk_1v2=276/144). At the completed cache size every condition
## has exactly 144 files, so nbk_0v1 and nbk_1v2 both become 288/144. We overwrite the EXPECTED
## global (created in the global environment when this script sources tuned_comparison.R) to
## the new counts rather than editing the original file, so a future rerun of
## `Rscript Codes/analysis/tuned_comparison.R` against old/partial data still checks against
## the original expectations.
## =====================================================================================

suppressMessages({
  source("Codes/analysis/tuned_comparison.R")   # functions only -- does not auto-run (sourced, not Rscript'd directly)
})

EXPECTED <- list(
  nbk_0v1           = list(n = 288, subj = 144),
  nbk_1v2           = list(n = 288, subj = 144),
  age_continuous    = list(n = 144, subj = 144),
  age_binary_median = list(n = 144, subj = 144)
)

OUT_MAIN    <- "results/tuned_comparison_v1003.csv"
OUT_SUMMARY <- "results/tuned_comparison_summary_v1003.csv"
OUT_PAIRED  <- "results/tuned_comparison_paired_tests_v1003.csv"
OUT_PERMUTED <- "results/tuned_comparison_permuted_v1003.csv"

main_v1003 <- function() {
  outcomes <- list(
    nbk_0v1           = "binary",
    nbk_1v2           = "binary",
    age_continuous     = "regression",
    age_binary_median  = "binary"
  )

  cat("\n================ v1003 JOB 3: tuned comparison, four outcomes ================\n")
  for (nm in names(outcomes)) {
    cat(sprintf("\n---- %s (%s) ----\n", nm, outcomes[[nm]]))
    run_outcome_tuned(nm, n_splits = 10, test_frac = 0.30, seed0 = 900, out_csv = OUT_MAIN)
  }

  main_df <- read.csv(OUT_MAIN)
  main_df <- main_df[main_df$permuted == FALSE, ]

  summ <- build_summary(main_df)
  write.csv(summ, OUT_SUMMARY, row.names = FALSE)
  cat("\nWrote", OUT_SUMMARY, "\n")
  print(summ, row.names = FALSE)

  pt <- do.call(rbind, lapply(names(outcomes), function(nm)
    paired_tests_for_outcome(main_df, nm, outcomes[[nm]])))
  write.csv(pt, OUT_PAIRED, row.names = FALSE)
  cat("\nWrote", OUT_PAIRED, "\n")
  print(pt, row.names = FALSE)

  cat("\n================ v1003 JOB 4: permutation (leakage) control ================\n")
  run_permutation_control(n_splits = 6, out_csv = OUT_PERMUTED)

  leak <- check_permutation_leak(perm_csv = OUT_PERMUTED)
  cat(sprintf("\n---- permuted age_continuous: one-sample t-test of Spearman vs 0, by arm (%d splits, n_test ~= %d subjects) ----\n",
              leak$n_splits, leak$n_test_ref))
  cat(sprintf("Theoretical null reference SD (1/sqrt(n_test-1)): %.4f -- compare to each arm's sd_spearman below.\n",
              leak$sd_null_ref))
  print(leak$table, row.names = FALSE)
  if (leak$leaked) {
    cat("\n*** LEAK DETECTED (v1003): at least one arm's mean permuted Spearman is POSITIVE and",
        sprintf("significantly greater than zero (one-sided p < %.2f). STOP -- investigate before", ALPHA_LEAK),
        "trusting", OUT_SUMMARY, "***\n")
  } else {
    cat("\nNo leak detected (v1003): no arm's mean permuted Spearman is significantly positive at",
        sprintf("p < %.2f.\n", ALPHA_LEAK))
  }

  invisible(list(summary = summ, paired_tests = pt, leak = leak))
}

if (sys.nframe() == 0) {
  main_v1003()
}
