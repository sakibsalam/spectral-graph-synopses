## =====================================================================================
## Rerun of the (outcome, variant) combinations the q >= n_train guard blocked in
## Codes/analysis/dense_variants.R, with that guard removed entirely.
##
## Coordinator correction (see task thread): the guard was wrong. graph_spectral_ard() /
## graph_spectral_ard_logit() are ARD shrinkage estimators -- ard_fit()'s inner solve is
## `alpha_noise * crossprod(Xa) + diag(la)`, and the added diagonal (the ARD prior) makes that
## system invertible at ANY q, including q > n_train. Unpenalized least squares is what breaks
## at q >= n; ARD is specifically the estimator that does not. The existing results already
## demonstrate this is not just a theoretical point: `kundu` (800 features) fits cleanly at
## n_train=101 for age_continuous and beats `current` (55 features, Spearman 0.797) at 0.849,
## p=0.002 -- i.e. a q >> n regime the guard was preventing our own variants from even entering.
##
## Blocked combinations, rerun here with the guard removed:
##   nbk_0v1:         onesided (1000), bigK20 (210), bigK30 (465), bigK40 (820)
##   age_continuous:  onesided (1000), bigK20 (210), bigK30 (465), bigK40 (820), multibasis (165)
## (onesided_rownorm, current, multibasis on nbk_0v1, ridge, kundu already fit under the guard
## and are NOT rerun here -- those results in results/dense_variants_sch100.csv stand.)
##
## Each fit is wrapped in tryCatch so a GENUINE numerical failure (singular solve that even the
## ARD-regularized system cannot resolve, non-convergence, NaN/Inf predictions) is caught and
## recorded with its actual error message, rather than silently skipped on dimension alone.
##
## Same splits, same seeds (900 + split), same subject-level 30% test scheme, same estimators
## (graph_spectral_ard / graph_spectral_ard_logit) as dense_variants.R -- only the guard differs.
##
## Output: results/dense_variants_sch100_nogard.csv (new rows only, resumable) and
## results/dense_variants_summary_nogard.csv (means + paired tests over the COMBINED guarded +
## no-guard data, via summarise_dense_variants() reused from dense_variants.R). The original
## results/dense_variants_sch100.csv and results/dense_variants_summary.csv are untouched.
## =====================================================================================

suppressMessages(source("Codes/analysis/dense_variants.R"))   # functions only -- guarded driver does not run

OUT_CSV_NG     <- "results/dense_variants_sch100_nogard.csv"
SUMMARY_CSV_NG <- "results/dense_variants_summary_nogard.csv"
COMBINED_CSV   <- "results/dense_variants_sch100_combined.csv"     # guarded rows + no-guard rows, for the summary
COMBINED_SUMMARY_CSV <- "results/dense_variants_summary_nogard.csv"

## ---- combinations to rerun without the guard ------------------------------------------
BLOCKED <- list(
  list(outcome = "nbk_0v1",        variant = "onesided"),
  list(outcome = "nbk_0v1",        variant = "bigK20"),
  list(outcome = "nbk_0v1",        variant = "bigK30"),
  list(outcome = "nbk_0v1",        variant = "bigK40"),
  list(outcome = "age_continuous", variant = "onesided"),
  list(outcome = "age_continuous", variant = "bigK20"),
  list(outcome = "age_continuous", variant = "bigK30"),
  list(outcome = "age_continuous", variant = "bigK40"),
  list(outcome = "age_continuous", variant = "multibasis")
)

get_outcome_set <- function(m, outcome) {
  if (outcome == "nbk_0v1") {
    set <- make_set_m(m, conds = c("nbk_0back", "nbk_1back"))
    y <- as.integer(set$meta$COND == "nbk_1back")
    type <- "binary"
  } else if (outcome == "age_continuous") {
    set <- make_set_m(m, conds = NULL, aggregate_by_subject = TRUE)
    y <- as.numeric(set$meta$AGE)
    type <- "regression"
  } else stop("unknown outcome: ", outcome)
  list(conn = set$conn, meta = set$meta, y = y, type = type, V = nrow(set$conn[[1]]))
}

## Fit ONE variant on ONE split, no dimension guard. Returns a one-row data.frame; on a
## genuine numerical failure, status="failed" and reason holds the actual error message
## (never a pre-emptive dimension skip -- that guard is gone).
fit_variant_nogard <- function(outcome, type, variant, conn, y, tr, te) {
  out <- tryCatch({
    t0 <- proc.time()[3]
    X <- build_variant_matrix(conn, tr, variant)
    q <- ncol(X)
    if (type == "regression") {
      f <- graph_spectral_ard(X[tr, , drop = FALSE], y[tr], X[te, , drop = FALSE])
      secs <- proc.time()[3] - t0
      pred <- f$pred
      if (any(!is.finite(pred))) stop(sprintf("non-finite predictions (%d of %d)", sum(!is.finite(pred)), length(pred)))
      data.frame(outcome = outcome, split = NA, variant = variant, dim = q, status = "ok",
                 reason = NA_character_, rmse = rmse(y[te], pred),
                 spearman = suppressWarnings(cor(y[te], pred, method = "spearman")),
                 auc = NA_real_, acc = NA_real_, seconds = secs,
                 n_selected = sum(abs(f$coef) > 1e-6), pred_sd = sd(pred),
                 stringsAsFactors = FALSE)
    } else {
      f <- graph_spectral_ard_logit(X[tr, , drop = FALSE], y[tr], X[te, , drop = FALSE])
      secs <- proc.time()[3] - t0
      prob <- f$prob
      if (any(!is.finite(prob))) stop(sprintf("non-finite predicted probabilities (%d of %d)", sum(!is.finite(prob)), length(prob)))
      data.frame(outcome = outcome, split = NA, variant = variant, dim = q, status = "ok",
                 reason = NA_character_, rmse = NA_real_, spearman = NA_real_,
                 auc = auc_pred(prob, y[te]), acc = mean((prob > 0.5) == y[te]),
                 seconds = secs, n_selected = f$nsel, pred_sd = sd(prob),
                 stringsAsFactors = FALSE)
    }
  }, error = function(e) {
    data.frame(outcome = outcome, split = NA, variant = variant, dim = NA_integer_, status = "failed",
               reason = paste0("ERROR: ", conditionMessage(e)),
               rmse = NA_real_, spearman = NA_real_, auc = NA_real_, acc = NA_real_,
               seconds = NA_real_, n_selected = NA_integer_, pred_sd = NA_real_,
               stringsAsFactors = FALSE)
  })
  out
}

run_nogard <- function(m, out_csv = OUT_CSV_NG) {
  ## group blocked combos by outcome so conn/y/splits are built once per outcome
  outcomes <- unique(sapply(BLOCKED, `[[`, "outcome"))
  for (oc in outcomes) {
    oset <- get_outcome_set(m, oc)
    conn <- oset$conn; y <- oset$y; type <- oset$type
    subj <- oset$meta$SUBJECT; usub <- unique(subj)
    variants_here <- sapply(BLOCKED[sapply(BLOCKED, `[[`, "outcome") == oc], `[[`, "variant")
    cat(sprintf("\n[%s] rerunning without guard: %s\n", oc, paste(variants_here, collapse = ", ")))

    done <- character(0)
    if (file.exists(out_csv)) {
      prev <- read.csv(out_csv, stringsAsFactors = FALSE)
      done <- paste(prev$outcome[prev$outcome == oc], prev$split[prev$outcome == oc],
                    prev$variant[prev$outcome == oc])
    }

    for (s in seq_len(N_SPLITS)) {
      set.seed(SEED0 + s)
      te_s <- sample(usub, max(1, round(TEST_FRAC * length(usub))))
      te <- which(subj %in% te_s); tr <- which(!(subj %in% te_s))
      if (type != "regression" && (length(unique(y[tr])) < 2 || length(unique(y[te])) < 2)) {
        cat(sprintf("  split %d: degenerate label split -- skipping\n", s)); next
      }
      for (variant in variants_here) {
        key <- paste(oc, s, variant)
        if (key %in% done) { cat(sprintf("  split %d %s already done -- skipping\n", s, variant)); next }
        r <- fit_variant_nogard(oc, type, variant, conn, y, tr, te)
        r$split <- s
        prev <- if (file.exists(out_csv)) read.csv(out_csv, stringsAsFactors = FALSE) else NULL
        dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
        write.csv(rbind(prev, r), out_csv, row.names = FALSE)
        if (r$status == "ok") {
          val <- if (type == "regression") sprintf("spearman=%.3f rmse=%.2f", r$spearman, r$rmse) else sprintf("auc=%.3f", r$auc)
          cat(sprintf("  split %d %-10s dim=%d %s (%.2fs)\n", s, variant, r$dim, val, r$seconds))
        } else {
          cat(sprintf("  split %d %-10s FAILED: %s\n", s, variant, r$reason))
        }
        flush.console()
      }
    }
  }
  invisible(TRUE)
}

m <- readRDS(MASTER_PATH)
cat(sprintf("Loaded %s: %d connectomes, V=%d, %d subjects\n",
            MASTER_PATH, length(m$conn), m$V, length(unique(m$subjects$SUBJECT))))

run_nogard(m)

## ---- combine guarded + no-guard results for the summary --------------------------------
## The original guarded run recorded the blocked combos as status="skipped" -- replace those
## rows with the new no-guard fit (status="ok" or "failed") so the combined table has exactly
## one row per (outcome, split, variant), then reuse dense_variants.R's own summariser.
d_orig <- read.csv(OUT_CSV, stringsAsFactors = FALSE)
d_ng    <- if (file.exists(OUT_CSV_NG)) read.csv(OUT_CSV_NG, stringsAsFactors = FALSE) else d_orig[0, ]

key_orig <- paste(d_orig$outcome, d_orig$split, d_orig$variant)
key_ng   <- paste(d_ng$outcome, d_ng$split, d_ng$variant)
d_combined <- rbind(d_orig[!(key_orig %in% key_ng), ], d_ng)
d_combined <- d_combined[order(d_combined$outcome, d_combined$split, d_combined$variant), ]
write.csv(d_combined, COMBINED_CSV, row.names = FALSE)
cat("\nWrote", COMBINED_CSV, "(", nrow(d_combined), "rows )\n")

summ <- summarise_dense_variants(out_csv = COMBINED_CSV, summary_csv = COMBINED_SUMMARY_CSV)
cat("\n==== Combined summary (guarded arms unchanged + no-guard reruns) ====\n")
print(summ, row.names = FALSE)

## ---- the two comparisons the coordinator asked to see explicitly -----------------------
cat("\n==== 1. onesided (1000, node identity kept) vs current (55) ====\n")
print(summ[summ$variant %in% c("current", "onesided"),
           c("outcome", "variant", "dim", "n_splits_ok", "mean_spearman", "mean_auc",
             "diff_vs_current", "wilcoxon_p_vs_current")], row.names = FALSE)

cat("\n==== 2. bigK40 (820, node identity discarded) vs onesided (1000, node identity kept) ====\n")
for (oc in c("nbk_0v1", "age_continuous")) {
  dbig <- d_combined[d_combined$outcome == oc & d_combined$variant == "bigK40" & d_combined$status == "ok", ]
  dones <- d_combined[d_combined$outcome == oc & d_combined$variant == "onesided" & d_combined$status == "ok", ]
  primary <- if (oc == "age_continuous") "spearman" else "auc"
  mrg <- merge(dbig[, c("split", primary)], dones[, c("split", primary)], by = "split", suffixes = c("_bigK40", "_onesided"))
  cat(sprintf("\n[%s] bigK40 n_ok=%d, onesided n_ok=%d, paired n=%d\n", oc, nrow(dbig), nrow(dones), nrow(mrg)))
  if (nrow(mrg) >= 3) {
    wt <- wilcox.test(mrg[[paste0(primary, "_bigK40")]], mrg[[paste0(primary, "_onesided")]], paired = TRUE)
    cat(sprintf("  mean bigK40 %s = %.4f, mean onesided %s = %.4f, diff = %.4f, Wilcoxon p = %.4f\n",
                primary, mean(mrg[[paste0(primary, "_bigK40")]]), primary, mean(mrg[[paste0(primary, "_onesided")]]),
                mean(mrg[[paste0(primary, "_bigK40")]] - mrg[[paste0(primary, "_onesided")]]), wt$p.value))
  } else {
    cat("  not enough paired fits to test\n")
  }
}
