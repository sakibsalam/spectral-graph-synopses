## =====================================================================================
## Dense-feature variants: is the graph-spectral reduction (q=55 at V=100) too aggressive?
##
## Our method computes C = U'AU (both sides projected onto the K-dim basis) then
## half-vectorizes to q=K(K+1)/2 features -- a very aggressive reduction that discards node
## identity (every feature is a mode-pair coefficient, not attached to any brain region).
## Kundu's method computes L = AU (one side only), giving V*R features that keep node
## identity. This file implements variants that interpolate between our 55-feature core and
## the raw p=4,950 edge representation, applying Kundu's one-sided move to OUR basis instead
## of his own adjacency embedding, and benchmarks them against `current` (the paper's method),
## `kundu` (Codes/benchmarks/competitor_benchmarks.R::kundu_latent_nonlinear, R=8), and `ridge`
## (raw-edge, Codes/core/abide_pipeline.R::raw_edge_ridge_woodbury /
## Codes/core/abide_classification.R::raw_edge_glmnet_logit -- the project's standard ridge
## comparator, unchanged).
##
## Variants (all built from the SAME train-only template basis as build_features(), via
## make_basis()):
##   current           -- vech(U'AU), K=10                              ->   55 features
##   onesided          -- AU flattened, K=10                            -> 1000 features
##   onesided_rownorm  -- row norms of AU (one scalar per node), K=10   ->  100 features
##   bigK20/30/40      -- vech(U'AU) at K=20/30/40                      -> 210/465/820 features
##   multibasis        -- vech(U'AU), K=10, concatenated across the     ->  165 features
##                         three bases {adjacency, nlaplacian, laplacian}
##
## Every variant is fit with the SAME estimator as the current method -- graph_spectral_ard()
## for regression, graph_spectral_ard_logit() for binary -- so only the feature map changes.
## GUARD: if a variant's feature count q >= n_train for a given split, that (outcome, split,
## variant) is logged as skipped (status="skipped", with a reason) rather than fit as an
## unidentified model. n_train here is ~201-202 for nbk_0v1 and ~101 for age_continuous
## (Schaefer-100, 144 subjects, 30% test, subject-level splits) -- small enough that most of
## the larger variants (onesided, bigK20/30/40) are skipped on every split; see the report.
##
## Outcomes: nbk_0v1 (binary, AUC) and age_continuous (regression, Spearman + RMSE), 10
## subject-level splits each, seed = 900 + split, 30% test -- identical split scheme to
## Codes/analysis/job2_nbk0v1_sch100.R / job3_age_continuous_sch100.R, so results are directly
## comparable to the existing Schaefer-100 tables.
##
## Uses data/osf_schaefer100/master.rds (V=100, p=4,950) throughout, per the task's directive
## to work at this resolution rather than V=229.
##
## Resumable: appends whole splits (all variant + reference rows) to
## results/dense_variants_sch100.csv; a split already fully written (9 rows: 7 variants + ridge
## + kundu) is skipped on rerun. Ownership: this file and its two output CSVs only -- does not
## touch CLAUDE.md, RESUME.md, notebooks/, docs/, Codes/analysis/model_df_comparison.R, or any
## existing results/*.csv.
## =====================================================================================

suppressMessages({
  source("Codes/analysis/multi_resolution_common.R")   # -> outcome_explorer.R -> abide_classification.R
                                                          # -> abide_pipeline.R -> network_bayes_simulation.R
                                                          # brings in: make_basis, vech_features, build_features,
                                                          # graph_spectral_ard[_logit], raw_edge_ridge_woodbury,
                                                          # raw_edge_glmnet_logit, auc_pred, rmse, make_set_m
                                                          # + comparators_cv.R, features_cv.R (build_edge_features)
  source("Codes/benchmarks/competitor_benchmarks.R")    # kundu_latent_nonlinear (binary reference arm)
})

OUT_CSV     <- "results/dense_variants_sch100.csv"
SUMMARY_CSV <- "results/dense_variants_summary.csv"
MASTER_PATH <- "data/osf_schaefer100/master.rds"

## ---- Job 1: feature-map variants -----------------------------------------------------

VARIANT_SPECS <- list(
  current          = list(kind = "twosided",  K = 10, basis = "adjacency"),
  onesided         = list(kind = "onesided",  K = 10, basis = "adjacency"),
  onesided_rownorm = list(kind = "rownorm",   K = 10, basis = "adjacency"),
  bigK20           = list(kind = "twosided",  K = 20, basis = "adjacency"),
  bigK30           = list(kind = "twosided",  K = 30, basis = "adjacency"),
  bigK40           = list(kind = "twosided",  K = 40, basis = "adjacency"),
  multibasis       = list(kind = "multibasis", K = 10, basis = c("adjacency", "nlaplacian", "laplacian"))
)

## Nominal feature count, computable without building the matrix (used for the q>=n_train
## guard so a variant that is certain to be skipped never pays for an eigendecomposition).
nominal_dim <- function(variant, V) {
  sp <- VARIANT_SPECS[[variant]]
  K <- sp$K
  switch(sp$kind,
         twosided   = K * (K + 1) / 2,
         onesided   = V * K,
         rownorm    = V,
         multibasis = length(sp$basis) * K * (K + 1) / 2)
}

## Train-only template basis Q (V x K) for one basis type -- identical construction to
## build_features()/make_basis() in Codes/core/abide_pipeline.R (mean |A| over TRAIN rows).
build_template_Q <- function(conn, train_idx, K, basis) {
  Wbar <- Reduce(`+`, lapply(conn[train_idx], abs)) / length(train_idx)
  diag(Wbar) <- 0
  make_basis(Wbar, K, basis)$Q
}

feat_twosided <- function(conn, Q) {
  K <- ncol(Q); N <- length(conn)
  X <- matrix(0, N, K * (K + 1) / 2)
  for (i in seq_len(N)) X[i, ] <- vech_features(t(Q) %*% conn[[i]] %*% Q)
  X
}

feat_onesided <- function(conn, Q) {
  N <- length(conn)
  X <- matrix(0, N, nrow(Q) * ncol(Q))
  for (i in seq_len(N)) X[i, ] <- as.numeric(conn[[i]] %*% Q)
  X
}

feat_rownorm <- function(conn, Q) {
  N <- length(conn)
  X <- matrix(0, N, nrow(Q))
  for (i in seq_len(N)) {
    L <- conn[[i]] %*% Q
    X[i, ] <- sqrt(rowSums(L^2))
  }
  X
}

## Build the full N x q feature matrix for one variant, given TRAIN row indices.
build_variant_matrix <- function(conn, train_idx, variant) {
  sp <- VARIANT_SPECS[[variant]]
  if (sp$kind == "multibasis") {
    mats <- lapply(sp$basis, function(b) {
      Q <- build_template_Q(conn, train_idx, sp$K, b)
      feat_twosided(conn, Q)
    })
    return(do.call(cbind, mats))
  }
  Q <- build_template_Q(conn, train_idx, sp$K, sp$basis)
  switch(sp$kind,
         twosided = feat_twosided(conn, Q),
         onesided = feat_onesided(conn, Q),
         rownorm  = feat_rownorm(conn, Q))
}

## ---- Kundu reference arm for the CONTINUOUS outcome -----------------------------------
## competitor_benchmarks.R::kundu_latent_nonlinear() clamps its output to [1e-6, 1-1e-6] for
## use as a binary AUC score -- appropriate for nbk_0v1 but not for age (range ~20-86), where
## clamping would destroy the prediction. This is the same node-latent embedding + RBF-kernel-
## ridge fit, with the clamp removed, for use as the regression reference arm.
kundu_regress_pred <- function(conn_train, y_train, conn_test, R = 8, ridge = 1e-2) {
  t0 <- proc.time()[3]
  Wbar <- Reduce(`+`, lapply(conn_train, abs)) / length(conn_train); diag(Wbar) <- 0
  eg <- eigen(Wbar, symmetric = TRUE); ord <- order(abs(eg$values), decreasing = TRUE)
  U <- eg$vectors[, ord[seq_len(R)], drop = FALSE]

  feat <- function(cl) t(sapply(cl, function(A) as.numeric(A %*% U)))
  Ftr <- feat(conn_train); Fte <- feat(conn_test)
  cm <- colMeans(Ftr); cs <- apply(Ftr, 2, sd); cs[cs == 0] <- 1
  Ftr <- sweep(sweep(Ftr, 2, cm, "-"), 2, cs, "/")
  Fte <- sweep(sweep(Fte, 2, cm, "-"), 2, cs, "/")

  sq <- function(A, B) outer(rowSums(A^2), rowSums(B^2), "+") - 2 * A %*% t(B)
  Dtr <- sq(Ftr, Ftr); sigma2 <- median(Dtr[upper.tri(Dtr)]); sigma2 <- ifelse(sigma2 > 0, sigma2, 1)
  Ktr <- exp(-Dtr / (2 * sigma2))
  Kte <- exp(-sq(Fte, Ftr) / (2 * sigma2))
  yc <- y_train - mean(y_train)
  alpha <- solve(Ktr + ridge * diag(nrow(Ktr)), yc)
  pred <- as.numeric(Kte %*% alpha) + mean(y_train)
  list(pred = pred, seconds = proc.time()[3] - t0)
}

## ---- Job 2: benchmark loop -------------------------------------------------------------

blank_row <- function(outcome, split, variant, dim_nom) {
  data.frame(outcome = outcome, split = split, variant = variant,
             dim = dim_nom, status = "skipped", reason = NA_character_,
             rmse = NA_real_, spearman = NA_real_, auc = NA_real_, acc = NA_real_,
             seconds = NA_real_, n_selected = NA_integer_, pred_sd = NA_real_,
             stringsAsFactors = FALSE)
}

eval_one_split <- function(outcome, type, conn, y, tr, te, V) {
  rows <- list()

  ## -- variants --
  for (variant in names(VARIANT_SPECS)) {
    q_nom <- nominal_dim(variant, V)
    n_train <- length(tr)
    if (q_nom >= n_train) {
      r <- blank_row(outcome, NA, variant, q_nom)
      r$reason <- sprintf("q>=n_train (q=%d, n_train=%d)", q_nom, n_train)
      rows[[length(rows) + 1]] <- r
      next
    }
    t0 <- proc.time()[3]
    X <- build_variant_matrix(conn, tr, variant)
    q <- ncol(X)
    if (q >= n_train) {   # re-check against the ACTUAL built dimension (should match q_nom)
      r <- blank_row(outcome, NA, variant, q)
      r$reason <- sprintf("q>=n_train after build (q=%d, n_train=%d)", q, n_train)
      rows[[length(rows) + 1]] <- r
      next
    }
    if (type == "regression") {
      f <- graph_spectral_ard(X[tr, , drop = FALSE], y[tr], X[te, , drop = FALSE])
      secs <- proc.time()[3] - t0
      pred <- f$pred
      r <- data.frame(outcome = outcome, split = NA, variant = variant, dim = q, status = "ok",
                       reason = NA_character_, rmse = rmse(y[te], pred),
                       spearman = suppressWarnings(cor(y[te], pred, method = "spearman")),
                       auc = NA_real_, acc = NA_real_, seconds = secs,
                       n_selected = sum(abs(f$coef) > 1e-6), pred_sd = sd(pred),
                       stringsAsFactors = FALSE)
    } else {
      f <- graph_spectral_ard_logit(X[tr, , drop = FALSE], y[tr], X[te, , drop = FALSE])
      secs <- proc.time()[3] - t0
      prob <- f$prob
      r <- data.frame(outcome = outcome, split = NA, variant = variant, dim = q, status = "ok",
                       reason = NA_character_, rmse = NA_real_, spearman = NA_real_,
                       auc = auc_pred(prob, y[te]), acc = mean((prob > 0.5) == y[te]),
                       seconds = secs, n_selected = f$nsel, pred_sd = sd(prob),
                       stringsAsFactors = FALSE)
    }
    rows[[length(rows) + 1]] <- r
  }

  ## -- reference arms: ridge (raw edge, project-standard) and kundu (node-latent + GP) --
  if (type == "regression") {
    ft <- build_edge_features_cached(conn)
    Xe <- ft$X_edge
    t0 <- proc.time()[3]
    fr <- raw_edge_ridge_woodbury(Xe[tr, , drop = FALSE], y[tr], Xe[te, , drop = FALSE])
    rows[[length(rows) + 1]] <- data.frame(
      outcome = outcome, split = NA, variant = "ridge", dim = ft$p, status = "ok",
      reason = NA_character_, rmse = rmse(y[te], fr$pred),
      spearman = suppressWarnings(cor(y[te], fr$pred, method = "spearman")),
      auc = NA_real_, acc = NA_real_, seconds = fr$seconds, n_selected = NA_integer_,
      pred_sd = sd(fr$pred), stringsAsFactors = FALSE)

    fk <- kundu_regress_pred(conn[tr], y[tr], conn[te], R = 8)
    rows[[length(rows) + 1]] <- data.frame(
      outcome = outcome, split = NA, variant = "kundu", dim = V * 8, status = "ok",
      reason = NA_character_, rmse = rmse(y[te], fk$pred),
      spearman = suppressWarnings(cor(y[te], fk$pred, method = "spearman")),
      auc = NA_real_, acc = NA_real_, seconds = fk$seconds, n_selected = NA_integer_,
      pred_sd = sd(fk$pred), stringsAsFactors = FALSE)
  } else {
    ft <- build_edge_features_cached(conn)
    Xe <- ft$X_edge
    fr <- raw_edge_glmnet_logit(Xe[tr, , drop = FALSE], y[tr], Xe[te, , drop = FALSE], alpha = 0)
    rows[[length(rows) + 1]] <- data.frame(
      outcome = outcome, split = NA, variant = "ridge", dim = ft$p, status = "ok",
      reason = NA_character_, rmse = NA_real_, spearman = NA_real_,
      auc = auc_pred(fr$prob, y[te]), acc = mean((fr$prob > 0.5) == y[te]),
      seconds = fr$seconds, n_selected = fr$nsel, pred_sd = sd(fr$prob), stringsAsFactors = FALSE)

    fk <- kundu_latent_nonlinear(conn[tr], y[tr], conn[te], R = 8)
    rows[[length(rows) + 1]] <- data.frame(
      outcome = outcome, split = NA, variant = "kundu", dim = V * 8, status = "ok",
      reason = NA_character_, rmse = NA_real_, spearman = NA_real_,
      auc = auc_pred(fk$prob, y[te]), acc = mean((fk$prob > 0.5) == y[te]),
      seconds = fk$seconds, n_selected = NA_integer_, pred_sd = sd(fk$prob), stringsAsFactors = FALSE)
  }

  do.call(rbind, rows)
}

## Cache the raw-edge design (depends only on conn, not on split) so it is built once per
## outcome rather than once per split.
.edge_cache <- new.env()
build_edge_features_cached <- function(conn) {
  key <- paste0(length(conn), "_", nrow(conn[[1]]))
  if (is.null(.edge_cache[[key]])) .edge_cache[[key]] <- build_edge_features(conn)
  .edge_cache[[key]]
}

N_SPLITS <- 10
SEED0 <- 900
TEST_FRAC <- 0.30
N_VARIANT_ROWS <- length(VARIANT_SPECS) + 2   # 7 variants + ridge + kundu

run_outcome_dense_variants <- function(m, outcome, type, conds, aggregate_by_subject, y_col,
                                       out_csv = OUT_CSV) {
  set <- make_set_m(m, conds = conds, aggregate_by_subject = aggregate_by_subject)
  conn <- set$conn; meta <- set$meta
  y <- if (outcome == "nbk_0v1") as.integer(meta$COND == "nbk_1back") else as.numeric(meta[[y_col]])
  V <- nrow(conn[[1]])
  subj <- meta$SUBJECT; usub <- unique(subj)
  cat(sprintf("\n[%s] type=%s N=%d rows, %d subjects, V=%d\n", outcome, type, length(conn), length(usub), V))

  done <- integer(0)
  if (file.exists(out_csv)) {
    prev <- read.csv(out_csv, stringsAsFactors = FALSE)
    tab <- table(prev$split[prev$outcome == outcome])
    done <- as.integer(names(tab)[tab >= N_VARIANT_ROWS])
  }

  for (s in seq_len(N_SPLITS)) {
    if (s %in% done) { cat(sprintf("  split %d already complete -- skipping\n", s)); next }
    set.seed(SEED0 + s)
    te_s <- sample(usub, max(1, round(TEST_FRAC * length(usub))))
    te <- which(subj %in% te_s); tr <- which(!(subj %in% te_s))
    if (type != "regression" && (length(unique(y[tr])) < 2 || length(unique(y[te])) < 2)) {
      cat(sprintf("  split %d: degenerate label split -- skipping\n", s)); next
    }
    r <- eval_one_split(outcome, type, conn, y, tr, te, V)
    r$split <- s
    prev <- if (file.exists(out_csv)) read.csv(out_csv, stringsAsFactors = FALSE) else NULL
    dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
    write.csv(rbind(prev, r), out_csv, row.names = FALSE)
    ok <- r[r$status == "ok", ]
    prim <- if (type == "regression") sprintf("%s=%.3f", ok$variant, ok$spearman) else sprintf("%s=%.3f", ok$variant, ok$auc)
    cat(sprintf("  split %d: %s\n", s, paste(prim, collapse = " ")))
    flush.console()
  }
  invisible(TRUE)
}

## ---- summary: means + paired Wilcoxon of every variant vs `current` and vs `kundu` ----
## (a general utility, parameterized by out_csv/summary_csv so other drivers -- e.g. the
## no-guard rerun in dense_variants_nogard.R -- can reuse it on their own combined data.)

summarise_dense_variants <- function(out_csv = OUT_CSV, summary_csv = SUMMARY_CSV) {
  d <- read.csv(out_csv, stringsAsFactors = FALSE)
  outcomes <- unique(d$outcome)
  variants <- unique(d$variant)
  rows <- list()

  for (oc in outcomes) {
    type <- if (oc == "age_continuous") "regression" else "binary"
    primary <- if (type == "regression") "spearman" else "auc"
    dsub <- d[d$outcome == oc, ]

    for (v in variants) {
      dv <- dsub[dsub$variant == v, ]
      n_ok <- sum(dv$status == "ok")
      n_skip <- sum(dv$status == "skipped")
      mean_dim <- if (n_ok > 0) mean(dv$dim[dv$status == "ok"]) else if (nrow(dv) > 0) dv$dim[1] else NA
      row <- data.frame(
        outcome = oc, variant = v, dim = mean_dim,
        n_splits_ok = n_ok, n_splits_skipped = n_skip,
        mean_rmse = if (type == "regression") mean(dv$rmse, na.rm = TRUE) else NA_real_,
        mean_spearman = if (type == "regression") mean(dv$spearman, na.rm = TRUE) else NA_real_,
        mean_auc = if (type == "binary") mean(dv$auc, na.rm = TRUE) else NA_real_,
        mean_acc = if (type == "binary") mean(dv$acc, na.rm = TRUE) else NA_real_,
        mean_seconds = mean(dv$seconds, na.rm = TRUE),
        mean_n_selected = mean(dv$n_selected, na.rm = TRUE),
        mean_pred_sd = mean(dv$pred_sd, na.rm = TRUE),
        stringsAsFactors = FALSE
      )

      ## paired tests vs current and vs kundu, matched by split, both arms status=="ok"
      pair_test <- function(ref) {
        dref <- dsub[dsub$variant == ref & dsub$status == "ok", c("split", primary)]
        dvv  <- dv[dv$status == "ok", c("split", primary)]
        names(dref)[2] <- "ref_val"; names(dvv)[2] <- "v_val"
        mrg <- merge(dvv, dref, by = "split")
        if (nrow(mrg) < 3 || v == ref) return(c(n_paired = nrow(mrg), mean_diff = NA_real_, wilcoxon_p = NA_real_))
        wt <- try(wilcox.test(mrg$v_val, mrg$ref_val, paired = TRUE), silent = TRUE)
        p <- if (inherits(wt, "try-error")) NA_real_ else wt$p.value
        c(n_paired = nrow(mrg), mean_diff = mean(mrg$v_val - mrg$ref_val), wilcoxon_p = p)
      }
      vs_current <- pair_test("current")
      vs_kundu <- pair_test("kundu")
      row$n_paired_vs_current <- vs_current["n_paired"]; row$diff_vs_current <- vs_current["mean_diff"]; row$wilcoxon_p_vs_current <- vs_current["wilcoxon_p"]
      row$n_paired_vs_kundu <- vs_kundu["n_paired"]; row$diff_vs_kundu <- vs_kundu["mean_diff"]; row$wilcoxon_p_vs_kundu <- vs_kundu["wilcoxon_p"]
      rows[[length(rows) + 1]] <- row
    }
  }
  summ <- do.call(rbind, rows)
  write.csv(summ, summary_csv, row.names = FALSE)
  cat("\nWrote", summary_csv, "\n")
  invisible(summ)
}

## ---- driver (only runs when this file is executed directly, e.g. `Rscript
## Codes/analysis/dense_variants.R` -- sourcing it from another script, as
## dense_variants_nogard.R does, only defines the functions above) ---------------------

if (sys.nframe() == 0) {
  m <- readRDS(MASTER_PATH)
  cat(sprintf("Loaded %s: %d connectomes, V=%d, %d subjects\n",
              MASTER_PATH, length(m$conn), m$V, length(unique(m$subjects$SUBJECT))))

  run_outcome_dense_variants(m, outcome = "nbk_0v1", type = "binary",
                             conds = c("nbk_0back", "nbk_1back"), aggregate_by_subject = FALSE,
                             y_col = NULL)

  run_outcome_dense_variants(m, outcome = "age_continuous", type = "regression",
                             conds = NULL, aggregate_by_subject = TRUE, y_col = "AGE")

  summ <- summarise_dense_variants()
  cat("\n==== Summary ====\n")
  print(summ, row.names = FALSE)
}
