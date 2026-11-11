## =====================================================================================
## Outcome explorer for the OSF Power-229 dataset (Rieck 2021 / BMGC data).
##
## Purpose: systematically enumerate and BENCHMARK every scalar outcome derivable from this
## dataset, on a COMMON metric surface, so we can identify which outcomes are (a) predictable
## at all and (b) favourable to the graph-spectral representation.
##
## Statistical design (applied identically to every outcome):
##   * SPLIT UNIT = SUBJECT. A subject contributes up to 7 condition-connectomes; all of a
##     subject's rows go to train OR test. This prevents subject-identity leakage, which would
##     otherwise inflate every within-subject (condition-decoding) outcome.
##   * PAIRED SPLITS. All methods see the identical splits, so method differences can be tested
##     with a paired test rather than compared as independent means.
##   * TRAIN-ONLY basis + standardization (build_features derives the template graph from the
##     training connectomes only).
##   * COMMON SIGNAL SCALE so heterogeneous outcomes are comparable:
##       regression      -> Spearman rho (rank correlation of predicted vs observed)
##       binary/multiclass -> Somers' D = 2*(AUC - 0.5)
##     Both are 0 under no signal and 1 under perfect ordering.
##
## Outcome types supported: "regression", "binary", "multiclass" (one-vs-rest).
## Results append to results/outcome_explorer.csv and are RESUMABLE per outcome name.
## =====================================================================================

suppressMessages({
  source("Codes/core/abide_classification.R")   # build_features, ARD logit, glmnet logit, auc_pred
  source("Codes/core/hybrid_model.R")           # graph_spectral_plus_sparse[_logit]
})

MASTER_RDS <- "data/osf_master/master.rds"
.master <- NULL
load_master <- function() {
  if (is.null(.master)) .master <<- readRDS(MASTER_RDS)
  .master
}

## Build an analysis set: optionally restrict to conditions, optionally average each subject's
## connectomes into one matrix (subject-level outcomes such as age).
make_set <- function(conds = NULL, aggregate_by_subject = FALSE) {
  m <- load_master(); s <- m$subjects
  keep <- if (is.null(conds)) seq_len(nrow(s)) else which(s$COND %in% conds)
  conn <- m$conn[keep]; meta <- s[keep, , drop = FALSE]
  if (aggregate_by_subject) {
    us <- unique(meta$SUBJECT)
    conn2 <- lapply(us, function(u) Reduce(`+`, conn[meta$SUBJECT == u]) / sum(meta$SUBJECT == u))
    meta2 <- meta[!duplicated(meta$SUBJECT), ]
    meta2 <- meta2[match(us, meta2$SUBJECT), , drop = FALSE]
    conn <- conn2; meta <- meta2
  }
  list(conn = conn, meta = meta)
}

## ---- metrics ------------------------------------------------------------------------
spearman <- function(a, b) suppressWarnings(cor(a, b, method = "spearman"))
r2_test  <- function(y, yhat) 1 - sum((y - yhat)^2) / sum((y - mean(y))^2)
brier    <- function(p, y) mean((p - y)^2)

## ---- one split of one outcome -------------------------------------------------------
## Returns a long data.frame: one row per (method, metric).
eval_split <- function(conn, y, groups, type, tr, te, K = 10, basis = "adjacency",
                       methods = c("gs_ard", "hybrid", "ridge", "lasso")) {
  ft <- build_features(conn, tr, K, basis = basis)
  Xs <- ft$X_spec; Xe <- ft$X_edge
  out <- list()
  add <- function(method, metric, value)
    out[[length(out) + 1]] <<- data.frame(method = method, metric = metric, value = value)

  if (type == "regression") {
    preds <- list()
    if ("gs_ard" %in% methods) preds$gs_ard <- graph_spectral_ard(Xs[tr, ], y[tr], Xs[te, ])$pred
    if ("hybrid" %in% methods) preds$hybrid <- graph_spectral_plus_sparse(Xs[tr, ], Xe[tr, ], y[tr], Xs[te, ], Xe[te, ])$pred
    if ("ridge"  %in% methods) preds$ridge  <- raw_edge_ridge_woodbury(Xe[tr, ], y[tr], Xe[te, ])$pred
    if ("lasso"  %in% methods) { l <- raw_edge_lasso(Xe[tr, ], y[tr], Xe[te, ]); if (!is.null(l)) preds$lasso <- l$pred }
    preds$baseline <- rep(mean(y[tr]), length(te))
    for (nm in names(preds)) {
      p <- preds[[nm]]
      add(nm, "rmse", rmse(y[te], p)); add(nm, "r2", r2_test(y[te], p))
      add(nm, "spearman", spearman(y[te], p)); add(nm, "signal", max(0, spearman(y[te], p)))
    }

  } else if (type == "binary") {
    probs <- list()
    if ("gs_ard" %in% methods) probs$gs_ard <- graph_spectral_ard_logit(Xs[tr, ], y[tr], Xs[te, ])$prob
    if ("hybrid" %in% methods) probs$hybrid <- graph_spectral_plus_sparse_logit(Xs[tr, ], Xe[tr, ], y[tr], Xs[te, ], Xe[te, ])$prob
    if ("ridge"  %in% methods) probs$ridge  <- raw_edge_glmnet_logit(Xe[tr, ], y[tr], Xe[te, ], alpha = 0)$prob
    if ("lasso"  %in% methods) probs$lasso  <- raw_edge_glmnet_logit(Xe[tr, ], y[tr], Xe[te, ], alpha = 1)$prob
    for (nm in names(probs)) {
      p <- probs[[nm]]; a <- auc_pred(p, y[te])
      add(nm, "auc", a); add(nm, "acc", mean((p > 0.5) == y[te]))
      add(nm, "brier", brier(p, y[te])); add(nm, "signal", max(0, 2 * (a - 0.5)))
    }

  } else if (type == "multiclass") {
    lev <- sort(unique(y)); K_ <- length(lev)
    scores <- list()
    for (nm in methods) scores[[nm]] <- matrix(NA_real_, length(te), K_)
    for (k in seq_along(lev)) {                       # one-vs-rest
      yk <- as.integer(y == lev[k])
      if (length(unique(yk[tr])) < 2) next
      if ("gs_ard" %in% methods) scores$gs_ard[, k] <- graph_spectral_ard_logit(Xs[tr, ], yk[tr], Xs[te, ])$prob
      if ("hybrid" %in% methods) scores$hybrid[, k] <- graph_spectral_plus_sparse_logit(Xs[tr, ], Xe[tr, ], yk[tr], Xs[te, ], Xe[te, ])$prob
      if ("ridge"  %in% methods) scores$ridge[, k]  <- raw_edge_glmnet_logit(Xe[tr, ], yk[tr], Xe[te, ], alpha = 0)$prob
      if ("lasso"  %in% methods) scores$lasso[, k]  <- raw_edge_glmnet_logit(Xe[tr, ], yk[tr], Xe[te, ], alpha = 1)$prob
    }
    for (nm in methods) {
      S <- scores[[nm]]; if (all(is.na(S))) next
      aucs <- sapply(seq_along(lev), function(k) auc_pred(S[, k], as.integer(y[te] == lev[k])))
      macro <- mean(aucs, na.rm = TRUE)
      pred_lab <- lev[apply(S, 1, which.max)]
      add(nm, "macro_auc", macro); add(nm, "acc", mean(pred_lab == y[te]))
      add(nm, "signal", max(0, 2 * (macro - 0.5)))
      ## Ordinal sensitivity: rank correlation of the score-weighted expected class with the truth.
      ## Only defined when the class labels are themselves numeric/ordered (e.g. n-back LOAD 0/1/2,
      ## age quartile 1..4); for nominal labels (DOMAIN, COND) there is no ordering, so it is skipped.
      ## NOTE: run_outcome() coerces y to character for type="multiclass", so `lev` is character —
      ## coerce explicitly here. (Doing `S %*% lev` on a character vector errors, and because
      ## eval_split is wrapped in try() that error would silently drop every multiclass split.)
      lev_num <- suppressWarnings(as.numeric(lev))
      if (!any(is.na(lev_num))) {
        exp_lab <- as.numeric(S %*% lev_num) / pmax(rowSums(S), 1e-9)
        add(nm, "spearman", spearman(suppressWarnings(as.numeric(y[te])), exp_lab))
      }
    }
    add("baseline", "acc", max(table(y[tr])) / length(tr)); add("baseline", "signal", 0)
  }
  do.call(rbind, out)
}

## ---- run one outcome over repeated splits (resumable) -------------------------------
run_outcome <- function(name, conds = NULL, aggregate_by_subject = FALSE, y_col = NULL, type,
                        n_splits = 10, test_frac = 0.30, K = 10, basis = "adjacency",
                        methods = c("gs_ard", "hybrid", "ridge", "lasso"), seed0 = 900,
                        out_csv = "results/outcome_explorer.csv", verbose = TRUE,
                        y_override_fn = NULL) {
  ## y_override_fn(meta): optional -- for outcomes whose label is not a stored master column
  ## (e.g. a binary "is COND == the second of two selected conditions" contrast). Takes the
  ## conds-filtered (pre-aggregate) meta data.frame and returns a vector aligned to its rows.
  ## Backward compatible: when NULL (default, every pre-existing call site), behavior is
  ## identical to before this parameter existed.
  set <- make_set(conds, aggregate_by_subject)
  meta <- set$meta; conn <- set$conn
  y <- if (!is.null(y_override_fn)) y_override_fn(meta) else meta[[y_col]]
  ok <- which(!is.na(y) & y != "")
  conn <- conn[ok]; meta <- meta[ok, , drop = FALSE]; y <- y[ok]
  if (type != "multiclass") y <- as.numeric(y) else y <- as.character(y)
  subj <- meta$SUBJECT; usub <- unique(subj)
  if (verbose) cat(sprintf("\n[%s] type=%s  N=%d rows, %d subjects, conds=%s\n", name, type,
                           length(conn), length(usub), if (is.null(conds)) "all" else paste(conds, collapse = "+")))

  done <- integer(0)
  if (file.exists(out_csv)) { prev <- read.csv(out_csv)
    done <- unique(prev$split[prev$outcome == name]) }
  acc <- list()
  for (s in seq_len(n_splits)) {
    if (s %in% done) next
    set.seed(seed0 + s)
    te_s <- sample(usub, max(1, round(test_frac * length(usub))))    # SUBJECT-level split
    te <- which(subj %in% te_s); tr <- which(!(subj %in% te_s))
    if (type != "regression" && (length(unique(y[tr])) < 2 || length(unique(y[te])) < 2)) next
    r <- try(eval_split(conn, y, subj, type, tr, te, K = K, basis = basis, methods = methods), silent = TRUE)
    if (inherits(r, "try-error")) { cat("   split", s, "failed:", as.character(r), "\n"); next }
    r$outcome <- name; r$split <- s; r$type <- type; r$n <- length(conn); r$n_subj <- length(usub)
    acc[[length(acc) + 1]] <- r
    dir.create("results", showWarnings = FALSE)
    prev <- if (file.exists(out_csv)) read.csv(out_csv) else NULL
    write.csv(rbind(prev, do.call(rbind, acc[length(acc)])), out_csv, row.names = FALSE)
    if (verbose) {
      sig <- r[r$metric == "signal", ]
      cat(sprintf("   split %d: %s\n", s, paste(sprintf("%s=%.2f", sig$method, sig$value), collapse = " ")))
      flush.console()
    }
  }
  invisible(TRUE)
}

## summary across everything computed so far
summarise_outcomes <- function(out_csv = "results/outcome_explorer.csv") {
  d <- read.csv(out_csv)
  agg <- aggregate(value ~ outcome + type + method + metric, d, mean)
  sig <- agg[agg$metric == "signal", ]
  w <- reshape(sig[, c("outcome", "type", "method", "value")], idvar = c("outcome", "type"),
               timevar = "method", direction = "wide")
  names(w) <- sub("value.", "", names(w))
  w[order(-w$hybrid), ]
}
