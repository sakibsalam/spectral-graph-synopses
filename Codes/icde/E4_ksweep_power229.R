## =====================================================================================
## E4 (optional): synopsis-size / denser-variant sweep at V = 229, mirroring
## Codes/analysis/dense_variants.R (which ran at Schaefer-100). Same variants, same 10 splits.
##   Rscript Codes/icde/E4_ksweep_power229.R [--smoke] [--no-onesided]
## WARNING: one-sided ARD (V*K = 2290 features) costs minutes to tens of minutes per fit at V = 229
## (ARD solves a q x q system per iteration, q^3 scaling). Use --no-onesided to skip it. Each fit is
## capped at 900 s (the cap is checked between interruptible steps, so a single long solve can overrun).
## Variants: current (K=10, 55), onesided (2290), onesided_rownorm (229), bigK20/30/40 (210/465/820),
## multibasis (165), plus references ridge (edge, as in the existing table) and kundu.
## Output columns follow results/dense_variants_sch100_combined.csv (minus acc / n_selected).
## =====================================================================================
source("Codes/icde/icde_common.R")
smoke <- has_flag("--smoke"); skip_one <- has_flag("--no-onesided")
n_splits <- if (smoke) 1 else 10
CAP <- 900
SPEC <- list(current = list(kind = "two", K = 10), onesided = list(kind = "one", K = 10),
             onesided_rownorm = list(kind = "row", K = 10), bigK20 = list(kind = "two", K = 20),
             bigK30 = list(kind = "two", K = 30), bigK40 = list(kind = "two", K = 40),
             multibasis = list(kind = "multi", K = 10))
if (smoke) SPEC <- SPEC[c("current", "bigK20")]
if (skip_one) SPEC <- SPEC[names(SPEC) != "onesided"]
build_variant <- function(conn, tr, sp) switch(sp$kind,
  two = feat_twosided(conn, template_Q(conn, tr, sp$K)),
  one = feat_onesided(conn, template_Q(conn, tr, sp$K)),
  row = feat_rownorm(conn, template_Q(conn, tr, sp$K)),
  multi = do.call(cbind, lapply(c("adjacency", "nlaplacian", "laplacian"),
                                function(b) feat_twosided(conn, template_Q(conn, tr, sp$K, b)))))

out <- out_path("ksweep_power229.csv", smoke)
done <- if (file.exists(out)) unique(with(read.csv(out), paste(outcome, split))) else character(0)
for (oc in c("nbk_0v1", "age_continuous")) {
  d <- load_set("power229", oc); conn <- d$conn; y <- d$y; subj <- d$subj
  Xe <- build_edge_features(conn)$X_edge
  for (s in seq_len(n_splits)) {
    if (paste(oc, s) %in% done) next
    sp <- make_split(subj, s); tr <- sp$tr; te <- sp$te
    rows <- list()
    rec <- function(variant, dim, status, val, secs, reason = NA_character_) {
      rows[[length(rows) + 1]] <<- data.frame(outcome = oc, split = s, variant = variant, dim = dim, status = status,
        reason = reason, rmse = if (d$type == "regression" && !is.na(val[1])) rmse(y[te], val) else NA_real_,
        spearman = if (d$type == "regression" && !is.na(val[1])) spearman(y[te], val) else NA_real_,
        auc = if (d$type == "binary" && !is.na(val[1])) auc_pred(val, y[te]) else NA_real_, seconds = secs)
    }
    run <- function(variant, dim, fn) {            # fn returns test predictions (prob or pred)
      t0 <- proc.time()[3]; setTimeLimit(elapsed = CAP, transient = FALSE)
      r <- tryCatch(fn(), error = function(e) conditionMessage(e))
      setTimeLimit(elapsed = Inf, cpu = Inf, transient = FALSE)
      if (is.character(r)) rec(variant, dim, if (grepl("time limit", r)) "timeout" else "error", NA, proc.time()[3] - t0, substr(r, 1, 100))
      else rec(variant, dim, "ok", r, proc.time()[3] - t0)
      cat(sprintf("  %s split %d %-17s dim=%5d %s\n", oc, s, variant, dim, tail(rows, 1)[[1]]$status)); flush.console()
    }
    for (v in names(SPEC)) {
      X <- tryCatch(build_variant(conn, tr, SPEC[[v]]), error = function(e) NULL)
      if (is.null(X)) { rec(v, NA, "error", NA, NA, "feature build failed"); next }
      run(v, ncol(X), function() if (d$type == "binary") graph_spectral_ard_logit(X[tr, ], y[tr], X[te, ])$prob
                                 else graph_spectral_ard(X[tr, ], y[tr], X[te, ])$pred)
    }
    run("ridge", ncol(Xe), function() if (d$type == "binary") raw_edge_glmnet_logit(Xe[tr, ], y[tr], Xe[te, ], alpha = 0)$prob
                                       else raw_edge_ridge_woodbury(Xe[tr, ], y[tr], Xe[te, ])$pred)   # alpha = 50, as in dense_variants.R
    run("kundu", d$V * 8, function() if (d$type == "binary") kundu_latent_nonlinear(conn[tr], y[tr], conn[te])$prob
                                     else kundu_regress_pred(conn[tr], y[tr], conn[te])$pred)
    append_csv(do.call(rbind, rows), out)
  }
}
cat("Wrote", out, "\n")
