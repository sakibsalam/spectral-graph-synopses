## =====================================================================================
## Figures for the ICDE paper (base R graphics, Okabe-Ito colorblind-safe palette).
##   Rscript Codes/icde/make_icde_figures.R [outdir]      # default outdir = results/icde/figs/
## Reads results/icde/*.csv plus the existing results/model_df_comparison.csv and
## results/dense_variants_sch100_combined.csv. A figure whose input is missing is skipped.
## Fig scaling: measured points solid; full-fit lines extrapolated from measured per-iteration
## cost (iterations x seconds_per_iter) are DASHED and labeled "projected". They are plots only,
## never a measured fit time.
## =====================================================================================
a <- grep("^--", commandArgs(trailingOnly = TRUE), invert = TRUE, value = TRUE)
outdir <- if (length(a)) a[1] else "results/icde/figs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
pal <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#999999")
med <- function(x) median(x, na.rm = TRUE)

## ---------------- Fig scaling: seconds vs V (log-log) -----------------------------------
f <- "results/icde/scaling.csv"
if (file.exists(f)) {
  s <- read.csv(f); s <- s[s$status == "ok", ]
  series <- list(                                  # label, design, learner, column, color, pch
    list("two-sided ARD (K=10)",      "two_sided", "ard_logit_K10",      "seconds",          pal[4], 16),
    list("one-sided ridge (R=10)",    "one_sided", "ridge_dual_cv_R10",  "seconds",          pal[3], 17),
    list("edge ridge (dual)",         "edge_store", "ridge_dual_cv",     "seconds",          pal[2], 15),
    list("edge lasso",                "edge_store", "lasso_glmnet",      "seconds",          pal[7], 18),
    list("edge SSVS / iter",          "edge_store", "ssvs_iter",         "seconds_per_iter", pal[5], 1),
    list("edge ARD / iter",           "edge_store", "ard_iter",          "seconds_per_iter", pal[6], 2),
    list("BMGC / iter",               "edge_store", "bmgc_iter",         "seconds_per_iter", pal[1], 0))
  proj <- list(list("edge ARD, 500 it (projected)", "ard_iter", 500, pal[6]),
               list("BMGC, 2000 it (projected)", "bmgc_iter", 2000, pal[1]),
               list("SSVS, 1200 it (projected)", "ssvs_iter", 1200, pal[5]))
  get <- function(des, lrn, col) { z <- s[s$design == des & s$learner == lrn & !is.na(s[[col]]), ]
    if (!nrow(z)) return(NULL); ag <- aggregate(z[[col]], list(V = z$V), med); names(ag)[2] <- "y"; ag[order(ag$V), ] }
  dat <- lapply(series, function(x) get(x[[2]], x[[3]], x[[4]]))
  pdat <- lapply(proj, function(x) { g <- get("edge_store", x[[2]], "seconds_per_iter"); if (!is.null(g)) g$y <- g$y * x[[3]]; g })
  allv <- unlist(lapply(c(dat, pdat), function(g) if (is.null(g)) NULL else g$y)); allx <- unlist(lapply(dat, function(g) g$V))
  if (length(allv)) {
    pdf(file.path(outdir, "fig_scaling.pdf"), width = 3.4, height = 3.2, pointsize = 7)
    par(mar = c(3.6, 3.8, 0.6, 0.6), mgp = c(2.2, 0.7, 0))
    plot(range(allx), range(allv[allv > 0]), type = "n", log = "xy", xlab = "number of vertices V",
         ylab = "seconds (fit, or per iteration)")
    for (i in seq_along(series)) if (!is.null(dat[[i]])) {
      lines(dat[[i]]$V, dat[[i]]$y, col = series[[i]][[5]], lty = 1); points(dat[[i]]$V, dat[[i]]$y, col = series[[i]][[5]], pch = series[[i]][[6]], cex = 0.8) }
    for (i in seq_along(proj)) if (!is.null(pdat[[i]])) lines(pdat[[i]]$V, pdat[[i]]$y, col = proj[[i]][[4]], lty = 2)
    legend("topleft", bty = "n", cex = 0.55, y.intersp = 0.9,
           legend = c(sapply(series, `[[`, 1), sapply(proj, `[[`, 1)),
           col = c(sapply(series, `[[`, 5), sapply(proj, `[[`, 4)),
           pch = c(sapply(series, `[[`, 6), rep(NA, length(proj))), lty = c(rep(1, length(series)), rep(2, length(proj))))
    dev.off(); cat("wrote fig_scaling.pdf\n")
  }
} else cat("skip fig_scaling: results/icde/scaling.csv missing\n")

## ---------------- Fig df: effective df vs nominal dimension, two resolutions, two outcomes -------
## nbk_0v1: Schaefer-100 from results/model_df_comparison.csv, Power-229 from results/icde/df_power229.csv.
## age_continuous: both from results/icde/df_power229.csv and results/icde/df_sch100.csv (E2, E2b).
f2 <- "results/icde/df_power229.csv"; f1 <- "results/model_df_comparison.csv"; f3 <- "results/icde/df_sch100.csv"
if (file.exists(f2)) {
  rd <- function(f, at) { z <- read.csv(f); data.frame(outcome = z$outcome, atlas = at, method = sub("^gs_ard$", "gs_ard_logit", z$method),
                                                      nominal = z$nominal_dim, eff = z$eff_df) }
  d <- rd(f2, "Power-229")
  if (file.exists(f3)) d <- rbind(d, rd(f3, "Schaefer-100"))
  if (file.exists(f1)) { sch <- read.csv(f1)
    d <- rbind(d, data.frame(outcome = "nbk_0v1", atlas = "Schaefer-100", method = sch$model, nominal = sch$nominal_params, eff = sch$effective_df)) }
  keep <- c("gs_ard_logit", "raw_edge_ridge", "kundu_kernel_ridge"); d <- d[d$method %in% keep, ]
  col <- c(gs_ard_logit = pal[4], raw_edge_ridge = pal[2], kundu_kernel_ridge = pal[5])
  lab <- c(gs_ard_logit = "two-sided ARD", raw_edge_ridge = "edge ridge (CV)", kundu_kernel_ridge = "kernel ridge (Kundu re-impl.)")
  pdf(file.path(outdir, "fig_df.pdf"), width = 6.8, height = 3.0, pointsize = 7)
  par(mfrow = c(1, 2), mar = c(3.6, 3.8, 1.6, 0.6), mgp = c(2.2, 0.7, 0))
  for (oc in c("nbk_0v1", "age_continuous")) {
    z <- d[d$outcome == oc, ]; if (!nrow(z)) { plot.new(); next }
    ag <- do.call(rbind, lapply(split(z, list(z$method, z$atlas), drop = TRUE), function(w)
      data.frame(method = w$method[1], atlas = w$atlas[1], nominal = w$nominal[1], m = mean(w$eff), sd = sd(w$eff))))
    sdz <- ifelse(is.na(ag$sd), 0, ag$sd)
    plot(range(ag$nominal) * c(0.6, 1.8), c(max(0.5, min(ag$m - sdz) * 0.7), max(ag$m + sdz) * 1.5),
         type = "n", log = "xy", xlab = "nominal dimension", ylab = "effective degrees of freedom", main = oc, cex.main = 0.9)
    abline(0, 1, lty = 3, col = "grey60")            # df = nominal dimension
    for (i in seq_len(nrow(ag))) { pch <- if (ag$atlas[i] == "Power-229") 17 else 16
      if (sdz[i] > 0) arrows(ag$nominal[i], max(ag$m[i] - sdz[i], 0.3), ag$nominal[i], ag$m[i] + sdz[i], angle = 90, code = 3, length = 0.03, col = col[ag$method[i]])
      points(ag$nominal[i], ag$m[i], pch = pch, col = col[ag$method[i]], cex = 1.1) }
    if (oc == "nbk_0v1") legend("topleft", bty = "n", cex = 0.65, legend = c(lab[names(col)], "Schaefer-100", "Power-229", "df = dimension"),
           col = c(col, "black", "black", "grey60"), pch = c(rep(15, 3), 16, 17, NA), lty = c(rep(NA, 5), 3))
  }
  dev.off(); cat("wrote fig_df.pdf\n")
} else cat("skip fig_df: need results/icde/df_power229.csv\n")

## ---------------- Fig kbyn: accuracy vs training fraction ---------------------------------
f <- "results/icde/k_by_n.csv"
if (file.exists(f)) {
  k <- read.csv(f); k$cfg <- ifelse(k$design == "two_sided_ard", paste0("2-sided K=", k$K_or_R),
                              ifelse(k$design == "one_sided_ridge", paste0("1-sided R=", k$K_or_R),
                              ifelse(k$design == "edge_ridge_cv", "edge ridge (CV)", "edge ridge (logistic)")))
  cfgs <- sort(unique(k$cfg)); nk <- sum(grepl("^2-sided", cfgs))
  cols <- setNames(rep(pal[8], length(cfgs)), cfgs)
  cols[grepl("^2-sided", cfgs)] <- colorRampPalette(c("#B7E3CF", "#006D4A"))(nk)
  cols[grepl("^1-sided", cfgs)] <- c(pal[3], pal[5])[seq_len(sum(grepl("^1-sided", cfgs)))]
  cols["edge ridge (CV)"] <- pal[2]; cols["edge ridge (logistic)"] <- pal[6]
  lty <- setNames(ifelse(grepl("^1-sided", cfgs), 2, ifelse(grepl("^edge", cfgs), 4, 1)), cfgs)
  panels <- unique(k[, c("atlas", "outcome")])
  pdf(file.path(outdir, "fig_kbyn.pdf"), width = 7, height = 5.2, pointsize = 8)
  par(mfrow = c(ceiling(nrow(panels) / 2), 2), mar = c(3.4, 3.6, 1.6, 0.6), mgp = c(2.1, 0.7, 0))
  for (i in seq_len(nrow(panels))) {
    z <- k[k$atlas == panels$atlas[i] & k$outcome == panels$outcome[i], ]
    ag <- aggregate(metric ~ cfg + train_frac, z, mean)
    plot(range(ag$train_frac), range(ag$metric, na.rm = TRUE), type = "n", xlab = "training fraction of subjects",
         ylab = paste("mean", z$metric_name[1]), main = paste(panels$atlas[i], panels$outcome[i]), cex.main = 0.9)
    for (cf in unique(ag$cfg)) { g <- ag[ag$cfg == cf, ]; g <- g[order(g$train_frac), ]
      lines(g$train_frac, g$metric, col = cols[cf], lty = lty[cf], lwd = 1.3); points(g$train_frac, g$metric, col = cols[cf], pch = 16, cex = 0.6) }
    if (i == 1) legend("bottomright", bty = "n", cex = 0.6, legend = names(cols), col = cols, lty = lty[names(cols)], lwd = 1.3)
  }
  dev.off(); cat("wrote fig_kbyn.pdf\n")
  ## best design per cell (text summary used for the guidance table)
  agc <- aggregate(metric ~ atlas + outcome + train_frac + cfg, k, mean)
  best <- do.call(rbind, lapply(split(agc, list(agc$atlas, agc$outcome, agc$train_frac), drop = TRUE),
                                function(z) z[which.max(z$metric), ]))
  write.csv(best, file.path(outdir, "kbyn_best_per_cell.csv"), row.names = FALSE)
} else cat("skip fig_kbyn: results/icde/k_by_n.csv missing\n")

## ---------------- Fig ksweep: accuracy vs synopsis dimension -----------------------------
f1 <- "results/dense_variants_sch100_combined.csv"; f2 <- "results/icde/ksweep_power229.csv"
srcs <- list(`Schaefer-100` = f1, `Power-229` = f2); srcs <- srcs[sapply(srcs, file.exists)]
if (length(srcs)) {
  pdf(file.path(outdir, "fig_ksweep.pdf"), width = 7, height = 3.0, pointsize = 8)
  par(mfrow = c(1, 2), mar = c(3.6, 3.8, 1.6, 0.6), mgp = c(2.2, 0.7, 0))
  for (oc in c("nbk_0v1", "age_continuous")) {
    mcol <- if (oc == "nbk_0v1") "auc" else "spearman"; pts <- list()
    for (nm in names(srcs)) { z <- read.csv(srcs[[nm]]); z <- z[z$outcome == oc & z$status == "ok" & !is.na(z[[mcol]]), ]
      if (nrow(z)) { ag <- aggregate(list(m = z[[mcol]], dim = z$dim), list(variant = z$variant), mean); ag$atlas <- nm; pts[[nm]] <- ag } }
    if (!length(pts)) { plot.new(); next }
    P <- do.call(rbind, pts)
    plot(P$dim, P$m, type = "n", log = "x", xlab = "feature dimension", ylab = paste("mean", mcol), main = oc, cex.main = 0.9)
    for (nm in names(pts)) { g <- pts[[nm]]; pch <- if (nm == "Power-229") 17 else 16
      points(g$dim, g$m, pch = pch, col = ifelse(grepl("onesided", g$variant), pal[6], ifelse(g$variant %in% c("ridge", "kundu"), pal[2], pal[4])))
      text(g$dim, g$m, g$variant, pos = 3, cex = 0.5, offset = 0.3) }
    if (oc == "nbk_0v1") legend("bottomright", bty = "n", cex = 0.6, pch = c(16, 17, 15, 15, 15), col = c("black", "black", pal[4], pal[6], pal[2]),
                                legend = c("Schaefer-100", "Power-229", "two-sided / multibasis", "one-sided", "reference"))
  }
  dev.off(); cat("wrote fig_ksweep.pdf\n")
} else cat("skip fig_ksweep: no input CSV\n")
