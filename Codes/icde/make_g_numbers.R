## =====================================================================================
## Paper numbers for the generic-reduction (G1), external-collection (G2), collection-size
## timing (G3) and grouped-fold (G5) experiments, and corrected resampled-t intervals for the
## paired comparisons of the training-size grid (E3).
##   Rscript Codes/icde/make_g_numbers.R <paper_dir>
## Writes <paper_dir>/numbers_g.tex, tab_generic.tex, tab_external.tex and data/g_*.dat.
## Every macro is computed here from results/icde/*.csv; nothing is typed by hand.
## Paired summaries over the 10 subject-level splits: mean difference, number of positive
## splits, exact two-sided Wilcoxon signed-rank p, and the 95% interval of the corrected
## resampled t statistic of Nadeau and Bengio (variance inflated by 1/J + n_test/n_train).
## =====================================================================================
args <- commandArgs(trailingOnly = TRUE); paper <- if (length(args)) args[1] else "paper"
R <- "results/icde"; dir.create(file.path(paper, "data"), showWarnings = FALSE, recursive = TRUE)
mac <- character(0)
M <- function(name, val) mac[[length(mac) + 1]] <<- sprintf("\\newcommand{\\%s}{%s}", name, val)
f3 <- function(x) formatC(x, format = "f", digits = 3)
f2 <- function(x) formatC(x, format = "f", digits = 2)
sg <- function(x, d = 3) { x <- round(x, d); paste0(ifelse(x == 0, "", ifelse(x < 0, "-", "+")), formatC(abs(x), format = "f", digits = d)) }
pf <- function(p) if (is.na(p)) "n/a" else formatC(p, format = "f", digits = 3)
RHO <- 0.30 / 0.70                                   # n_test / n_train for 30% of subjects held out
paired <- function(x, y) {                           # x, y ordered by split
  d <- x - y; J <- length(d); m <- mean(d)
  se <- sqrt((1 / J + RHO) * var(d)); h <- qt(0.975, J - 1) * se
  list(m = m, wins = sum(d > 0), J = J, lo = m - h, hi = m + h,
       p = if (all(d == 0)) NA else suppressWarnings(wilcox.test(x, y, paired = TRUE, exact = TRUE)$p.value),
       pnb = 2 * pt(-abs(m / se), J - 1))
}
emit_pair <- function(name, pr, d = 3) {
  M(paste0(name, "D"), sprintf("\\ensuremath{%s}", sg(pr$m, d))); M(paste0(name, "P"), pf(pr$p)); M(paste0(name, "W"), pr$wins)
  M(paste0(name, "Abs"), formatC(abs(pr$m), format = "f", digits = d))
  M(paste0(name, "CI"), sprintf("\\ensuremath{[%s,\\,%s]}", sg(pr$lo, d), sg(pr$hi, d)))
  M(paste0(name, "Sig"), if (pr$lo > 0 || pr$hi < 0) "yes" else "no")
}
AT <- c(sch100 = "Sch", power229 = "Pow"); OC <- c(age_continuous = "A", nbk_0v1 = "T")
FR <- c("0.25" = "Q", "0.5" = "H", "0.75" = "R", "1" = "F")

## ---------------- E3 grid: corrected intervals for the comparisons already in the paper ----
kb <- read.csv(file.path(R, if (file.exists(file.path(R, "k_by_n_grouped.csv"))) "k_by_n_grouped.csv" else "k_by_n.csv"))
kget <- function(at, oc, design, kr, fr) { x <- kb[kb$atlas == at & kb$outcome == oc & kb$design == design &
  kb$train_frac == fr & (if (is.na(kr)) is.na(kb$K_or_R) else kb$K_or_R %in% kr), ]; x$metric[order(x$split)] }
for (at in names(AT)) for (oc in names(OC)) {
  tg <- paste0(AT[at], OC[oc])
  two <- kget(at, oc, "two_sided_ard", 10, 1); one <- kget(at, oc, "one_sided_ridge", 10, 1); edge <- kget(at, oc, "edge_ridge_cv", NA, 1)
  emit_pair(paste0("nbTwoEdge", tg), paired(two, edge)); emit_pair(paste0("nbOneEdge", tg), paired(one, edge))
  emit_pair(paste0("nbOneTwo", tg), paired(one, two))
  emit_pair(paste0("nbKTwentyTen", tg), paired(kget(at, oc, "two_sided_ard", 20, 1), two))
  M(paste0("kbEdgeAll", tg), f3(mean(edge)))
}

## ---------------- G1: generic reductions and matched learners -----------------------------
g <- do.call(rbind, lapply(list.files(R, "^generic_(sch100|power229)_.*\\.csv$", full.names = TRUE), read.csv))
g$arm <- paste0(g$design, ifelse(is.na(g$size), "", g$size), "_", g$learner)
ga <- aggregate(cbind(metric, rmse) ~ atlas + outcome + arm + train_frac + split, g, mean, na.action = na.pass)
gget <- function(at, oc, arm, fr, col = "metric") { x <- ga[ga$atlas == at & ga$outcome == oc & ga$arm == arm & ga$train_frac == fr, ]
  x[[col]][order(x$split)] }
ARM <- c(two_sided10_ard = "TwoArd", two_sided10_ridge = "TwoRidge", pca55_ard = "PcaArd", pca55_ridge = "PcaRidge",
         rp55_ard = "RpSArd", rp55_ridge = "RpS", rp210_ridge = "RpM", rp1000_ridge = "RpL",
         one_sided10_ridge = "One", edge_ridge = "Edge")
PAIRS <- list(c("pca55_ridge", "edge_ridge"), c("pca55_ridge", "two_sided10_ard"), c("pca55_ridge", "two_sided10_ridge"),
              c("two_sided10_ridge", "two_sided10_ard"), c("two_sided10_ridge", "edge_ridge"), c("two_sided10_ridge", "rp55_ridge"),
              c("two_sided10_ard", "rp55_ard"), c("rp1000_ridge", "edge_ridge"), c("one_sided10_ridge", "edge_ridge"),
              c("pca55_ridge", "one_sided10_ridge"), c("pca55_ard", "two_sided10_ard"), c("one_sided10_ridge", "two_sided10_ridge"))
for (at in names(AT)) for (oc in names(OC)) for (fr in c(0.25, 0.5, 0.75, 1)) {
  tg <- paste0(AT[at], OC[oc], FR[as.character(fr)])
  for (a in names(ARM)) { v <- gget(at, oc, a, fr); if (!length(v)) next
    M(paste0("g", tg, ARM[a]), f3(mean(v)))
    if (oc == "age_continuous") M(paste0("g", tg, ARM[a], "Rmse"), f2(mean(gget(at, oc, a, fr, "rmse")))) }
  if (fr %in% c(0.25, 1)) for (pr in PAIRS) { x <- gget(at, oc, pr[1], fr); y <- gget(at, oc, pr[2], fr)
    if (length(x) == 10 && length(y) == 10) emit_pair(paste0("gd", tg, ARM[pr[1]], "v", ARM[pr[2]]), paired(x, y)) }
}
## between-draw spread of the random projection (sd over the 5 draws, averaged over splits)
for (at in names(AT)) { x <- g[g$atlas == at & g$outcome == "age_continuous" & g$arm == "rp55_ridge" & g$train_frac == 1, ]
  M(paste0("gRpDrawSd", AT[at]), f3(mean(tapply(x$metric, x$split, sd)))) }
M("gNdraws", length(unique(g$draw[g$design == "rp"])))
## table: one row per arm, both tasks, both resolutions, full training size
lab <- c(two_sided10_ard = "T, ARD", two_sided10_ridge = "T, ridge", pca55_ard = "PCA, ARD", pca55_ridge = "PCA, ridge",
         rp55_ard = "RP, ARD", rp55_ridge = "RP, ridge", rp1000_ridge = "RP, ridge", one_sided10_ridge = "O, ridge", edge_ridge = "E, ridge")
nf <- c(two_sided10_ard = "55", two_sided10_ridge = "55", pca55_ard = "55", pca55_ridge = "55", rp55_ard = "55", rp55_ridge = "55",
        rp1000_ridge = "1{,}000", one_sided10_ridge = "$VR$", edge_ridge = "$p$")
tab <- character(0)
for (a in names(lab)) {
  cell <- function(at, oc, col = "metric", fmt = f3) { v <- gget(at, oc, a, 1, col); if (length(v)) fmt(mean(v)) else "--" }
  tab <- c(tab, sprintf("%s & %s & %s & %s & %s & %s & %s & %s \\\\", lab[a], nf[a],
    cell("sch100", "nbk_0v1"), cell("power229", "nbk_0v1"), cell("sch100", "age_continuous"), cell("power229", "age_continuous"),
    cell("sch100", "age_continuous", "rmse", f2), cell("power229", "age_continuous", "rmse", f2)))
  if (a %in% c("two_sided10_ridge", "rp1000_ridge", "one_sided10_ridge")) tab <- c(tab, "\\midrule")
}
tab[length(tab)] <- sub(" \\\\\\\\$", "", tab[length(tab)])          # the .tex file supplies the last row end
writeLines(tab, file.path(paper, "tab_generic.tex"))
## figure data: accuracy against training size
for (at in names(AT)) for (oc in names(OC)) for (a in names(ARM)) {
  x <- ga[ga$atlas == at & ga$outcome == oc & ga$arm == a, ]; if (!nrow(x)) next
  n <- aggregate(n_train ~ train_frac, g[g$atlas == at & g$outcome == oc, ], function(v) round(mean(v)))
  m <- merge(aggregate(metric ~ train_frac, x, mean), n)
  write.table(data.frame(ntrain = m$n_train, metric = round(m$metric, 4)), file.path(paper, "data", sprintf("g_%s_%s_%s.dat", AT[at], OC[oc], ARM[a])),
              row.names = FALSE, quote = FALSE)
}

## ---------------- G2: external collections --------------------------------------------------
if (file.exists(file.path(R, "external_collections.csv"))) {
  e <- read.csv(file.path(R, "external_collections.csv"))
  e$arm <- paste0(e$design, ifelse(is.na(e$size), "", e$size), "_", e$learner)
  ea <- aggregate(cbind(spearman, rmse) ~ collection + arm + split, e, mean)
  eget <- function(co, arm, col) { x <- ea[ea$collection == co & ea$arm == arm, ]; x[[col]][order(x$split)] }
  EA <- c(two_sided10_ard = "TwoArd", two_sided10_ridge = "TwoRidge", pca55_ridge = "PcaRidge", pca55_ard = "PcaArd", rp55_ridge = "RpS",
          rp1000_ridge = "RpL", one_sided10_ridge = "One", edge_ridge = "Edge", edge_ridge_fixed50 = "EdgeFixed", baseline_train_mean = "Base")
  CO <- c(abide = "Abide", acpi = "Acpi"); RHO <- 0.30 / 0.70
  rows <- character(0)
  for (co in names(CO)) {
    for (a in names(EA)) { M(paste0("x", CO[co], EA[a], "Rmse"), f2(mean(eget(co, a, "rmse")))); M(paste0("x", CO[co], EA[a], "Sp"), f3(mean(eget(co, a, "spearman")))) }
    for (a in setdiff(names(EA), c("edge_ridge", "baseline_train_mean")))
      emit_pair(paste0("xd", CO[co], EA[a], "vEdge"), paired(eget(co, a, "rmse"), eget(co, "edge_ridge", "rmse")), 2)
    emit_pair(paste0("xd", CO[co], "PcaRidgevTwoArd"), paired(eget(co, "pca55_ridge", "rmse"), eget(co, "two_sided10_ard", "rmse")), 2)
    emit_pair(paste0("xd", CO[co], "TwoArdvRpS"), paired(eget(co, "two_sided10_ard", "rmse"), eget(co, "rp55_ridge", "rmse")), 2)
    M(paste0("x", CO[co], "Ntrain"), e$n_train[e$collection == co][1])
    ms <- function(a) sprintf("%s", f2(mean(eget(co, a, "rmse"))))
    rows <- c(rows, sprintf("%s & %s & %s & %s & %s & %s & %s & %s \\\\", toupper(co), ms("two_sided10_ard"), ms("two_sided10_ridge"),
      ms("pca55_ridge"), ms("rp55_ridge"), ms("one_sided10_ridge"), ms("edge_ridge"), ms("baseline_train_mean")))
  }
  rows[length(rows)] <- sub(" \\\\\\\\$", "", rows[length(rows)])
  writeLines(rows, file.path(paper, "tab_external.tex"))
}

## ---------------- G3: cost against the number of training graphs ----------------------------
if (file.exists(file.path(R, "scaling_n.csv"))) {
  s <- read.csv(file.path(R, "scaling_n.csv")); s$key <- paste(s$design, s$stage)
  med <- aggregate(seconds ~ n_train + key, s, median)
  KEY <- c("two_sided fit_ard_logit" = "TwoArd", "two_sided fit_ridge_dual_cv" = "TwoRidge", "one_sided fit_ridge_dual_cv" = "One",
           "edge_store fit_ridge_dual_cv" = "Edge", "two_sided build" = "BuildTwo", "one_sided build" = "BuildOne",
           "template mean_plus_eigen" = "Template", "edge_store build" = "BuildEdge")
  nm <- c("200" = "A", "400" = "B", "800" = "C", "1600" = "D", "3200" = "E")
  for (k in names(KEY)) { x <- med[med$key == k, ]; x <- x[order(x$n_train), ]
    write.table(data.frame(n = x$n_train, sec = signif(x$seconds, 4)), file.path(paper, "data", sprintf("g_n_%s.dat", KEY[k])), row.names = FALSE, quote = FALSE)
    for (i in seq_len(nrow(x))) M(paste0("sn", KEY[k], nm[as.character(x$n_train[i])]),
                                 if (x$seconds[i] >= 100) formatC(x$seconds[i], format = "f", digits = 0, big.mark = "{,}") else if (x$seconds[i] >= 10) formatC(x$seconds[i], format = "f", digits = 1) else f2(x$seconds[i]))
    if (nrow(x) >= 3) M(paste0("sn", KEY[k], "Slope"), f2(coef(lm(log(seconds) ~ log(n_train), x[x$n_train >= 800, ]))[2]))
  }
  top <- max(s$n_train); e <- med[med$key == "edge_store fit_ridge_dual_cv", ]; t <- med[med$key == "two_sided fit_ard_logit", ]; o <- med[med$key == "one_sided fit_ridge_dual_cv", ]
  M("snTop", formatC(top, big.mark = "{,}", format = "d"))
  M("snEdgeTwoRatioTop", formatC(e$seconds[e$n_train == top] / t$seconds[t$n_train == top], format = "f", digits = 0, big.mark = "{,}"))
  M("snEdgeOneRatioTop", formatC(e$seconds[e$n_train == top] / o$seconds[o$n_train == top], format = "f", digits = 0))
  M("snEdgeTwoRatioA", formatC(e$seconds[e$n_train == 200] / t$seconds[t$n_train == 200], format = "f", digits = 0))
  M("snEdgeMinTop", f2(e$seconds[e$n_train == top] / 60))
  mm <- function(d) 8 * top * max(s$n_features[s$design == d & s$n_train == top], na.rm = TRUE)   # bytes, float64
  M("snEdgeGBTop", f2(mm("edge_store") / 1e9)); M("snOneMBTop", formatC(mm("one_sided") / 1e6, format = "f", digits = 1)); M("snTwoMBTop", f2(mm("two_sided") / 1e6))
}

## ---------------- G10: primal ridge on the synopses against the number of training graphs ---------
if (file.exists(file.path(R, "scaling_n_primal.csv"))) {
  sp <- read.csv(file.path(R, "scaling_n_primal.csv")); mp <- aggregate(seconds ~ n_train + design, sp, median)
  fs <- function(x) if (x >= 100) formatC(x, format = "f", digits = 0, big.mark = "{,}") else if (x >= 10) formatC(x, format = "f", digits = 1) else f2(x)
  for (dsg in c("two_sided", "one_sided")) { x <- mp[mp$design == dsg, ]; x <- x[order(x$n_train), ]; nm <- if (dsg == "two_sided") "TwoPrimal" else "OnePrimal"
    write.table(data.frame(n = x$n_train, sec = signif(x$seconds, 4)), file.path(paper, "data", sprintf("g_n_%s.dat", nm)), row.names = FALSE, quote = FALSE)
    M(paste0("sn", nm, "A"), fs(x$seconds[x$n_train == 200])); M(paste0("sn", nm, "E"), fs(x$seconds[x$n_train == max(x$n_train)]))
    M(paste0("sn", nm, "Slope"), f2(coef(lm(log(seconds) ~ log(n_train), x[x$n_train >= 800, ]))[2])) }
  s0 <- read.csv(file.path(R, "scaling_n.csv")); top <- max(s0$n_train)
  e <- median(s0$seconds[s0$design == "edge_store" & s0$stage == "fit_ridge_dual_cv" & s0$n_train == top])
  M("snEdgeTwoPrimalRatioTop", formatC(e / mp$seconds[mp$design == "two_sided" & mp$n_train == top], format = "f", digits = 0, big.mark = "{,}"))
  M("snEdgeOnePrimalRatioTop", formatC(e / mp$seconds[mp$design == "one_sided" & mp$n_train == top], format = "f", digits = 0))
}
## ---------------- G5: inner folds over graphs or over subjects (binary task) ----------------
gf <- list.files(R, "^grouped_cv_.*\\.csv$", full.names = TRUE)
if (length(gf)) {
  q <- do.call(rbind, lapply(gf, read.csv))
  GA <- c(edge = "Edge", one_sided_5 = "OneFive", one_sided_10 = "One", two_sided_10 = "TwoRidge", pca_55 = "PcaRidge")
  qget <- function(at, arm, fold, fr) { x <- q[q$atlas == at & q$arm == arm & q$folds == fold & q$train_frac == fr, ]; x$auc[order(x$split)] }
  for (at in names(AT)) for (a in names(GA)) {
    xs <- qget(at, a, "subject", 1); xg <- qget(at, a, "graph", 1); if (length(xs) != 10) next
    M(paste0("q", AT[at], GA[a], "Subj"), f3(mean(xs))); M(paste0("q", AT[at], GA[a], "Graph"), f3(mean(xg)))
    emit_pair(paste0("qd", AT[at], GA[a]), paired(xs, xg))
    two <- gget(at, "nbk_0v1", "two_sided10_ard", 1)
    if (length(two) == 10) emit_pair(paste0("qd", AT[at], "TwoArdv", GA[a], "Subj"), paired(two, xs))
  }
}
## ---------------- G6: effective df on the binary task with inner folds by subject ------------
df_files <- list.files(R, "^df_grouped_.*\\.csv$", full.names = TRUE)
if (length(df_files)) {
  z <- do.call(rbind, lapply(df_files, read.csv))
  for (at in names(AT)) {
    zz <- z[z$atlas == at, ]; if (length(unique(zz$split)) != 10) next
    v <- function(m, fo, col = "eff_df") { x <- zz[zz$method == m & zz$folds == fo, ]; x[[col]][order(x$split)] }
    two <- v("two_sided_ard", "none"); rs <- v("edge_ridge_logit", "subject"); rg <- v("edge_ridge_logit", "graph"); ls <- v("edge_lasso_logit", "subject")
    tg <- AT[at]; f1 <- function(x) formatC(x, format = "f", digits = 1)
    M(paste0("dg", tg, "Two"), f1(mean(two))); M(paste0("dg", tg, "TwoMin"), f1(min(two))); M(paste0("dg", tg, "TwoMax"), f1(max(two)))
    M(paste0("dg", tg, "Ridge"), f1(mean(rs))); M(paste0("dg", tg, "RidgeMin"), f1(min(rs))); M(paste0("dg", tg, "RidgeMax"), f1(max(rs)))
    M(paste0("dg", tg, "RidgeGraph"), f1(mean(rg))); M(paste0("dg", tg, "Ratio"), f1(mean(rs) / mean(two)))
    M(paste0("dg", tg, "Lasso"), f1(mean(ls))); M(paste0("dg", tg, "LassoMin"), min(ls)); M(paste0("dg", tg, "LassoMax"), max(ls))
    M(paste0("dg", tg, "RidgeAuc"), f3(mean(v("edge_ridge_logit", "subject", "auc")))); M(paste0("dg", tg, "RidgeGraphAuc"), f3(mean(v("edge_ridge_logit", "graph", "auc"))))
    M(paste0("dg", tg, "LassoAuc"), f3(mean(v("edge_lasso_logit", "subject", "auc")))); M(paste0("dg", tg, "TwoAuc"), f3(mean(v("two_sided_ard", "none", "auc"))))
    emit_pair(paste0("dgd", tg, "TwovRidgeAuc"), paired(v("two_sided_ard", "none", "auc"), v("edge_ridge_logit", "subject", "auc")))
    M(paste0("dg", tg, "Ntrain"), 201)
    p <- if (at == "sch100") 4950 else 26106
    write.table(data.frame(nominal = p, df = round(rs, 3)), file.path(paper, "data", sprintf("g_df_%s_ridge_T.dat", tg)), row.names = FALSE, quote = FALSE)
    write.table(data.frame(nominal = p, df = pmax(ls, 1)), file.path(paper, "data", sprintf("g_df_%s_lasso_T.dat", tg)), row.names = FALSE, quote = FALSE)
    write.table(data.frame(nominal = 55, df = round(two, 3)), file.path(paper, "data", sprintf("g_df_%s_two_T.dat", tg)), row.names = FALSE, quote = FALSE)
  }
}
## ---------------- range of the effective-df ratio (edge ridge over two-sided ARD), four settings ----
dfr <- function(f, ridge, two) { x <- read.csv(file.path(R, f)); x <- x[x$outcome == "age_continuous", ]
  mean(x$eff_df[x$method == ridge]) / mean(x$eff_df[x$method == two]) }
zr <- do.call(rbind, lapply(list.files(R, "^df_grouped_.*\\.csv$", full.names = TRUE), read.csv))
rat <- c(dfr("df_sch100.csv", "raw_edge_ridge", "gs_ard"), dfr("df_power229.csv", "raw_edge_ridge", "gs_ard"),
         sapply(split(zr, zr$atlas), function(x) mean(x$eff_df[x$method == "edge_ridge_logit" & x$folds == "subject"]) / mean(x$eff_df[x$method == "two_sided_ard"])))
M("dfRatioLow", round(min(rat))); M("dfRatioHigh", round(max(rat)))
gapshare <- function(at) 1 - mean(gget(at, "age_continuous", "edge_ridge", 1) - gget(at, "age_continuous", "two_sided10_ridge", 1)) /
                              mean(gget(at, "age_continuous", "edge_ridge", 1) - gget(at, "age_continuous", "two_sided10_ard", 1))
M("gLearnerShareSch", round(100 * gapshare("sch100"))); M("gLearnerSharePow", round(100 * gapshare("power229")))
## ---------------- edge ARD cells of the vertex sweep (E1): measured, capped, memory-bound -----------
scv <- read.csv(file.path(R, "scaling.csv")); ar <- scv[scv$learner == "ard_iter", ]
okv <- sort(unique(ar$V[ar$status == "ok"])); vmid <- max(okv)
spi <- function(v) median(ar$seconds_per_iter[ar$V == v & ar$status == "ok"])
M("scArdMidV", vmid); M("scArdIterMid", formatC(spi(vmid), format = "f", digits = 0))
M("scArdMidPred", formatC(spi(100) * ((vmid * (vmid - 1) / 2) / 4950)^3, format = "f", digits = 0))
M("scArdCapV", min(ar$V[ar$status == "timeout"])); M("scArdMemGB", formatC(6 * 26106^2 * 8 / 2^30, format = "f", digits = 1))
M("scArdProjPowDays", formatC(spi(100) * (26106 / 4950)^3 * 500 / 86400, format = "f", digits = 0))
## ---------------- arithmetic and the earlier one-sided ARD and kernel ridge runs at V = 229 ---
M("gPcaLoadMB", formatC(26106 * 55 * 8 / 1e6, format = "f", digits = 1)); M("gBasisKB", formatC(229 * 10 * 8 / 1e3, format = "f", digits = 1))
ks <- read.csv(file.path(R, "ksweep_power229.csv"))
kv <- function(oc, variant, col) { x <- ks[ks$outcome == oc & ks$variant == variant & ks$status == "ok", ]; x[[col]][order(x$split)] }
emit_pair("nbOneArdTwoPowA", paired(kv("age_continuous", "onesided", "spearman"), kv("age_continuous", "current", "spearman")))
emit_pair("nbOneArdTwoPowT", paired(kv("nbk_0v1", "onesided", "auc"), kv("nbk_0v1", "current", "auc")))
M("krPowT", f3(mean(kv("nbk_0v1", "kundu", "auc")))); M("krPowA", f3(mean(kv("age_continuous", "kundu", "spearman"))))
M("krPowARmse", f2(mean(kv("age_continuous", "kundu", "rmse")))); M("oneArdPowARmse", f2(mean(kv("age_continuous", "onesided", "rmse"))))
M("twoArdPowARmse", f2(mean(kv("age_continuous", "current", "rmse"))))
envf <- file.path(R, "scaling_n_env.txt")
if (file.exists(envf)) { ev <- readLines(envf); gv <- function(k) sub(paste0("^", k, ":\\s*"), "", grep(paste0("^", k, ":"), ev, value = TRUE)[1])
  M("snPlatformText", sprintf("%s desktop (macOS, 12 cores), %s GB RAM, %s, R reference BLAS (single-threaded)", gv("cpu"), round(as.numeric(gv("ram_gb"))),
                             sub(" \\(.*", "", sub("^R version", "R", ev[1])))) }
## ---------------- G7: queries answered from a cached posterior (pilot, Age) ---------------------
if (file.exists(file.path(R, "projpost_accuracy.csv"))) {
  pa <- read.csv(file.path(R, "projpost_accuracy.csv")); pq <- read.csv(file.path(R, "projpost_queries.csv"))
  pv <- function(at, me, k) { x <- pa[pa$atlas == at & pa$method == me & pa$k == k, ]; x$spearman[order(x$split)] }
  KN <- c("5" = "Five", "10" = "Ten", "20" = "Twenty"); rowsP <- character(0)
  for (at in names(AT)) { tg <- AT[at]; x <- pq[pq$atlas == at, ]
    M(paste0("ppDense", tg), f3(mean(pv(at, "dense_posterior_mean", 55))))
    M(paste0("ppCacheKB", tg), round(mean(x$cache_bytes) / 1e3)); M(paste0("ppEdgeMB", tg), formatC(mean(x$edge_store_bytes) / 1e6, format = "f", digits = 1))
    M(paste0("ppCover", tg), round(100 * mean(x$pi90_coverage))); M(paste0("ppWidth", tg), round(mean(x$pi90_mean_width)))
    M(paste0("ppNodeSec", tg), f2(mean(x$node_query_seconds)))
    tb <- table(unlist(lapply(strsplit(x$pip10_idx, " "), as.integer))); M(paste0("ppStablePairs", tg), sum(tb >= 8))
    tn <- table(unlist(lapply(strsplit(x$top10_nodes, " "), as.integer))); M(paste0("ppStableNodes", tg), sum(tn >= 8))
    for (k in names(KN)) { M(paste0("ppProj", KN[k], tg), f3(mean(pv(at, "projection_of_mean", as.numeric(k)))))
      M(paste0("ppRefit", KN[k], tg), f3(mean(pv(at, "lasso_refit_on_cores", as.numeric(k)))))
      emit_pair(paste0("ppd", KN[k], tg), paired(pv(at, "projection_of_mean", as.numeric(k)), pv(at, "lasso_refit_on_cores", as.numeric(k)))) }
  }
  for (k in names(KN)) rowsP <- c(rowsP, sprintf("%s & %s & %s & %s & %s \\\\", k, f3(mean(pv("sch100", "projection_of_mean", as.numeric(k)))),
    f3(mean(pv("sch100", "lasso_refit_on_cores", as.numeric(k)))), f3(mean(pv("power229", "projection_of_mean", as.numeric(k)))), f3(mean(pv("power229", "lasso_refit_on_cores", as.numeric(k))))))
  rowsP <- c(rowsP, sprintf("55 (dense) & %s & & %s &", f3(mean(pv("sch100", "dense_posterior_mean", 55))), f3(mean(pv("power229", "dense_posterior_mean", 55)))))
  M("ppDraws", 500)
}
## ---------------- G8: cached-posterior queries on four collections (table and macros) ------------
if (file.exists(file.path(R, "projpost_ext_accuracy.csv"))) {
  xa <- read.csv(file.path(R, "projpost_ext_accuracy.csv")); xq <- read.csv(file.path(R, "projpost_ext_queries.csv"))
  CL <- c(sch100 = "Main, $V{=}100$", power229 = "Main, $V{=}229$", abide = "ABIDE", acpi = "ACPI")
  xv <- function(co, me, k) { x <- xa[xa$collection == co & xa$method == me & xa$k == k, ]; x$spearman[order(x$split)] }
  xm <- function(co, me, g) { x <- xa[xa$collection == co & xa$method == me & !is.na(xa$g) & xa$g == g, ]; x$spearman[order(x$split)] }
  rowsX <- character(0); nall <- 0; nzero <- 0
  for (co in names(CL)) { q <- xq[xq$collection == co, ]
    rowsX <- c(rowsX, sprintf("%s & %s & %s & %s & %s & %s & %d & %d \\\\", CL[co], f3(mean(q$dense_spearman)),
      f3(mean(xv(co, "projection_l1", 20))), f3(mean(xv(co, "lasso_refit", 20))), f3(mean(xm(co, "projection_modes", 5))),
      f3(mean(xm(co, "ridge_refit_same_modes", 5))), round(100 * mean(q$cache_cover90)), round(mean(q$cache_bytes) / 1e3)))
    for (k in c(5, 10, 20)) { pr <- paired(xv(co, "projection_l1", k), xv(co, "lasso_refit", k)); nall <- nall + 1; nzero <- nzero + (pr$lo <= 0 && pr$hi >= 0) }
    for (gg in c(3, 5, 7)) { pr <- paired(xm(co, "projection_modes", gg), xm(co, "ridge_refit_same_modes", gg)); nall <- nall + 1; nzero <- nzero + (pr$lo <= 0 && pr$hi >= 0) }
  }
  rowsX[length(rowsX)] <- sub(" \\\\\\\\$", "", rowsX[length(rowsX)])
  writeLines(rowsX, file.path(paper, "tab_projpost.tex"))
  M("pxPairsN", nall); M("pxPairsTie", nzero)
  M("pxCoverMin", round(100 * min(tapply(xq$cache_cover90, xq$collection, mean)))); M("pxCoverMax", round(100 * max(tapply(xq$cache_cover90, xq$collection, mean))))
  M("pxCacheMinKB", round(min(xq$cache_bytes) / 1e3)); M("pxCacheMaxKB", round(max(xq$cache_bytes) / 1e3))
  M("pxEdgeMinMB", formatC(min(tapply(xq$edge_store_bytes, xq$collection, mean)) / 1e6, format = "f", digits = 1))
  M("pxEdgeMaxMB", formatC(max(tapply(xq$edge_store_bytes, xq$collection, mean)) / 1e6, format = "f", digits = 1))
  bc <- tapply(xq$boot_cover90, xq$collection, mean); M("pxBootCoverMin", round(100 * min(bc))); M("pxBootCoverMax", round(100 * max(bc)))
  bs <- tapply(xq$boot_seconds, xq$collection, mean); M("pxBootMinS", formatC(min(bs), format = "f", digits = 1)); M("pxBootMaxS", formatC(max(bs), format = "f", digits = 1)); M("pxBootB", max(xq$boot_B))
}
## ---------------- G9: two-sided size sweep with ridge ----------------------------------------
if (file.exists(file.path(R, "ksweep_ridge.csv"))) {
  kr <- read.csv(file.path(R, "ksweep_ridge.csv"))
  kg <- function(at, oc, K, fr) { x <- kr[kr$atlas == at & kr$outcome == oc & kr$K == K & kr$train_frac == fr, ]; x$metric[order(x$split)] }
  nbig <- 0; nsig <- 0; nsmall <- 0; ncell <- 0
  for (at in names(AT)) for (oc in names(OC)) {
    for (K in c(15, 20, 30)) { pr <- paired(kg(at, oc, K, 1), kg(at, oc, 10, 1)); nbig <- nbig + 1; nsig <- nsig + (pr$lo > 0 || pr$hi < 0)
      if (K == 20) { emit_pair(paste0("krTwentyTen", AT[at], OC[oc]), pr); M(paste0("krTwenty", AT[at], OC[oc]), f3(mean(kg(at, oc, 20, 1)))) } }
    M(paste0("krThirty", AT[at], OC[oc]), f3(mean(kg(at, oc, 30, 1)))); M(paste0("krTen", AT[at], OC[oc]), f3(mean(kg(at, oc, 10, 1))))
    for (fr in c(0.25, 0.5, 0.75, 1)) { mk <- sapply(c(5, 10, 15, 20, 30), function(K) mean(kg(at, oc, K, fr))); ncell <- ncell + 1
      nsmall <- nsmall + (c(5, 10, 15, 20, 30)[which.max(mk)] <= 10) }
  }
  M("krBigN", nbig); M("krBigSig", nsig); M("krBestSmall", nsmall); M("krBestN", ncell)
}
writeLines(c("% Generated by Codes/icde/make_g_numbers.R. Do not edit by hand.", unlist(mac)), file.path(paper, "numbers_g.tex"))
cat("wrote", length(mac), "macros to", file.path(paper, "numbers_g.tex"), "\n")
