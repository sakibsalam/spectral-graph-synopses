#!/usr/bin/env python3
"""Builds numbers.tex, data/*.dat and tab_*.tex for the ICDE paper from the repo result CSVs.
Usage: python3 make_numbers.py [REPO_DIR] [ICDE_RESULTS_DIR] [OUTPUT_DIR]
OUTPUT_DIR defaults to the folder of this script.
Stdlib only. Pending experiment outputs (results/icde/*.csv) give \\pending{name} macros when absent."""
import csv, os, sys, statistics as st, itertools
REPO = sys.argv[1] if len(sys.argv) > 1 else "."
HERE = sys.argv[3] if len(sys.argv) > 3 else os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(HERE, "data"), exist_ok=True)
RES = os.path.join(REPO, "results")
IRES = sys.argv[2] if len(sys.argv) > 2 else os.path.join(RES, "icde")
M = {}          # macro name -> text

def rd(path):
    if not os.path.exists(path): return None
    with open(path, newline="") as f: return list(csv.DictReader(f))
def fl(x):
    try: return float(x)
    except (TypeError, ValueError): return None
def commas(n): return f"{int(round(n)):,}".replace(",", "{,}")
def f3(x): return f"{x:.3f}"
def f2(x): return f"{x:.2f}"
def f1(x): return f"{x:.1f}"
def sec(x):
    return f"{x:.2f}" if x < 1 else (f"{x:.1f}" if x < 100 else f"{x:.0f}")
def pv(p): return f"{p:.3f}"
def put(name, val): M[name] = str(val)
def pending(name): M[name] = r"\pending{%s}" % name

def wilcoxon_exact(d):
    d = [x for x in d if abs(x) > 1e-12]
    n = len(d)
    if n == 0: return 1.0
    ab = sorted(abs(x) for x in d)
    # average ranks, doubled to integers
    rank = {}
    i = 0
    while i < n:
        j = i
        while j + 1 < n and abs(ab[j + 1] - ab[i]) < 1e-12: j += 1
        for k in range(i, j + 1): rank[k] = (i + j + 2)  # 2 * mean rank
        i = j + 1
    rk = []
    used = {}
    for x in d:
        idx = min(k for k in range(n) if abs(ab[k] - abs(x)) < 1e-12)
        rk.append(rank[idx])
    tot = sum(rk); obs = sum(r for r, x in zip(rk, d) if x > 0)
    counts = {0: 1}
    for r in rk:
        new = dict(counts)
        for s, c in counts.items(): new[s + r] = new.get(s + r, 0) + c
        counts = new
    allc = sum(counts.values()); dev = abs(obs - tot / 2)
    return sum(c for s, c in counts.items() if abs(s - tot / 2) >= dev - 1e-9) / allc

def paired(a, b):
    ks = sorted(set(a) & set(b)); d = [a[k] - b[k] for k in ks]
    return st.mean(d), sum(1 for x in d if x > 0), len(d), wilcoxon_exact(d)

def col(rows, flt, key, by="split"):
    return {int(r[by]): fl(r[key]) for r in rows if flt(r) and fl(r[key]) is not None}
def mean(dct): return st.mean(dct.values())
def sd(dct): return st.stdev(dct.values())

# ------------------------------------------------------------------ dataset facts
subj = rd(os.path.join(REPO, "data/osf_master/subjects.csv"))
if subj:
    put("nGraphs", commas(len(subj))); put("nSubjects", len(set(r["SUBJECT"] for r in subj)))
    ages = [fl(r["AGE"]) for r in subj if fl(r["AGE"]) is not None]
    put("ageMin", int(min(ages))); put("ageMax", int(max(ages)))
    put("nCond", len(set(r["COND"] for r in subj)))
pS, pP = 100 * 99 // 2, 229 * 228 // 2
q = 10 * 11 // 2
put("pSch", commas(pS)); put("pPow", commas(pP)); put("qTen", q)
put("ratioSch", f"{pS/q:.0f}"); put("ratioPow", f"{pP/q:.0f}")
put("bytesEdgeSch", commas(8 * pS)); put("bytesEdgePow", commas(8 * pP))
put("bytesOneSch", commas(8 * 100 * 10)); put("bytesOnePow", commas(8 * 229 * 10))
put("bytesTwo", commas(8 * q))
N = int(M["nGraphs"].replace("{,}", "")) if subj else 1003
put("collEdgePowMB", f1(8 * pP * N / 1e6)); put("collTwoKB", f1((8 * q * N + 8 * 229 * 10) / 1e3))
put("collOnePowMB", f1((8 * 229 * 10 * N) / 1e6))
put("minWilcoxP", "0.002")

# ------------------------------------------------------------------ capacity at V=100
df = rd(os.path.join(RES, "model_df_comparison.csv"))
if df:
    g = lambda m, k: col(df, lambda r: r["model"] == m, k)
    for tag, m in [("Ours", "gs_ard_logit"), ("Ridge", "raw_edge_ridge"), ("Lasso", "raw_edge_lasso"),
                   ("Ssvs", "raw_edge_ssvs"), ("Kr", "kundu_kernel_ridge")]:
        e = g(m, "effective_df")
        put(f"dfSch{tag}", f2(mean(e))); put(f"dfSch{tag}Sd", f2(sd(e)))
        put(f"dfSch{tag}Min", f2(min(e.values()))); put(f"dfSch{tag}Max", f2(max(e.values())))
        put(f"secDf{tag}", sec(mean(g(m, "seconds"))))
    d, w, n, p = paired(g("gs_ard_logit", "effective_df"), g("raw_edge_ridge", "effective_df"))
    put("dfSchDiff", f"{d:+.2f}"); put("dfSchP", pv(p)); put("dfSchWins", w)
    put("dfSchLassoMinI", int(min(g("raw_edge_lasso", "effective_df").values())))
    put("dfSchLassoMaxI", int(max(g("raw_edge_lasso", "effective_df").values())))
    put("dfSchKrNom", commas(800)); put("dfSchRatioKr", f"{mean(g('kundu_kernel_ridge','effective_df'))/mean(g('gs_ard_logit','effective_df')):.0f}")
    for tag, m in [("Ours", "gs_ard_logit"), ("Ridge", "raw_edge_ridge"), ("Kr", "kundu_kernel_ridge")]:
        put(f"aucDf{tag}", f3(mean(g(m, "fitted_auc"))))
    # data files for Fig 3(a)
    for tag, m in [("ours", "gs_ard_logit"), ("ridge", "raw_edge_ridge"), ("lasso", "raw_edge_lasso"),
                   ("ssvs", "raw_edge_ssvs"), ("kr", "kundu_kernel_ridge")]:
        with open(os.path.join(HERE, "data", f"df_sch_{tag}.dat"), "w") as f:
            f.write("nominal df\n")
            for r in df:
                if r["model"] == m: f.write(f"{r['nominal_params']} {fl(r['effective_df']):.4f}\n")
        e = g(m, "effective_df")
        with open(os.path.join(HERE, "data", f"df_sch_{tag}_mean.dat"), "w") as f:
            nom = [r["nominal_params"] for r in df if r["model"] == m][0]
            f.write(f"nominal df\n{nom} {mean(e):.4f}\n")

# ------------------------------------------------------------------ accuracy at V=100
fm = rd(os.path.join(RES, "five_method_sch100.csv"))
if fm:
    a = lambda m, k="auc": col(fm, lambda r: r["method"] == m, k)
    ours, ridge, lasso, ssvs = a("gs_ard"), a("raw_edge_ridge"), a("raw_edge_lasso"), a("raw_edge_ssvs")
    put("aucSchOurs", f3(mean(ours))); put("aucSchRidge", f3(mean(ridge)))
    put("aucSchLasso", f3(mean(lasso))); put("aucSchSsvs", f3(mean(ssvs)))
    put("ridgeConstSch", sum(1 for v in ridge.values() if abs(v - 0.5) < 1e-9))
    for nm, x, y in [("OursRidge", ours, ridge), ("OursLasso", ours, lasso), ("OursSsvs", ours, ssvs), ("SsvsRidge", ssvs, ridge)]:
        d, w, n, p = paired(x, y)
        put(f"dSchT{nm}", f"{d:+.3f}"); put(f"wSchT{nm}", w); put(f"pSchT{nm}", pv(p))
    put("nTrainT", int(round(mean(a("gs_ard", "n_train")))))
    put("ssvsNnz", int(max(a("raw_edge_ssvs", "n_nonzero").values())))
nb = rd(os.path.join(RES, "nbk_0v1_sch100.csv"))
if nb:
    put("aucSchRidgeWide", f3(st.mean(fl(r["ridge_wide_auc"]) for r in nb)))
if df:
    put("aucSchRidgeDev", M["aucDfRidge"])
ag = rd(os.path.join(RES, "age_continuous_sch100.csv"))
if ag:
    b = lambda m, k="spearman": col(ag, lambda r: r["method"] == m and r["metric"] == k, "value")
    put("spSchOurs", f3(mean(b("gs_ard")))); put("spSchRidge", f3(mean(b("ridge")))); put("spSchLasso", f3(mean(b("lasso"))))
    put("rmseSchOurs", f2(mean(b("gs_ard", "rmse")))); put("rmseSchRidge", f2(mean(b("ridge", "rmse"))))
    put("rmseSchBase", f2(mean(b("baseline", "rmse"))))
    d, w, n, p = paired(b("gs_ard"), b("ridge")); put("dSchAOursRidge", f"{d:+.3f}"); put("wSchAOursRidge", w); put("pSchAOursRidge", pv(p))
    d, w, n, p = paired(b("gs_ard"), b("lasso")); put("dSchAOursLasso", f"{d:+.3f}"); put("pSchAOursLasso", pv(p))
    put("nSubjAge", len(set(r["n_subj"] for r in ag)) and list(set(r["n_subj"] for r in ag))[0])

# ------------------------------------------------------------------ accuracy at V=229
tc = rd(os.path.join(RES, "tuned_comparison_v1003.csv"))
if tc:
    t = lambda oc, arm, k: col(tc, lambda r: r["outcome"] == oc and r["arm"] == arm and r["permuted"] in ("FALSE", "False", "F"), k)
    for oc, k, tag in [("nbk_0v1", "auc", "T"), ("age_continuous", "spearman", "A")]:
        o = t(oc, "gs_ard_fixed", k)
        for nm, arm in [("Ridge", "ridge_old"), ("RidgeDev", "ridge_cv"), ("Lasso", "lasso_old"), ("LassoDev", "lasso_cv")]:
            put(f"{k}Pow{tag}{nm}", f3(mean(t(oc, arm, k))))
        put(f"{k}Pow{tag}Ours", f3(mean(o)))
        for nm, arm in [("Ridge", "ridge_old"), ("RidgeDev", "ridge_cv")]:
            d, w, n, p = paired(o, t(oc, arm, k)); put(f"dPow{tag}Ours{nm}", f"{d:+.3f}"); put(f"wPow{tag}Ours{nm}", w); put(f"pPow{tag}Ours{nm}", pv(p))
        put(f"nConstPow{tag}", sum(1 for v in t(oc, "ridge_cv", k).values() if abs(v - 0.5) < 1e-9) if k == "auc" else 0)
        sdd = lambda arm: mean(col(tc, lambda r: r["outcome"] == oc and r["arm"] == arm, "seconds"))
        put(f"secPow{tag}Ours", sec(sdd("gs_ard_fixed"))); put(f"secPow{tag}RidgeCv", sec(sdd("ridge_cv")))
        put(f"secPow{tag}RidgeOld", sec(sdd("ridge_old"))); put(f"secPow{tag}LassoOld", sec(sdd("lasso_old")))
    put("rmsePowAOurs", f2(mean(t("age_continuous", "gs_ard_fixed", "rmse")))); put("rmsePowARidge", f2(mean(t("age_continuous", "ridge_cv", "rmse"))))
    # inner CV selection (V=229)
    for oc, k, tag in [("nbk_0v1", "auc", "T"), ("age_continuous", "spearman", "A")]:
        fx, cv = t(oc, "gs_ard_fixed", k), t(oc, "gs_ard_cv", k)
        d, w, n, p = paired(cv, fx); put(f"dCvFixed{tag}", f"{d:+.3f}"); put(f"pCvFixed{tag}", pv(p))
        put(f"cvMean{tag}", f3(mean(cv))); put(f"fixMean{tag}", f3(mean(fx)))
    ks = [r for r in tc if r["arm"] == "gs_ard_cv" and r["permuted"] in ("FALSE", "False", "F")]
    for oc, tag in [("nbk_0v1", "T"), ("age_continuous", "A"), ("nbk_1v2", "C"), ("age_binary_median", "B")]:
        kk = [int(float(r["K"])) for r in ks if r["outcome"] == oc and r["K"] != "NA"]
        if kk: put(f"cvKmode{tag}", max(set(kk), key=kk.count)); put(f"cvKmin{tag}", min(kk)); put(f"cvKmax{tag}", max(kk))
    # max shift over four tasks, fixed vs cv, as accuracy change
    mx = 0
    for oc, k in [("nbk_0v1", "auc"), ("age_continuous", "spearman"), ("nbk_1v2", "auc"), ("age_binary_median", "auc")]:
        mx = max(mx, abs(mean(t(oc, "gs_ard_cv", k)) - mean(t(oc, "gs_ard_fixed", k))))
    put("cvMaxShift", f3(mx))

# ------------------------------------------------------------------ synopsis size and vertex identity (V=100)
dv = rd(os.path.join(RES, "dense_variants_sch100_combined.csv"))
if dv:
    def v(oc, var, k): return col(dv, lambda r: r["outcome"] == oc and r["variant"] == var and r["status"] == "ok", k)
    for oc, k, tag in [("nbk_0v1", "auc", "T"), ("age_continuous", "spearman", "A")]:
        for var, nm in [("current", "KTen"), ("bigK20", "KTwenty"), ("bigK30", "KThirty"), ("bigK40", "KForty"),
                        ("multibasis", "Multi"), ("onesided", "One"), ("onesided_rownorm", "Row"), ("kundu", "Kr"), ("ridge", "Ridge")]:
            put(f"{k}Var{tag}{nm}", f3(mean(v(oc, var, k))))
        base = v(oc, "current", k)
        for var, nm in [("bigK20", "KTwenty"), ("bigK30", "KThirty"), ("bigK40", "KForty"), ("multibasis", "Multi"), ("onesided", "One"), ("onesided_rownorm", "Row"), ("kundu", "Kr"), ("ridge", "Ridge")]:
            d, w, n, p = paired(v(oc, var, k), base)
            put(f"dVar{tag}{nm}", f"{d:+.3f}"); put(f"wVar{tag}{nm}", w); put(f"pVar{tag}{nm}", pv(p))
        d, w, n, p = paired(v(oc, "onesided", k), v(oc, "bigK40", k)); put(f"dVar{tag}OneVsForty", f"{d:+.3f}"); put(f"wVar{tag}OneVsForty", w); put(f"pVar{tag}OneVsForty", pv(p))
        d, w, n, p = paired(v(oc, "onesided", k), v(oc, "kundu", k)); put(f"dVar{tag}OneVsKr", f"{d:+.3f}"); put(f"pVar{tag}OneVsKr", pv(p))
        d, w, n, p = paired(base, v(oc, "kundu", k)); put(f"dVar{tag}TenVsKr", f"{d:+.3f}"); put(f"pVar{tag}TenVsKr", pv(p)); put(f"wVar{tag}TenVsKr", w)
        d, w, n, p = paired(base, v(oc, "ridge", k)); put(f"dVar{tag}TenVsRidge", f"{d:+.3f}"); put(f"pVar{tag}TenVsRidge", pv(p))
        # dat files for Fig 4
        with open(os.path.join(HERE, "data", f"ksweep_{tag}.dat"), "w") as f:
            f.write("K metric\n")
            for var, K in [("current", 10), ("bigK20", 20), ("bigK30", 30), ("bigK40", 40)]:
                for s, val in sorted(v(oc, var, k).items()): f.write(f"{K} {val:.4f}\n")
        with open(os.path.join(HERE, "data", f"ksweep_{tag}_mean.dat"), "w") as f:
            f.write("K metric\n")
            for var, K in [("current", 10), ("bigK20", 20), ("bigK30", 30), ("bigK40", 40)]:
                f.write(f"{K} {mean(v(oc, var, k)):.4f}\n")
    sc2 = lambda oc, var: mean(v(oc, var, "seconds"))
    for var, nm in [("current", "KTen"), ("bigK20", "KTwenty"), ("bigK30", "KThirty"), ("bigK40", "KForty"), ("onesided", "One"), ("kundu", "Kr"), ("ridge", "Ridge")]:
        put(f"secVarT{nm}", sec(sc2("nbk_0v1", var))); put(f"secVarA{nm}", sec(sc2("age_continuous", var)))
    for var, nm in [("bigK20", "KTwenty"), ("bigK30", "KThirty"), ("bigK40", "KForty")]:
        put(f"ratioT{nm}", f"{sc2('nbk_0v1', var)/sc2('nbk_0v1', 'current'):.0f}")
        bt, st_ = v("nbk_0v1", var, "auc"), v("nbk_0v1", "current", "auc")
        put(f"lowT{nm}", sum(1 for k_ in set(bt) & set(st_) if bt[k_] < st_[k_]))
    put("ratioOneT", f"{sc2('nbk_0v1','onesided')/sc2('nbk_0v1','current'):.0f}")
    put("ratioOneA", f"{sc2('age_continuous','onesided')/sc2('age_continuous','current'):.0f}")
    put("ratioKrT", f"{sc2('nbk_0v1','kundu')/sc2('nbk_0v1','current'):.2f}")
    put("dimOne", commas(1000)); put("dimForty", commas(820)); put("dimKr", commas(800))

# ------------------------------------------------------------------ derived constants (exact arithmetic, no timings)
cub = (pP / pS) ** 3
put("cubicFactor", f"{cub:.0f}")
put("nSplits", 10)
put("nTrainA", int(round(0.7 * int(M["nSubjects"])))) if "nSubjects" in M else pending("nTrainA")
bigV, bigN = 400, 40000
bigP = bigV * (bigV - 1) // 2
put("bigV", bigV); put("bigN", commas(bigN)); put("bigP", commas(bigP))
put("bigEdgeGB", f1(8 * bigP * bigN / 1e9)); put("bigTwoMB", f1(8 * q * bigN / 1e6))

# ------------------------------------------------------------------ compressibility table
cp = rd(os.path.join(RES, "spectral_compressibility.csv"))
if cp:
    names = []
    for r in cp:
        if r["dataset"] not in names: names.append(r["dataset"])
    lines = []
    for nm in names:
        rs = [r for r in cp if r["dataset"] == nm]
        cells = []
        for th in ("0.5", "0.8", "0.9"):
            rr = [r for r in rs if abs(float(r["threshold"]) - float(th)) < 1e-9][0]
            cells.append("n.r." if rr["K"] == "NA" else f"{int(float(rr['K']))}")
        lab = nm.replace("&", r"\&").replace("(func, V=200)", "(functional)").replace("(struct, ", "(structural, ")
        lines.append(f"{lab} & {rs[0]['V']} & {rs[0]['N']} & " + " & ".join(cells))
    with open(os.path.join(HERE, "tab_compress.tex"), "w") as f: f.write(" \\\\\n".join(lines) + "\n")

# ------------------------------------------------------------------ E1: scaling
sc = rd(os.path.join(IRES, "scaling.csv"))
def med(design, learner, V, key="seconds"):
    if not sc: return None
    xs = [fl(r[key]) for r in sc if r["design"] == design and r["learner"] == learner and int(r["V"]) == V
          and r["status"] == "ok" and fl(r[key]) is not None]
    return st.median(xs) if xs else None
E1 = [("scTwoA", "two_sided", "ard_logit_K10", 50), ("scTwoC", "two_sided", "ard_logit_K10", 100), ("scTwoB", "two_sided", "ard_logit_K10", 229),
      ("scTwoD", "two_sided", "ard_logit_K10", 400), ("scOneB", "one_sided", "ridge_dual_cv_R10", 229),
      ("scEdgeRidgeB", "edge_store", "ridge_dual_cv", 229), ("scEdgeRidgeD", "edge_store", "ridge_dual_cv", 400),
      ("scLassoB", "edge_store", "lasso_glmnet", 229), ("scBuildTwoB", "two_sided", "build", 229), 
      ("scTemplateB", "template", "mean_plus_eigen", 229), ("scTemplateD", "template", "mean_plus_eigen", 400)]
for nm, d, l, V in E1:
    x = med(d, l, V); put(nm, sec(x)) if x is not None else pending(nm)
for nm, d, l, V, key in [("scArdIterC", "edge_store", "ard_iter", 100, "seconds_per_iter"),
                          ("scArdIterB", "edge_store", "ard_iter", 229, "seconds_per_iter"), ("scBmgcIterB", "edge_store", "bmgc_iter", 229, "seconds_per_iter"),
                          ("scBmgcIterC", "edge_store", "bmgc_iter", 100, "seconds_per_iter")]:
    x = med(d, l, V, key); put(nm, sec(x)) if x is not None else pending(nm)
x = med("edge_store", "bmgc_iter", 229, "seconds_per_iter")
if x is not None: put("scBmgcProjLoH", f"{x*2000/3600:.0f}"); put("scBmgcProjHiH", f"{x*10000/3600:.0f}")
else: pending("scBmgcProjLoH"); pending("scBmgcProjHiH")
x = med("edge_store", "ard_iter", 100, "seconds_per_iter")
if x is not None: put("scArdProjH", f1(x * 500 / 3600))
else: pending("scArdProjH")
mem = [fl(r["peak_mb"]) for r in (sc or []) if r["design"] == "two_sided" and r["learner"] == "ard_logit_K10" and int(r["V"]) == 229 and fl(r["peak_mb"]) is not None]
put("scMemTwoB", f"{st.median(mem):.0f}") if mem else pending("scMemTwoB")
mem = [fl(r["peak_mb"]) for r in (sc or []) if r["design"] == "edge_store" and r["learner"] == "ridge_dual_cv" and int(r["V"]) == 229 and fl(r["peak_mb"]) is not None]
put("scMemEdgeB", f"{st.median(mem):.0f}") if mem else pending("scMemEdgeB")

# ---- E1-derived cost macros: measured values and clearly labeled projections (no hand-typed timings)
def medv(design, learner, V, key="seconds"): return med(design, learner, V, key)
def putm(name, x, fmt=sec):
    put(name, fmt(x)) if x is not None else pending(name)
for nm, d, l, V, key in [("scSsvsIterC", "edge_store", "ssvs_iter", 100, "seconds_per_iter"), ("scSsvsIterB", "edge_store", "ssvs_iter", 229, "seconds_per_iter"),
                          ("scBmgcIterD", "edge_store", "bmgc_iter", 200, "seconds_per_iter"), ("scOneC", "one_sided", "ridge_dual_cv_R10", 100, "seconds"),
                          ("scEdgeRidgeC", "edge_store", "ridge_dual_cv", 100, "seconds"), ("scLassoC", "edge_store", "lasso_glmnet", 100, "seconds")]:
    putm(nm, medv(d, l, V, key))
x = medv("edge_store", "ssvs_iter", 100, "seconds_per_iter")
putm("scSsvsProjC", x * 1200 if x is not None else None)      # PROJECTED: per-iteration x 1,200 iterations
x = medv("edge_store", "ssvs_iter", 229, "seconds_per_iter")
putm("scSsvsProjB", x * 1200 if x is not None else None)
x = medv("edge_store", "ard_iter", 100, "seconds_per_iter")
putm("scArdProjIterPowMin", x * cub / 60 if x is not None else None, f1)   # PROJECTED: cubic in p from the measured V=100 iteration
putm("scArdProjPowH", x * cub * 500 / 3600 if x is not None else None, f1)   # PROJECTED: 500 iterations
bad = sorted(int(r["V"]) for r in (sc or []) if r["learner"] == "ard_iter" and r["status"] == "infeasible_memory")
put("scArdInfeasV", bad[0]) if bad else pending("scArdInfeasV")
if sc:
    bt = sorted(set(r["blas_threads"] for r in sc)); put("scBlasOne", "1" if all("OPENBLAS_NUM_THREADS=1" in b for b in bt) else "0")
# platform of the 2026-11-11 rerun (system_profiler, sw_vers, sessionInfo; see results/icde/LOG.md)
put("platformText", r"Apple Mac mini (M6, 12 cores, 32 GB RAM), macOS 27.2, R 4.6.1 with its reference BLAS and LAPACK, one thread")
# ---- additional E1 macros and pgfplots data for the scaling figure
for nm, d, l, V in [("scOneD", "one_sided", "ridge_dual_cv_R10", 400), ("scOneA", "one_sided", "ridge_dual_cv_R10", 50)]:
    putm(nm, medv(d, l, V))
Vs = sorted(set(int(r["V"]) for r in (sc or [])))
tw = [medv("two_sided", "ard_logit_K10", V_) for V_ in Vs]
if sc and all(x is not None for x in tw): put("scTwoMin", sec(min(tw))); put("scTwoMax", sec(max(tw)))
else: pending("scTwoMin"); pending("scTwoMax")
a_, b_, c_ = medv("edge_store", "ridge_dual_cv", 229), medv("one_sided", "ridge_dual_cv_R10", 229), medv("two_sided", "ard_logit_K10", 229)
if None not in (a_, b_, c_): put("ratioEdgeOneFitB", f1(a_ / b_)); put("ratioOneTwoFitB", f1(b_ / c_))
else: pending("ratioEdgeOneFitB"); pending("ratioOneTwoFitB")
if sc:
    def wdat(name, rows, hdr="V sec"):
        with open(os.path.join(HERE, "data", name), "w") as f:
            f.write(hdr + "\n")
            for V_, x_ in rows: f.write(f"{V_} {x_:.6g}\n")
    def series(d, l, key="seconds", mult=1.0):
        return [(V_, medv(d, l, V_, key) * mult) for V_ in Vs if medv(d, l, V_, key) is not None]
    wdat("sc_two.dat", series("two_sided", "ard_logit_K10")); wdat("sc_one.dat", series("one_sided", "ridge_dual_cv_R10"))
    wdat("sc_edge_ridge.dat", series("edge_store", "ridge_dual_cv")); wdat("sc_edge_lasso.dat", series("edge_store", "lasso_glmnet"))
    wdat("sc_ssvs_proj.dat", series("edge_store", "ssvs_iter", "seconds_per_iter", 1200))   # PROJECTED: 1,200 iterations
    wdat("sc_bmgc_proj.dat", series("edge_store", "bmgc_iter", "seconds_per_iter", 2000))   # PROJECTED: 2,000 iterations
    wdat("sc_bmgc_proj_hi.dat", series("edge_store", "bmgc_iter", "seconds_per_iter", 10000))   # PROJECTED: 10,000 iterations
    x100 = medv("edge_store", "ard_iter", 100, "seconds_per_iter")
    if x100 is not None:
        pS_ = 100 * 99 / 2
        rows = [(V_, x100 * 500 * ((V_ * (V_ - 1) / 2) / pS_) ** 3) for V_ in Vs if V_ >= 100]   # PROJECTED: cubic in p, 500 iterations
        wdat("sc_ard_proj.dat", rows)
        wdat("sc_ard_stop.dat", [(r_[0], r_[1]) for r_ in rows if r_[0] in bad])            # cells not run (memory)
# break-even number of fits m at V=229: template + synopsis build + m two-sided fits vs m dual edge-ridge fits
tm, bd, tw, rg = (medv("template", "mean_plus_eigen", 229), medv("two_sided", "build", 229),
                  medv("two_sided", "ard_logit_K10", 229), medv("edge_store", "ridge_dual_cv", 229))
if None in (tm, bd, tw, rg): pending("breakEvenText")
elif rg > tw:
    mstar_ = int(-(-(tm + bd) // (rg - tw)))
    put("breakEvenText", "With the medians of Fig.~\\ref{fig:scaling} at $V=229$, the template and synopsis build (\\scTemplateB{} s plus \\scBuildTwoB{} s) are repaid, relative to dual edge ridge, within $m^\\ast=%d$ fit%s." % (mstar_, "s" if mstar_ > 1 else ""))
else:
    put("breakEvenText", "With the medians of Fig.~\\ref{fig:scaling} at $V=229$, a dual-form edge ridge fit (\\scEdgeRidgeB{} s) is not slower than a two-sided ARD fit (\\scTwoB{} s), so no break-even number of fits exists on compute alone and the advantage of the synopsis is storage and memory.")

# ------------------------------------------------------------------ E3 data (read first: E2 uses n_train)
# Ridge arms of the binary task use inner folds by subject (k_by_n_grouped.csv, built by
# Codes/icde/make_k_by_n_grouped.R); fall back to the original grid if it is absent.
kb = rd(os.path.join(IRES, "k_by_n_grouped.csv")) or rd(os.path.join(IRES, "k_by_n.csv"))
def kbc(atlas, oc, design, KR=None, frac=1.0):
    """{split: metric} for one cell of the E3 grid (edge rows: design edge_ridge_cv only)."""
    if not kb: return {}
    out = {}
    for r in kb:
        if (r["atlas"] == atlas and r["outcome"] == oc and r["design"] == design and fl(r["metric"]) is not None
                and abs(float(r["train_frac"]) - frac) < 1e-9
                and (KR is None or (r["K_or_R"] != "NA" and int(float(r["K_or_R"])) == KR))):
            out[int(r["split"])] = fl(r["metric"])
    return out
def kbn(atlas, oc, frac):
    xs = [int(r["n_train"]) for r in (kb or []) if r["atlas"] == atlas and r["outcome"] == oc and abs(float(r["train_frac"]) - frac) < 1e-9]
    return st.mean(xs) if xs else None
ATL = [("Sch", "sch100"), ("Pow", "power229")]
TSK = [("T", "nbk_0v1"), ("A", "age_continuous")]
FRS = [0.25, 0.5, 0.75, 1.0]
KS_T = [5, 10, 15, 20, 30]
dfMaxA = int(kbn("power229", "age_continuous", 1.0)) if kb else None
if dfMaxA is not None: put("dfMaxA", dfMaxA)
else: pending("dfMaxA")

# ------------------------------------------------------------------ E2: df at Power-229 and Schaefer-100
e2 = (rd(os.path.join(IRES, "df_power229.csv")) or []) + (rd(os.path.join(IRES, "df_sch100.csv")) or [])
e2 = e2 or None
def e2col(atlas, oc, method, key):
    return col(e2, lambda r: r["atlas"] == atlas and r["outcome"] == oc and r["method"] == method, key) if e2 else {}
for atlas, oc, tag, ours in [("power229", "nbk_0v1", "PowT", "gs_ard_logit"), ("power229", "age_continuous", "PowA", "gs_ard"),
                              ("sch100", "age_continuous", "SchA", "gs_ard")]:
    o, r_, k = e2col(atlas, oc, ours, "eff_df"), e2col(atlas, oc, "raw_edge_ridge", "eff_df"), e2col(atlas, oc, "kundu_kernel_ridge", "eff_df")
    for nm, dd in [("Ours", o), ("Ridge", r_), ("Kr", k)]:
        if dd: put(f"df{tag}{nm}", f2(mean(dd))); put(f"df{tag}{nm}Sd", f2(sd(dd)) if len(dd) > 1 else "0")
        else: pending(f"df{tag}{nm}"); pending(f"df{tag}{nm}Sd")
    for nm, mth in [("Ours", ours), ("Ridge", "raw_edge_ridge")]:
        mm = e2col(atlas, oc, mth, "metric")
        put(f"fitMetric{tag}{nm}", f3(mean(mm))) if mm else pending(f"fitMetric{tag}{nm}")
    if o and r_:
        d, w, n, p = paired(o, r_); put(f"df{tag}Diff", f"{d:+.2f}"); put(f"df{tag}P", pv(p))
        put(f"df{tag}Ratio", f"{mean(r_)/mean(o):.2f}")
    else: pending(f"df{tag}Diff"); pending(f"df{tag}P"); pending(f"df{tag}Ratio")
    if e2 and r_:
        cs = col(e2, lambda r: r["atlas"] == atlas and r["outcome"] == oc and r["method"] == "raw_edge_ridge" and fl(r["pred_sd"]) is not None and fl(r["pred_sd"]) < 1e-8, "split")
        put(f"nConst{tag}", len(cs))
    else: pending(f"nConst{tag}")
    if oc == "age_continuous":
        # ridge df at the ceiling n_train: boundary of the penalty grid (caveat in capacity section and threat (x))
        if r_ and dfMaxA is not None:
            put(f"dfBound{tag}", sum(1 for v_ in r_.values() if v_ >= dfMaxA - 0.01))
        else: pending(f"dfBound{tag}")
if e2:
    for tag, m in [("ours", "gs_ard_logit"), ("ours", "gs_ard"), ("ridge", "raw_edge_ridge"), ("kr", "kundu_kernel_ridge")]:
        for atl, at in [("pow", "power229"), ("sch", "sch100")]:
            for oc in sorted(set(r["outcome"] for r in e2)):
                rows = [r for r in e2 if r["atlas"] == at and r["method"] == m and r["outcome"] == oc]
                if rows and not (at == "sch100" and oc == "nbk_0v1"):
                    with open(os.path.join(HERE, "data", f"df_{atl}_{tag}_{oc}.dat"), "w") as f:
                        f.write("nominal df\n")
                        for r in rows: f.write(f"{r['nominal_dim']} {fl(r['eff_df']):.4f}\n")

# ------------------------------------------------------------------ E3: synopsis size x training size
def sdisp(c):
    """display label of a configuration key"""
    d, k = c
    return (f"T, $K{{=}}{k}$" if d == "T" else f"O, $R{{=}}{k}$" if d == "O" else "E")
def cfgs(atlas, oc, frac):
    """{(design letter, size): {split: metric}} for the designs reported in the paper.
    Edge ridge enters on T2 only (T1 edge arms belong to the concurrent paper); glmnet-wide is never used."""
    g = {}
    for k_ in KS_T: g[("T", k_)] = kbc(atlas, oc, "two_sided_ard", k_, frac)
    for r_ in (5, 10): g[("O", r_)] = kbc(atlas, oc, "one_sided_ridge", r_, frac)
    if oc == "age_continuous": g[("E", 0)] = kbc(atlas, oc, "edge_ridge_cv", None, frac)
    return {k_: v_ for k_, v_ in g.items() if v_}
def rng(xs):
    xs = sorted(set(int(round(x)) for x in xs))
    return str(xs[0]) if len(xs) == 1 else f"{xs[0]} to {xs[-1]}"
tab = []
if kb:
    atlases = [("sch100", "$V{=}100$"), ("power229", "$V{=}229$")]
    outs = [("nbk_0v1", "State (AUC)"), ("age_continuous", "Age (Spearman)")]
    for at, atl in atlases:
        for oc, ocl in outs:
            cells = []
            for f_ in FRS:
                groups = cfgs(at, oc, f_)
                if len(groups) < 2: cells.append("n/a"); continue
                ranked = sorted(groups.items(), key=lambda kv: -st.mean(kv[1].values()))
                best, runner = ranked[0], ranked[1]
                d, w, n, p = paired(best[1], runner[1])
                # dagger: 95% interval of the corrected resampled t (Nadeau and Bengio) excludes zero
                dd_ = [best[1][s_] - runner[1][s_] for s_ in best[1] if s_ in runner[1]]
                se_ = ((1.0 / len(dd_) + 0.30 / 0.70) * st.variance(dd_)) ** 0.5
                dag = r"$^\dagger$" if len(dd_) == 10 and abs(st.mean(dd_)) > 2.262157 * se_ else ""
                line1 = f"{sdisp(best[0])} (${d:+.3f}$){dag}"
                if best[0][0] == "O":
                    tk = max((kv for kv in groups.items() if kv[0][0] == "T"), key=lambda kv: st.mean(kv[1].values()))[0]
                    line1 += r"\newline " + f"best T: $K{{=}}{tk[1]}$"
                cells.append(r"\begin{tabular}[c]{@{}l@{}}" + line1.replace(r"\newline ", r"\\ ") + r"\end{tabular}")
            tab.append(f"{atl} & {ocl} & " + " & ".join(cells) + r" \\")
    ntrow = {}
    for oc, lab in (("nbk_0v1", "State, training graphs"), ("age_continuous", "Age, training subjects")):
        ntrow[oc] = r"\multicolumn{2}{@{}l}{" + lab + "} & " + " & ".join(rng([kbn(a_, oc, f_) for _, a_ in ATL]) for f_ in FRS) + r" \\"
    with open(os.path.join(HERE, "tab_guidance_data.tex"), "w") as f:
        f.write("\\begin{tabular}{@{}ll" + "l" * len(FRS) + "@{}}\n\\toprule\n")
        f.write("Resolution & Task & " + " & ".join(f"{int(round(100*f_))}\\% of training data" for f_ in FRS) + " \\\\\n\\midrule\n")
        f.write("\n".join(tab) + "\n\\midrule\n" + ntrow["nbk_0v1"] + "\n" + ntrow["age_continuous"] + "\n\\bottomrule\n\\end{tabular}\n")
    put("guidanceReady", "1")
else:
    pending("guidanceReady")
    if os.path.exists(os.path.join(HERE, "tab_guidance_data.tex")): os.remove(os.path.join(HERE, "tab_guidance_data.tex"))

if kb:
    # one-sided R=10 level (Table II row, dotted references)
    for an, at in ATL:
        for tn, oc in TSK:
            m_ = kbc(at, oc, "one_sided_ridge", 10); put(f"kbOne{an}{tn}", f3(mean(m_)))
        put(f"kbEdge{an}A", f3(mean(kbc(at, "age_continuous", "edge_ridge_cv"))))
    # one-sided vs edge ridge on age, and two-sided vs edge ridge on age
    cells_oe = []
    for an, at in ATL:
        E_ = kbc(at, "age_continuous", "edge_ridge_cv")
        d, w, n, p = paired(kbc(at, "age_continuous", "one_sided_ridge", 10), E_)
        put(f"kbOneEdge{an}A", f"{d:+.3f}"); put(f"kbOneEdge{an}AP", pv(p)); put(f"kbOneEdge{an}AW", w)
        d, w, n, p = paired(kbc(at, "age_continuous", "two_sided_ard", 10), E_)
        put(f"kbTwoEdge{an}A", f"{d:+.3f}"); put(f"kbTwoEdge{an}AAbs", f"{abs(d):.3f}"); put(f"kbTwoEdge{an}AP", pv(p)); put(f"kbTwoEdge{an}AW", w)
        d, w, n, p = paired(kbc(at, "age_continuous", "one_sided_ridge", 10), kbc(at, "age_continuous", "two_sided_ard", 10))
        put(f"kbOneTwo{an}A", f"{d:+.3f}"); put(f"kbOneTwo{an}AP", pv(p)); put(f"kbOneTwo{an}AW", w)
        g25 = {k_: kbc(at, "age_continuous", "two_sided_ard", k_, 0.25) for k_ in KS_T}
        kbest = max(g25, key=lambda k_: mean(g25[k_]))
        d, w, n, p = paired(g25[kbest], kbc(at, "age_continuous", "edge_ridge_cv", None, 0.25))
        put(f"kbTwoEdgeLow{an}A", f"{d:+.3f}"); put(f"kbTwoEdgeLow{an}AAbs", f"{abs(d):.3f}")
        for f_ in FRS:
            d, w, n, p = paired(kbc(at, "age_continuous", "one_sided_ridge", 10, f_), kbc(at, "age_continuous", "edge_ridge_cv", None, f_))
            cells_oe.append(d)
    put("kbOneEdgeMinA", f"{min(cells_oe):+.3f}"); put("kbOneEdgeMaxA", f"{max(cells_oe):+.3f}")
    # weak contrast: best one-sided vs best two-sided, T1, 8 cells
    ns = 0; ncell = 0
    for an, at in ATL:
        for f_ in FRS:
            gg = cfgs(at, "nbk_0v1", f_)
            bo = max((kv for kv in gg.items() if kv[0][0] == "O"), key=lambda kv: mean(kv[1]))[1]
            bt = max((kv for kv in gg.items() if kv[0][0] == "T"), key=lambda kv: mean(kv[1]))[1]
            d, w, n, p = paired(bo, bt); ncell += 1; ns += (p >= 0.05)
    put("kbOneTwoTns", ns); put("kbOneTwoTn", ncell)
    # best two-sided K over the 16 cells; larger K at full training
    small = 0; ncell = 0; bign = lower = sig = 0
    for an, at in ATL:
        for tn, oc in TSK:
            for f_ in FRS:
                gg = {k_: kbc(at, oc, "two_sided_ard", k_, f_) for k_ in KS_T}
                kbest = max(gg, key=lambda k_: mean(gg[k_])); ncell += 1; small += (kbest in (5, 10))
            for k_ in (15, 20, 30):
                d, w, n, p = paired(kbc(at, oc, "two_sided_ard", k_), kbc(at, oc, "two_sided_ard", 10))
                bign += 1; lower += (d < 0); sig += (p < 0.05)
    put("kbBestKSmall", small); put("kbBestKN", ncell); put("kbKbigN", bign); put("kbKbigLower", lower); put("kbKbigSig", sig)
    def nlist(oc, fmt_last_range=True):
        xs = [rng([kbn(a_, oc, f_) for _, a_ in ATL]) for f_ in FRS]
        return ", ".join(xs[:-1]) + ", and " + xs[-1]
    put("kbNtrainA", nlist("age_continuous")); put("kbNtrainT", nlist("nbk_0v1"))
    # data files: size sweep at full training, and the training-size curves
    for an, at in ATL:
        for tn, oc in TSK:
            with open(os.path.join(HERE, "data", f"kb_ksweep_{an}_{tn}.dat"), "w") as f:
                f.write("K metric\n")
                for k_ in KS_T:
                    for s_, v_ in sorted(kbc(at, oc, "two_sided_ard", k_).items()): f.write(f"{k_} {v_:.4f}\n")
            with open(os.path.join(HERE, "data", f"kb_ksweep_{an}_{tn}_mean.dat"), "w") as f:
                f.write("K metric\n")
                for k_ in KS_T: f.write(f"{k_} {mean(kbc(at, oc, 'two_sided_ard', k_)):.4f}\n")
            arms = [("T5", "two_sided_ard", 5), ("T10", "two_sided_ard", 10), ("T30", "two_sided_ard", 30),
                    ("O5", "one_sided_ridge", 5), ("O10", "one_sided_ridge", 10)]
            if oc == "age_continuous": arms.append(("E", "edge_ridge_cv", None))
            for lab, dsg, k_ in arms:
                with open(os.path.join(HERE, "data", f"kbyn_{an}_{tn}_{lab}.dat"), "w") as f:
                    f.write("ntrain metric\n")
                    for f_ in FRS: f.write(f"{kbn(at, oc, f_):.1f} {mean(kbc(at, oc, dsg, k_, f_)):.4f}\n")
    put("ratioEdgeOnePow", f1(8 * pP / (8 * 229 * 10)))
else:
    for nm in ("kbOneSchT", "kbOnePowT", "kbOneSchA", "kbOnePowA", "kbOneEdgeMinA", "kbOneEdgeMaxA", "kbOneTwoTns", "kbBestKSmall",
               "kbKbigN", "kbNtrainA", "kbNtrainT"): pending(nm)

# ------------------------------------------------------------------ E4: K sweep and vertex identity at V=229 (filled only when all 10 splits are present)
k4 = rd(os.path.join(IRES, "ksweep_power229.csv"))
def e4(oc, var, key):
    return col(k4, lambda r: r["outcome"] == oc and r["variant"] == var and r["status"] == "ok", key) if k4 else {}
for tn, oc in TSK:
    n_ok = len(e4(oc, "current", "auc" if tn == "T" else "spearman"))
    put(f"eFourN{tn}", n_ok)
def e4complete(oc, var):
    key = "auc" if oc == "nbk_0v1" else "spearman"
    return k4 is not None and len(e4(oc, var, key)) == 10 and len(e4(oc, "current", key)) == 10
for nm, var in [("KTen", "current"), ("KTwenty", "bigK20"), ("KThirty", "bigK30"), ("KForty", "bigK40"), ("OneArd", "onesided")]:
    if e4complete("nbk_0v1", var): put(f"aucPow{nm}", f3(mean(e4("nbk_0v1", var, "auc"))))
    else: pending(f"aucPow{nm}")
for nm, var in [("KTen", "current"), ("OneArd", "onesided")]:
    if e4complete("age_continuous", var): put(f"spPow{nm}", f3(mean(e4("age_continuous", var, "spearman"))))
    else: pending(f"spPow{nm}")
if e4complete("nbk_0v1", "onesided"):
    a_, b_ = e4("nbk_0v1", "onesided", "seconds"), e4("nbk_0v1", "current", "seconds")
    r_ = st.median([a_[s_] / b_[s_] for s_ in a_ if s_ in b_])
    r_ = float(f"{r_:.2g}")                                   # two significant figures (within-run ratio)
    put("ratioOnePowT", f"{r_:,.0f}".replace(",", "{,}"))
    d, w, n, p = paired(e4("nbk_0v1", "onesided", "auc"), e4("nbk_0v1", "current", "auc"))
    put("dPowTOneArd", f"{d:+.3f}"); put("pPowTOneArd", pv(p)); put("wPowTOneArd", w)
else:
    for nm in ("ratioOnePowT", "dPowTOneArd", "pPowTOneArd", "wPowTOneArd"): pending(nm)
if e4complete("age_continuous", "onesided"):
    d, w, n, p = paired(e4("age_continuous", "onesided", "spearman"), e4("age_continuous", "current", "spearman"))
    put("dPowAOneArd", f"{d:+.3f}"); put("pPowAOneArd", pv(p)); put("wPowAOneArd", w)
else:
    for nm in ("dPowAOneArd", "pPowAOneArd", "wPowAOneArd"): pending(nm)

# ------------------------------------------------------------------ extra collections (age, existing results)
for tag, fn in (("Abide", "abide_repeated.csv"), ("Acpi", "acpi_age_repeated.csv")):
    rows = rd(os.path.join(RES, fn))
    if not rows:
        for k in ("nSplit", "ArdM", "ArdS", "RidgeM", "RidgeS", "LassoM", "LassoS", "BaseM", "BaseS", "RidgeWins", "pRidge"): pending("x" + tag + k)
        continue
    g = {k: [fl(r[k]) for r in rows] for k in ("gs_ard", "ridge", "lasso", "baseline")}
    put("x" + tag + "nSplit", len(rows))
    for k, nm in (("gs_ard", "Ard"), ("ridge", "Ridge"), ("lasso", "Lasso"), ("baseline", "Base")):
        put("x" + tag + nm + "M", f2(st.mean(g[k]))); put("x" + tag + nm + "S", f2(st.stdev(g[k])))
    put("x" + tag + "RidgeWins", sum(a < b for a, b in zip(g["ridge"], g["gs_ard"])))
    put("x" + tag + "pRidge", pv(wilcoxon_exact([a - b for a, b in zip(g["gs_ard"], g["ridge"])])))
put("xV", 200); put("xP", commas(200 * 199 / 2)); put("xQ", 55)

# ------------------------------------------------------------------ write numbers.tex
with open(os.path.join(HERE, "numbers.tex"), "w") as f:
    f.write("% Generated by make_numbers.py from the result CSVs. Do not edit by hand.\n")
    f.write("\\providecommand{\\pending}[1]{\\textcolor{red}{[#1]}}\n")
    for k in sorted(M):
        f.write("\\newcommand{\\%s}{%s}\n" % (k, M[k]))
print(f"wrote numbers.tex with {len(M)} macros; pending: {sorted(k for k,v in M.items() if 'pending' in v)}")
