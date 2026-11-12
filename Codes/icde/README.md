# Codes/icde: experiments and number generation

Purpose: measured cost versus number of vertices (E1), effective degrees of freedom at 229 vertices (E2),
synopsis size by training-set size grid (E3), optional synopsis-size sweep at 229 vertices (E4), and figures.
Run everything from the repo root. Nothing under `Codes/core`, `Codes/analysis`, or existing `results/*.csv` is modified;
all outputs go to `results/icde/`.

## Inputs
`data/osf_master/master.rds` (Power-229), `data/osf_schaefer100/master.rds` (Schaefer-100), existing
`results/model_df_comparison.csv` and `results/dense_variants_sch100_combined.csv` (figures only).
R packages: glmnet, igraph. E1 BMGC cell also needs gdata, pROC, mvtnorm, pscl, MCMCpack, LaplacesDemon, GIGrvg, BayesLogit, MASS.

## Files
- `icde_common.R`: shared helpers (data loading, splits seed 900+s, template/synopsis features, df estimators copied from `model_df_comparison.R`).
- `E1_scaling.R` -> `results/icde/scaling.csv` (V,p,design,learner,repeat,seconds,seconds_per_iter,n_iter_timed,peak_mb,bytes_per_graph,status,measured,blas_threads,note).
  Status in ok / timeout / infeasible_memory / error. Edge ARD is skipped (recorded, `measured=FALSE`, note says PROJECTED) when
  6*p*p*8 bytes exceed 70% of RAM (`--ram_gb=N` overrides detection) or when the p^3 projection of its 2 timed iterations from the
  previous V exceeds the 600 s cap (setTimeLimit cannot stop a LAPACK solve, so the guard must act before the call). A cell slower
  than 120 s is not repeated (single measurement). bytes_per_graph = 8 x number of stored features (float64 payload).
  Template row = mean |A| + eigendecomposition; synopsis build rows exclude it (one-sided reuses the template eigenvectors).
  An untimed warmup runs first; report the median of the 3 repeats. `peak_mb` is the R-heap increment (gc max used), not BLAS scratch or OS-level memory.
  BMGC: `bmgc_scaling_probe()`; `L <- 2` is hardcoded in `guha_bmgc.R`, harmless for one layer. No projections are computed here.
- `E2_df_power229.R` (E2) -> `results/icde/df_power229.csv`; with `--atlas=sch100 --outcomes=age_continuous` (E2b) ->
  `results/icde/df_sch100.csv`. Columns: atlas,outcome,split,method,nominal_dim,eff_df,metric_name,metric,seconds,pred_sd,penalty.
  nbk_0v1 (metric auc): gs_ard_logit, raw_edge_ridge (GLM ridge, deviance CV), kundu_kernel_ridge, estimators as in
  `model_df_comparison.R`. age_continuous (metric spearman, subject-mean connectomes): gs_ard (1 + MacKay sum gamma from
  ard_fit), raw_edge_ridge (raw_edge_ridge_cv; 1 + sum d/(d + alpha) over train Gram eigenvalues), kundu_kernel_ridge
  (1 + trace K(K + lambda I)^-1). On nbk_0v1 at V = 229 ridge collapses (df near 1, pred_sd near 0), so age gives the
  non-degenerate comparison at both resolutions. Schaefer-100 nbk_0v1 is not rerun (already in model_df_comparison.csv).
- `E3_k_by_n.R` -> `results/icde/k_by_n.csv` (atlas,outcome,design,K_or_R,train_frac,split,n_train,metric_name,metric,seconds).
  Resumable by (atlas, outcome, split). `--no-glmnet-ref` drops the logistic wide-path ridge reference.
- `E4_ksweep_power229.R` (optional) -> `results/icde/ksweep_power229.csv`. Reference arms as in dense_variants.R: ridge =
  glmnet logistic (AUC CV) for nbk_0v1 and Woodbury ridge alpha = 50 for age; kundu as in competitor_benchmarks.R. `--no-onesided` skips the very slow one-sided ARD fit.
- `make_icde_figures.R [outdir]` -> `fig_scaling.pdf`, `fig_df.pdf`, `fig_kbyn.pdf`, `fig_ksweep.pdf`, `kbyn_best_per_cell.csv`
  (default outdir `results/icde/figs/`). Dashed lines in fig_scaling are projections from measured per-iteration cost and are labeled so.
- `run_icde.sh [--smoke]`: E1, E2, E2b, E3 (and E4 if `RUN_E4=1`), then figures; log in `results/icde/LOG.md` with sessionInfo.

## Experiments G1 to G10
- `icde_common.R` also holds `subject_folds()` (inner CV folds by subject), `ridge_dual_cv_folds()` (dual-Gram CV ridge with
  given folds; same grid and refit as `raw_edge_ridge_cv`), and `pca_scores()` (edge-space PCA from training graphs only).
- `G1_generic_synopses.R --atlas=<sch100|power229> --outcome=<age_continuous|nbk_0v1>` -> `results/icde/generic_<atlas>_<outcome>.csv`
  (atlas,outcome,design,size,learner,draw,train_frac,split,n_train,n_features,metric_name,metric,rmse,build_seconds,fit_seconds).
  Designs two_sided (K = 10), one_sided (R = 10), edge, pca (55), rp (55, 210, 1000; 5 Gaussian draws); learners ard and ridge.
  Ridge uses inner folds by subject. Average the rp draws within split before any paired comparison.
- `G2_external_collections.R` -> `external_collections.csv`: the same arms on ABIDE and ACPI age (seed 100 + split).
- `G3_scaling_n.R [--nmax=3200]` -> `scaling_n.csv`, `scaling_n_env.txt`: time against the number of training graphs at V = 229 on
  resampled collections. Idle machine. Timings from this file are never combined with `scaling.csv` (different platform).
- `G5_grouped_cv_t1.R --atlas=<...>` -> `grouped_cv_<atlas>.csv`: ridge arms on nbk_0v1 with inner folds over graphs or subjects.
- `G6_df_grouped_t1.R --atlas=<...>` -> `df_grouped_<atlas>.csv`: effective df of logistic edge ridge and the size of edge lasso
  under both fold schemes, with two-sided ARD refitted on the same splits.
- `make_k_by_n_grouped.R` -> `k_by_n_grouped.csv`: the E3 grid with nbk_0v1 ridge rows replaced by the subject-fold values of G5.
- `make_g_numbers.R <paper_dir>` -> `numbers_g.tex`, `tab_generic.tex`, `tab_external.tex`, `data/g_*.dat` in the paper folder.
  Paired summaries: mean difference, wins, exact Wilcoxon p, and the corrected resampled-t interval (variance factor 1/J + 3/7).

Why the folds matter: on nbk_0v1 each subject has two graphs with opposite labels. `raw_edge_ridge_cv` and the glmnet calls with
`nfolds` form folds over graphs, which splits the pair and selects too much shrinkage for edge-level ridge and lasso. The E2 and
E3 ridge and lasso values on nbk_0v1 are affected and are superseded by G5 and G6; age rows (one graph per subject) are not.

## Threads
`run_icde.sh` exports OMP/OpenBLAS/vecLib/MKL thread counts from `BLAS_THREADS` (default 1) and logs them with
sessionInfo (which names the BLAS library). Every scaling.csv row carries the thread setting. Report it in the paper.

## Run
    bash Codes/icde/run_icde.sh --smoke     # a few minutes (E1 at V=30,60; figures skipped); writes results/icde/smoke_*.csv
    bash Codes/icde/run_icde.sh             # full run, machine idle for E1
    Rscript Codes/icde/E3_k_by_n.R          # any script alone

## Expected runtimes (estimates, not measured)
E1: 1 to 3 hours (edge ARD, BMGC and SSVS dominate; capped cells take 10 minutes each). E2: 5 to 15 minutes.
E3: 20 to 60 minutes (two-sided K=30 logistic ARD and the logistic ridge reference are the slow arms).
E4: several hours with one-sided ARD, about 1 hour without.
