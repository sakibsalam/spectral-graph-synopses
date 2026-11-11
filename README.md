# Spectral synopses for collections of dense weighted graphs: benchmark code and results

Code, result files, and number-generation scripts for the paper

> S. Salam and A. Banerjee. *Spectral Synopses for Learning over Collections of Dense Weighted Graphs:
> Capacity, Cost, and Vertex Identity.* Submitted to ICDE 2027 (Experiment, Analysis, and Benchmark track).

The paper compares ways of storing a collection of dense weighted graphs on a shared vertex set
(one brain connectome per subject) before fitting a predictive model:

| Design | Stored per graph | Features |
|---|---|---|
| Edge store (E) | upper triangle of `A` | `p = V(V-1)/2` |
| One-sided synopsis (O) | `A U_R`, one row per vertex | `V R` |
| Two-sided spectral synopsis (T) | `U_K' A U_K`, a `K x K` core | `q = K(K+1)/2` |
| PCA, random projection (references) | projection of the edge vector | `q` |

`U` holds leading eigenvectors of a template graph (mean absolute adjacency of the training graphs).
Each design is fitted with dual-form ridge; the reductions to 55 numbers are also fitted with
automatic relevance determination (ARD).

## Layout

```
Codes/icde/        experiments of the paper (E1 to E4, G1 to G6), shared helpers, number generation
Codes/core/        model code: ARD (Gaussian and logistic), dual ridge with CV, template bases
Codes/analysis/    earlier runs whose result files the paper also uses (synopsis-size sweep, operator study)
Codes/benchmarks/  kernel ridge re-implementation of a latent-node method; wrapper for the BMGC sampler
Codes/simulation/  spectral compressibility across collections
Codes/fetch/       download and assembly scripts for the public data
paper/             make_numbers.py: builds the LaTeX macros, table rows, and plot data from results/
results/           result CSVs (results/icde for the paper's experiments)
data/              not included; README files say where each collection comes from and how to fetch it
environment/       R session information for the runs
```

All scripts use paths relative to the repository root. Run every command from the root.

## Data

No data are redistributed here. The collections are public:

- Lifespan cognitive-control cohort (Rieck et al., 2021), Open Science Framework project `m5crs`:
  Power-229 and Schaefer-100 connectivity matrices. `python3 Codes/fetch/osf_guha_fetch.py`, then
  `Rscript Codes/fetch/build_osf_master.R` (see `data/osf_master/README.txt`, `data/osf_schaefer100/README.txt`).
- ABIDE and ACPI (Preprocessed Connectomes Project): `Codes/fetch/abide_fetch.py`, `Codes/fetch/acpi_fetch.py`.
- OASIS-3 structural connectomes from braingraph.org (compressibility table only).

## Experiments

| Script | What it measures | Output in `results/icde/` | Paper |
|---|---|---|---|
| `E1_scaling.R` | fit and build time, memory versus number of vertices | `scaling.csv` | Fig. 2, Table I |
| `E2_df_power229.R` | effective degrees of freedom (Age; State with graph-level folds, superseded by G6) | `df_power229.csv`, `df_sch100.csv` | Fig. 5 |
| `E3_k_by_n.R` | two-sided size by training-set size grid | `k_by_n.csv` | Table VI |
| `E4_ksweep_power229.R` | synopsis-size sweep and one-sided ARD at V = 229 (earlier Linux run, not repeated) | `ksweep_power229.csv` | Fig. 6 |
| `G1_generic_synopses.R` | design by learner: T, PCA, random projection, O, E with ridge and ARD | `generic_<atlas>_<task>.csv` | Table II, Fig. 4 |
| `G2_external_collections.R` | the Age comparison on ABIDE and ACPI | `external_collections.csv` | Table III |
| `G3_scaling_n.R` | fit time versus number of graphs (resampled collections) | `scaling_n.csv` | Fig. 3 |
| `G7_projection_posterior_pilot.R` | fit one posterior, cache it, answer sparse-summary, interval, and vertex queries from the cache | `projpost_accuracy.csv`, `projpost_queries.csv` | Sec. V-F, Table IV |
| `G8_projection_posterior_extended.R` | cached-posterior queries on four collections: sparse and mode-level summaries, predictive intervals | `projpost_ext_*.csv` | Sec. V-F, Table IV |
| `G9_ksweep_ridge.R` | two-sided size sweep with ridge | `ksweep_ridge.csv` | Sec. V-D |
| `G10_scaling_n_primal.R` | primal ridge on the synopses against the number of graphs | `scaling_n_primal.csv` | Fig. 3 |
| `G5_grouped_cv_t1.R` | inner folds over graphs versus over subjects, ridge arms on State | `grouped_cv_<atlas>.csv` | Sec. V-C |
| `G6_df_grouped_t1.R` | effective df on State with subject-level folds | `df_grouped_<atlas>.csv` | Sec. V-C, Fig. 5 |
| `make_k_by_n_grouped.R` | merges G5 into the E3 grid | `k_by_n_grouped.csv` | Table VI |
| `make_g_numbers.R`, `paper/make_numbers.py` | every number, table row, and plotted point | `numbers*.tex`, `tab_*.tex`, `data/*.dat` | all |

Protocol, in short: 10 subject-level splits (30% of subjects held out, seed 900 + split); template,
PCA loadings, and standardization from training graphs only; inner cross-validation folds formed over
subjects; paired comparisons summarized by the mean difference and a corrected resampled-t interval
(Nadeau and Bengio, 2003).

```
bash Codes/icde/run_icde.sh --smoke      # minutes; checks that everything runs
bash Codes/icde/run_icde.sh              # E1 to E3 (E1 needs an idle machine); RUN_E4=1 adds E4
for a in sch100 power229; do for o in age_continuous nbk_0v1; do
  Rscript Codes/icde/G1_generic_synopses.R --atlas=$a --outcome=$o; done; done
Rscript Codes/icde/G2_external_collections.R
Rscript Codes/icde/G3_scaling_n.R        # idle machine
for a in sch100 power229; do Rscript Codes/icde/G5_grouped_cv_t1.R --atlas=$a; Rscript Codes/icde/G6_df_grouped_t1.R --atlas=$a; done
Rscript Codes/icde/make_k_by_n_grouped.R
python3 paper/make_numbers.py . results/icde paper/out     # numbers.tex, tables, plot data
Rscript Codes/icde/make_g_numbers.R paper/out              # numbers_g.tex, tables, plot data
```

`make_numbers.py` also reads `data/osf_master/subjects.csv` (subject count, age range, number of conditions), which the fetch
step creates; without it those few macros are marked pending.

## Environment

R 4.3 or newer with `glmnet` and `igraph`; the BMGC timing cell of E1 also needs `gdata`, `pROC`,
`mvtnorm`, `pscl`, `MCMCpack`, `LaplacesDemon`, `GIGrvg`, `BayesLogit`. Python 3 (standard library only)
for `make_numbers.py`; the fetch scripts need `numpy`. All timings in the paper come from one machine (Apple Mac mini, M6, 32 GB, R 4.6.1 with the reference BLAS, one thread),
recorded in `results/icde/LOG.md`. `results/icde/earlier_linux_run/` keeps an earlier run of E1 to E4 on a Linux virtual machine; its accuracy
values agree with the current files to within 1e-15.

## Third-party code

The BMGC sampler of Guha and co-authors is not included. `Codes/benchmarks/guha_bmgc_benchmark.R` expects
a copy of the authors' `binary_multi_network_function_with_results.R`, with the two small changes described in
that file's header, at `Codes/benchmarks/guha_bmgc.R`. Without it, E1 skips the BMGC cell and runs everything else.
The kernel ridge arm is our own re-implementation of the latent-node approach of Ma, Kundu, and Stevens;
no public code exists for it.

## License

To be chosen by the authors before public release.
