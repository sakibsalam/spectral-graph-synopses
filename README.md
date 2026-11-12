# Spectral graph synopses

Code and result files for benchmarking ways of storing a collection of dense weighted graphs on a shared
vertex set (for example, one brain connectome per subject) before fitting a predictive model: the full edge
store, one-sided and two-sided spectral synopses, and generic reductions (PCA, random projection).

## Layout

```
Codes/icde/        experiments, shared helpers, and number generation (see Codes/icde/README.md)
Codes/core/        model code: ARD, dual ridge with cross-validation, template bases
Codes/analysis/    earlier runs whose result files are also used
Codes/benchmarks/  kernel ridge re-implementation; wrapper for the BMGC sampler
Codes/simulation/  spectral compressibility across collections
Codes/fetch/       download and assembly scripts for the public data
paper/             make_numbers.py: builds tables and plot data from results/
results/           result CSVs (results/icde for the main experiments)
data/              not included; README files say where each collection comes from
environment/       R session information for the runs
```

Run all commands from the repository root.

## Data

No data are redistributed. The collections are public:

- Lifespan cognitive-control cohort (OSF project `m5crs`): `python3 Codes/fetch/osf_guha_fetch.py`, then
  `Rscript Codes/fetch/build_osf_master.R` (see `data/osf_master/README.txt`, `data/osf_schaefer100/README.txt`).
- ABIDE and ACPI (Preprocessed Connectomes Project): `Codes/fetch/abide_fetch.py`, `Codes/fetch/acpi_fetch.py`.
- OASIS-3 structural connectomes from braingraph.org.

## Running

```
bash Codes/icde/run_icde.sh --smoke      # quick check that everything runs
bash Codes/icde/run_icde.sh              # main experiments
```

Individual experiments and their outputs are described in `Codes/icde/README.md`. Tables and plot data are built with

```
python3 paper/make_numbers.py . results/icde paper/out
Rscript Codes/icde/make_g_numbers.R paper/out
```

## Environment

R 4.3 or newer with `glmnet` and `igraph`; some benchmark cells need further packages listed in
`Codes/icde/README.md`. Python 3 for `make_numbers.py` (standard library only); the fetch scripts need `numpy`.
Session information is in `environment/` and `results/icde/LOG.md`.

## Third-party code

The BMGC sampler of Guha and co-authors is not included. `Codes/benchmarks/guha_bmgc_benchmark.R` expects a
copy of the authors' `binary_multi_network_function_with_results.R` at `Codes/benchmarks/guha_bmgc.R`;
without it, the BMGC timing cell is skipped.
