#!/bin/bash
# Run the ICDE experiments in order from the repo root.
#   bash Codes/icde/run_icde.sh [--smoke]          (RUN_E4=1 to also run the optional E4)
# E1 timings are only valid on an otherwise idle machine: quit browsers, IDEs, training jobs.
cd "$(dirname "$0")/../.." || exit 1
RS=${RSCRIPT:-/usr/local/bin/Rscript}
# Fix the BLAS/OpenMP thread count so every design is timed under the same setting (default 1;
# override with BLAS_THREADS=8 bash ...). The value and the BLAS library are written to the log.
T=${BLAS_THREADS:-1}
export OMP_NUM_THREADS=$T OPENBLAS_NUM_THREADS=$T VECLIB_MAXIMUM_THREADS=$T MKL_NUM_THREADS=$T
ARG="$1"
mkdir -p results/icde
LOG=results/icde/LOG.md
stamp() { date '+%Y-%m-%d %H:%M:%S %Z'; }
{ echo; echo "## Run $(stamp) ${ARG}"; echo "host: $(uname -a)"; echo "threads: OMP=$OMP_NUM_THREADS OPENBLAS=$OPENBLAS_NUM_THREADS VECLIB=$VECLIB_MAXIMUM_THREADS MKL=$MKL_NUM_THREADS"; } >> "$LOG"
$RS -e 'cat(capture.output(sessionInfo()), sep="\n")' >> "$LOG" 2>&1
echo "WARNING: E1 needs an idle machine (CPU and RAM). Starting in 10 s; Ctrl-C to abort."; sleep 10
run() {   # run <name> <script> [args]
  name=$1; shift
  echo "- $(stamp) start $name" | tee -a "$LOG"
  $RS "$@" > "results/icde/${name}.out" 2>&1; rc=$?
  echo "- $(stamp) end $name (exit $rc)" | tee -a "$LOG"
}
run E1 Codes/icde/E1_scaling.R $ARG
run E2 Codes/icde/E2_df_power229.R $ARG
run E2b Codes/icde/E2_df_power229.R --atlas=sch100 --outcomes=age_continuous $ARG
run E3 Codes/icde/E3_k_by_n.R $ARG
if [ "$RUN_E4" = "1" ]; then run E4 Codes/icde/E4_ksweep_power229.R $ARG; fi
if [ "$ARG" = "--smoke" ]; then echo "smoke: figures skipped (they read the full-run CSVs)"; else run figs Codes/icde/make_icde_figures.R; fi
echo "Done. Log: $LOG"
