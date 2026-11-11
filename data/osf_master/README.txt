DATASET: OSF Power-229 master assembly (all connectomes + all derived outcome variables)
=========================================================================================

WHAT THIS FOLDER IS
The single assembled dataset that every OSF-derived analysis in this project reads: all 1003
connectomes from data/osf_guha, joined with age and every outcome variable derived from
condition/domain/load, in one R object (master.rds) plus a plain-CSV mirror of the metadata
(subjects.csv). This is what Codes/analysis/outcome_explorer.R, mode_network_mapping.R, and
rieck_derived_measures.R all actually load.

EXACT SOURCE
  Built entirely from files already in this project -- no separate network download:
    Connectomes: data/osf_guha/cache/*.npy (see data/osf_guha/README.txt for how those were
    fetched from OSF project 5zyws, https://osf.io/5zyws/).
    Ages: data/osf_guha/participant_ages.csv (from OSF component dm23b,
    https://osf.io/download/u3xpd/).

HOW TO (RE)BUILD
  Rscript Codes/fetch/build_osf_master.R
  Reads data/osf_guha/subjects.csv + cache/, joins participant_ages.csv, derives DOMAIN/LOAD/
  AGE_BIN/AGE_Q/IS_TASK (rules below), and writes master.rds, subjects.csv, and meta.csv into
  this folder. No saved script for this step could be found under Codes/ when it was first
  needed (it had evidently been run ad hoc in an earlier session); it was reconstructed by
  reverse-engineering the existing master.rds/subjects.csv schema and verified to reproduce
  every metadata column exactly and the connectome values to about 5e-7 (float32 storage
  precision, not a real discrepancy).

FILES IN THIS FOLDER
  master.rds  -- an R list with three elements: conn (a list of 1003 229x229 numeric
                 matrices, one per row of subjects.csv, in the same order), subjects (a data
                 frame, 1003 rows x 9 columns -- identical content to subjects.csv below), and
                 V (integer, 229). Load with readRDS("data/osf_master/master.rds"); slice
                 connectome i with d$conn[[i]] (1-indexed in R).
  subjects.csv  -- the same 1003 x 9 metadata table as master.rds$subjects, in CSV form (see
                    dictionary below). There is no connectomes_full.csv in this folder --
                    connectome data lives only inside master.rds (as a list of matrices, not
                    a stacked CSV) and in data/osf_guha/connectomes_full.csv.
  meta.csv  -- one header + one data row: N=1003, V=229, atlas=Power229, pipeline=osf_master,
               strategy=all_conditions.

VARIABLE DICTIONARY FOR subjects.csv (all 9 columns, values confirmed by reading the file)
  FILE_ID   string, e.g. Power229_10_sub-001_gng_inhibition (no .txt extension, unlike
            data/osf_guha/subjects.csv's FILE_ID).
  SUBJECT   integer, the bare subject number (e.g. 1 for sub-001). 144 unique values.
  COND      string, condition name. Confirmed 7 distinct values: gng_inhibition,
            gng_initiation, nbk_0back, nbk_1back, nbk_2back, rest, tsw_shifting.
  DOMAIN    string, the cognitive-control domain the condition belongs to. Confirmed mapping
            (each COND maps to exactly one DOMAIN): gng_inhibition -> inhibition;
            gng_initiation -> initiation; nbk_0back, nbk_1back, nbk_2back -> wm (working
            memory); tsw_shifting -> shifting; rest -> rest. Note DOMAIN therefore has FIVE
            distinct values (inhibition, initiation, wm, shifting, rest), not four -- rest
            gets its own domain level rather than being folded into one of the four
            cognitive-control domains.
  LOAD      integer or NA. Confirmed: defined only for the three n-back rows
            (nbk_0back -> 0, nbk_1back -> 1, nbk_2back -> 2); NA for every other condition
            (rest, both go/no-go conditions, task-switching). Treat as a 3-level ordinal
            experimental factor, not a general continuous covariate -- it is only ever 0, 1,
            2, or NA.
  AGE       integer, years, joined from participant_ages.csv by SUBJECT. Confirmed range in
            this sample: 20 to 86. Constant within subject (same value on every row for a
            given SUBJECT).
  AGE_BIN   integer, 0 or 1. Confirmed: median split at 47.5 (computed once per subject, not
            per row); 1 means "older" (AGE >= 48), 0 means "younger" (AGE <= 47). Confirmed
            exactly balanced at the subject level: 72 subjects coded 1, 72 coded 0 (mirrored
            across each subject's rows, so the per-row counts are close to but not exactly
            50/50 -- they follow each subject's own row count).
  AGE_Q     integer, 1 to 4. Confirmed: quartile bin of AGE, 1 = youngest quartile, 4 =
            oldest quartile. Per-subject counts: 39 subjects in quartile 1, 33 in quartile 2,
            36 in quartile 3, 36 in quartile 4 (quartile boundaries are not perfectly even
            because AGE has repeated values near the cut points).
  IS_TASK   integer, 0 or 1. Confirmed: 0 only for COND=rest, 1 for every task condition
            (all six non-rest conditions).

COUNTS
  1003 rows, 144 distinct subjects, V=229 (Power atlas). Edge dimension
  p = V(V-1)/2 = 26,106. Conditions-per-subject: 139 subjects have all 7 conditions, 5
  subjects have 6 (none have fewer). Per-condition row counts: nbk_0back 144, nbk_1back 144,
  nbk_2back 144, gng_inhibition 144, gng_initiation 143, tsw_shifting 143, rest 141.

KNOWN CAVEATS
  AGE is the only participant-level variable available anywhere in this dataset -- there is
  no sex, education, MMSE, motion, or behavioral-score column, confirmed by checking every
  file OSF actually released for this project (see data/osf_guha/README.txt). Any analysis
  claiming to use a different demographic on this data is working from something not present
  in the source.
  LOAD encodes the n-back task's EXPERIMENTALLY ASSIGNED load level (0/1/2), not a measured
  behavioral response (e.g. not accuracy or reaction time) -- treating it as a continuous
  regression target assumes the 0->1 and 1->2 gaps are equal in effect size, which is a
  modeling choice, not a property confirmed in the data.
  This master file was rebuilt on 2026-09-10 against the corrected, complete 1003-file
  osf_guha cache. Any results file elsewhere in this project without a "_v1003" suffix in its
  filename was computed against an earlier, incomplete 958-file version of this master table
  and should be treated as stale for OSF-derived analyses (see results/ and
  data/_archive_v958/README.txt).

CITE: Rieck et al. 2021 (Neurobiology of Aging, source data collection); Rodriguez-Acosta,
Guha & Dinov 2024 (BMGC, the competitor benchmark that motivated using this dataset).
