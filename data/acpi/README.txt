DATASET: ACPI (Addiction Connectome Preprocessed Initiative) -- resting-state functional connectomes
=========================================================================================

WHAT THIS FOLDER IS
150 resting-state functional-connectivity matrices from the ACPI dataset (same fcp-indi S3
family as ABIDE), used as a second, independent functional dataset to check whether the
ABIDE age finding (dense raw-edge signal) replicates elsewhere. Fully open data, no
data-use agreement required.

EXACT SOURCE
  Aggregated phenotype CSV:
    https://s3.amazonaws.com/fcp-indi/data/Projects/ACPI/PhenotypicData/acpi_aggregated_phenotypic_data.csv
  Site phenotype files (cannabis-use group, not merged into subjects.csv here):
    https://s3.amazonaws.com/fcp-indi/data/Projects/ACPI/PhenotypicData/mta_1_phenotypic_data.csv
    (and nyu_1_, du_1_ variants)
  ROI time series (per subject-session), pattern:
    https://s3.amazonaws.com/fcp-indi/data/Projects/ACPI/OutputsTimeSeries/<site>/<ID>-session_1/ts_cc200_rois.csv
    (site in {mta_1, nyu_1}; atlas in {cc200, aal, ho}. Layout is ROIs x [metadata column +
    timepoints] -- the transpose of ABIDE's layout; the fetcher handles this.)

HOW TO (RE)DOWNLOAD / REBUILD
  python3 Codes/fetch/acpi_fetch.py --n 150 --atlas cc200
  (Note: the script now lives under Codes/fetch/, not directly under Codes/ as an earlier
  note in this file said. ID-join note: the phenotype ID is an int, e.g. 28046; the folder ID
  is zero-padded, e.g. 0028046 -- the fetcher matches on the integer value.)

FILES IN THIS FOLDER
  connectomes_full.csv  -- N*V rows by V columns (150*200 = 30,000 rows by 200 columns), no
                            header, one 200x200 Fisher-z connectivity matrix stacked per
                            subject, in the same order as subjects.csv. Subject i (0-indexed)
                            is rows [i*200 : (i+1)*200), all 200 columns.
  subjects.csv  -- one row per subject, in stacking order (see dictionary below).
  meta.csv  -- one header + one data row: N=150, V=200, atlas=cc200, pipeline=acpi,
               strategy=OutputsTimeSeries.

VARIABLE DICTIONARY FOR subjects.csv
  Columns (confirmed from the file): FILE_ID, SITE_ID, AGE_AT_SCAN, SEX, BMI.
  FILE_ID       string. Zero-padded ACPI subject ID, e.g. 0028046.
  SITE_ID       string. Confirmed constant value in this sample: "ACPI" for all 150 rows --
                the per-site source (mta_1/nyu_1/du_1) is not preserved in this column.
  AGE_AT_SCAN   numeric, years. Confirmed range in this sample: 20.3 to 50.73.
  SEX           integer, ABIDE-style convention carried over: confirmed values 1 and 2.
                Split in this sample: 27 coded 1, 123 coded 2 (123F/27M, i.e. 1 = male,
                2 = female).
  BMI           intended to be numeric body-mass index. CONFIRMED IN THIS SAMPLE: every one
                of the 150 rows holds the placeholder character "#", not a numeric value --
                BMI is effectively 100% missing in this downloaded file, even though it is
                a column header and is listed as an available outcome in the aggregated
                ACPI phenotype file upstream. Do not use this column as-is; if BMI is
                needed, re-check the source aggregated CSV directly and fix the fetcher's
                column mapping first.

COUNTS
  Subjects: 150. Rows in connectomes_full.csv: 30,000 (= 150 x 200). Nodes V = 200 (CC200).
  Edge dimension p = V(V-1)/2 = 19,900.

KNOWN CAVEATS
  BMI is present as a column but is unusable -- see above; treat this folder as carrying
  only age and sex as usable outcomes, not "age, sex, BMI" as an earlier note in this file stated.
  Cannabis/marijuana-use group is NOT in this file at all -- it lives only in the separate
  site phenotype files (mta_1_/nyu_1_/du_1_phenotypic_data.csv) linked above, which were
  never downloaded or joined here.
  SITE_ID collapses to a single constant "ACPI" -- if a real per-site variable is ever
  needed, it must be re-derived from the site file names, not from this column.

CITE: derived from the same Preprocessed Connectomes Project pipeline family as ABIDE.
