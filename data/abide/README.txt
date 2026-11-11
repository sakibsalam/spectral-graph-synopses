DATASET: ABIDE Preprocessed (PCP) -- resting-state functional connectomes (CC200 atlas)
=========================================================================================

WHAT THIS FOLDER IS
200 resting-state functional-connectivity matrices from the ABIDE (Autism Brain Imaging Data
Exchange) Preprocessed initiative, CC200 parcellation, used as the primary real-data example
in this project (age regression and autism-vs-control classification). Fully open data, no
data-use agreement required.

EXACT SOURCE
  Phenotype CSV:
    https://s3.amazonaws.com/fcp-indi/data/Projects/ABIDE_Initiative/Phenotypic_V1_0b_preprocessed1.csv
  ROI time series (per subject), pattern:
    https://s3.amazonaws.com/fcp-indi/data/Projects/ABIDE_Initiative/Outputs/cpac/filt_global/rois_cc200/<FILE_ID>_rois_cc200.1D
    (pipeline=cpac, strategy=filt_global, atlas=rois_cc200; <FILE_ID> e.g. NYU_0051091)
  Project page: http://preprocessed-connectomes-project.org/abide/

HOW TO (RE)DOWNLOAD / REBUILD
  python3 Codes/fetch/abide_fetch.py --n 200 --atlas cc200 --outdir data/abide
  (Note: the script now lives under Codes/fetch/, not directly under Codes/ as an earlier
  note in this file said. Pulls the phenotype CSV, downloads each subject's ROI time series,
  builds a Fisher-z correlation connectome per subject, and writes connectomes_full.csv,
  subjects.csv, and meta.csv below.)

FILES IN THIS FOLDER
  Phenotypic_V1_0b_preprocessed1.csv  -- the raw ABIDE phenotype table as downloaded (all
                                         ABIDE sites/subjects, most of which are not used here).
  connectomes_full.csv  -- the built connectome data. No header row. N*V rows by V columns
                            (200*200 = 40,000 rows by 200 columns), i.e. one 200x200 subject
                            connectivity matrix stacked directly on top of the next, in the
                            same order as the rows of subjects.csv. To slice out subject i
                            (0-indexed): rows [i*200 : (i+1)*200), all 200 columns.
  subjects.csv  -- one row per subject actually used, in the order the matrices are stacked
                    in connectomes_full.csv (see the dictionary below).
  meta.csv  -- one header row + one data row: N, V, atlas, pipeline, strategy. Confirmed
               values here: N=200, V=200, atlas=cc200, pipeline=cpac, strategy=filt_global.

VARIABLE DICTIONARY FOR subjects.csv
  Columns (confirmed from the file): FILE_ID, SITE_ID, DX_GROUP, AGE_AT_SCAN, SEX, FIQ, VIQ, PIQ.
  FILE_ID       string. ABIDE subject identifier, e.g. NYU_0050952. Also the ROI-file stem.
  SITE_ID       string. Confirmed values in this 200-subject sample: NYU, UM_1 only (of the
                many ABIDE sites -- this download used just these two).
  DX_GROUP      integer, ABIDE convention. Confirmed values: 1 = autism, 2 = control.
                Split in this sample: 100 / 100.
  AGE_AT_SCAN   numeric, years. Confirmed range in this sample: 6.47 to 39.1.
  SEX           integer, ABIDE convention. Confirmed values: 1 = male, 2 = female.
                Split in this sample: 160 male / 40 female.
  FIQ, VIQ, PIQ  numeric IQ scores (full-scale, verbal, performance). ABIDE's missing-value
                sentinel is -9999. Confirmed in this sample: FIQ has 1 missing (-9999), VIQ
                has 0 missing, PIQ has 1 missing; valid-value ranges are FIQ 76-148, VIQ
                73-145, PIQ 64-149. Treat -9999 as NA before using these columns, it is not
                a real score.

COUNTS
  Subjects: 200. Rows in connectomes_full.csv: 40,000 (= 200 subjects x 200 nodes). Nodes
  V = 200 (CC200 atlas). Edge dimension p = V(V-1)/2 = 19,900.

KNOWN CAVEATS
  Only the outcomes actually present in the phenotype file are usable: diagnosis (DX_GROUP),
  age, sex, FIQ/VIQ/PIQ. There is no sex-balanced or site-balanced guarantee -- SEX is
  4:1 male:female and there are only two sites.
  -9999 in FIQ/VIQ/PIQ is a missing-value code, not a real score -- one subject each for
  FIQ and PIQ.
  This is the CC200 atlas variant (p=19,900 edges); raw-edge ARD and spike-and-slab are not
  computationally feasible at this resolution in this project's code -- that is the AAL
  variant's reason for existing (see data/abide_aal/README.txt).
  Real-data findings from this folder (age is a distributed signal favoring dense raw-edge
  ridge; autism-vs-control classification is where the graph-spectral method's parsimony
  shows) are recorded in CLAUDE.md section on key findings and in the 2026-06-22 workbook
  entries -- not repeated here.

CITE: Craddock et al. 2013 (Preprocessed Connectomes Project); Di Martino et al. 2014 (ABIDE).
