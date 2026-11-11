DATASET: OSF "Schaefer 100 parcel 17 networks" -- task-fMRI functional connectomes (V=100
replication of the Power-229 / Guha-Kundu benchmark dataset)
=========================================================================================

WHAT THIS FOLDER IS
A coarser-parcellation replication of the same 144-subject, 7-condition task-fMRI dataset
used throughout this project as data/osf_guha/ (Power atlas, V=229). Here the SAME subjects
and conditions are re-parcellated at V=100 (Schaefer 100-parcel, 17-network scheme), giving
edge dimension p = V(V-1)/2 = 4,950 -- small enough that raw-edge ARD and raw-edge spike-
and-slab SSVS (network_bayes_simulation.R::raw_edge_ard / spike_slab_ssvs), infeasible at
V=229's p=26,106, become tractable. This is the first real-data run of the full five-method
comparison from the simulations (graph-spectral ARD, hybrid, raw-edge ridge, raw-edge lasso,
raw-edge ARD, raw-edge spike-slab SSVS). Open OSF download, no data-use agreement.

EXACT SOURCE
  Project page: https://osf.io/hnj7s/  (title: "Schaefer 100 parcel 17 networks", a sibling
    OSF component of node 5zyws [Power-229] under parent project m5crs)
  File listing (API, paginated; each file's download URL is in data[i].links.download):
    https://api.osf.io/v2/nodes/hnj7s/files/osfstorage/
  File-name pattern: Schaefer100_17_sub-<ID>_<cond>.txt
  Age + network labels (same sibling OSF component dm23b, "Atlas and Participant info", as
    used for data/osf_guha/ -- shared across the Power-229 and both Schaefer atlases):
    https://osf.io/download/u3xpd/  ->  participant_ages.csv (not stored separately in this
      folder -- AGE is already joined into subjects.csv, see below)
    Schaefer_100parcel_17network_labels.csv (stored in this folder) -- node-to-network
      assignment for the 100 ROIs, one row per node.

CONFIRMED ENUMERATION QUIRK (same as Power-229): the OSF API listing for node hnj7s returns
1007 raw entries but only 1003 unique file IDs -- one filename is listed multiple times under
distinct file IDs (the same quirk documented for node 5zyws in data/osf_guha/README.txt).
The fetcher below de-duplicates by file ID, so the cache holds the correct 1003 unique files.

HOW TO (RE)DOWNLOAD / REBUILD (resumable -- caches to data/osf_schaefer100/cache/, skips
what is already there)
  python3 Codes/fetch/osf_guha_fetch.py --node hnj7s --prefix Schaefer100_17 \
    --outdir data/osf_schaefer100 --workers 6
  (Codes/fetch/osf_guha_fetch.py is the SAME generalized fetcher used for Power-229 and
  Schaefer-200 -- --node/--prefix/--outdir select the OSF component, filename prefix, and
  destination folder; the pagination + de-dup-by-file-ID logic is unchanged from the
  Power-229 fix.) After the cache is complete, rebuild the tabular files (subjects.csv,
  connectomes_full.csv, meta.csv, master.rds) with:
    Rscript -e 'source("Codes/fetch/build_osf_master_generic.R");
                build_osf_master_generic(cache_dir = "data/osf_schaefer100/cache",
                                          prefix = "Schaefer100_17",
                                          out_dir = "data/osf_schaefer100")'
  (build_osf_master_generic.R is the atlas-parameterized twin of
  Codes/fetch/build_osf_master.R -- a callable function, not a CLI script -- used because
  that script is hardcoded to the Power-229 folder/filename pattern. Ages are joined from
  data/osf_guha/participant_ages.csv, the same 144 subjects across all three parcellations.)

FILES IN THIS FOLDER
  cache/  -- 1003 cached per-file matrices (one per condition-file), the resumable cache
             osf_guha_fetch.py writes and reads back. Confirmed count: 1003, matching the
             1003 unique file IDs on the OSF API.
  connectomes_full.csv  -- N*V rows by V columns (1003*100 = 100,300 rows by 100 columns), no
                            header, one 100x100 matrix stacked per file, in the same row
                            order as subjects.csv. File i (0-indexed) occupies rows
                            [i*100 : (i+1)*100), all 100 columns.
  subjects.csv  -- one row per file, 9 columns (see dictionary below). Confirmed 1003 data
                    rows, 144 unique SUBJECT values.
  meta.csv  -- one header + one data row: N=1003, V=100, atlas=Schaefer100_17,
               pipeline=osf_master_generic, strategy=all_conditions.
  master.rds  -- R list with elements conn (list of 1003 100x100 numeric matrices), subjects
                 (the same data.frame as subjects.csv), V=100. This is the file every
                 analysis script in Codes/analysis/*_sch100.R and
                 Codes/analysis/multi_resolution_common.R reads via readRDS().
  Schaefer_100parcel_17network_labels.csv  -- 100 rows (+ header), columns ROI_number,
                            ROI_label, Network, Network_color: which of the 17 Schaefer
                            functional networks (e.g. VisCent, DorsAttnA, DefaultA) each of
                            the 100 nodes belongs to, plus a hex/name color for plotting.
                            Confirmed 100 data rows.

VARIABLE DICTIONARY FOR subjects.csv (identical schema to data/osf_master/'s Power-229
master -- same 144 subjects, same derivation logic, just re-parcellated)
  FILE_ID   string, the original OSF filename without extension, e.g.
            Schaefer100_17_sub-001_gng_inhibition.
  SUBJECT   integer, subject id (144 unique values, 1-144 as stored -- NOT zero-padded like
            the Power-229 folder's string SUBJECT; confirmed by direct inspection).
  COND      string, one of: rest, nbk_0back, nbk_1back, nbk_2back, gng_initiation,
            gng_inhibition, tsw_shifting (confirmed, 7 distinct values). Per-condition file
            counts (confirmed by direct count): nbk_1back 144, tsw_shifting 144,
            gng_inhibition 143, gng_initiation 143, nbk_0back 143, nbk_2back 143, rest 143.
  DOMAIN    string, one of: inhibition, initiation, shifting, wm, rest (5 levels; wm =
            working memory, covers all 3 n-back loads pooled). Confirmed counts: inhibition
            143, initiation 143, rest 143, shifting 144, wm 430 (= 143+144+143, the three
            n-back loads combined -- wm is 3x over-represented relative to the other domains
            if used unpooled; see CLAUDE.md's domain_4class_balanced / domain_5class_all
            outcomes and their documented balancing judgement call).
  LOAD      integer 0/1/2 for the three n-back conditions, NA for every other condition (573
            NAs, confirmed = 1003 - 430 wm rows). Only meaningful together with
            COND %in% c("nbk_0back","nbk_1back","nbk_2back").
  AGE       integer, years, range 20-86 (confirmed, 0 NAs). One value per SUBJECT, repeated
            across that subject's condition rows -- this is a participant-level, not a
            condition-level, variable.
  AGE_BIN   integer 0/1, median split (0 = younger, 1 = older). Confirmed counts: 0->500,
            1->503 rows (144 subjects, so not exactly balanced at the subject level once
            expanded to condition-rows because subjects with slightly more or fewer complete
            conditions contribute more or fewer rows to whichever side of the median they
            fall on).
  AGE_Q     integer 1-4, quartile split (1 = youngest quartile ... 4 = oldest). Confirmed
            counts: 272 / 228 / 252 / 251 rows.
  IS_TASK   integer 0/1, 0 = rest, 1 = any task condition. Confirmed counts: 0->143 (matches
            the rest row count exactly), 1->860.

COUNTS
  1003 files, 144 distinct subjects, V=100 (Schaefer 100-parcel, 17-network atlas). Edge
  dimension p = V(V-1)/2 = 4,950. Conditions-per-subject: 139 subjects have all 7 conditions,
  5 subjects have 6 (none have fewer) -- confirmed identical to the Power-229 folder's
  coverage, since both are re-parcellations of the exact same underlying scans.

KNOWN CAVEATS
  Age is the ONLY participant-level variable released anywhere on this OSF project (no sex,
  education, MMSE, or head-motion data) -- confirmed for this project's dm23b component
  listing at the time the Power-229 folder was audited (see data/osf_guha/README.txt); the
  same participant info component is shared across all three atlases, so this holds here too.
  The apparent "1007 files" figure is an OSF listing artifact identical to the one documented
  for Power-229 (data/osf_guha/README.txt) -- one filename listed under 5 distinct file IDs,
  not 45+ genuinely missing files; the fetcher's page-number enumeration + de-dup-by-file-ID
  logic (fixed for Power-229 on 2026-09-10) already handles this correctly here, confirmed by
  the enumeration self-check recorded in FETCH_STATE.txt (1007 raw entries, 1003 unique file
  IDs, 1003 unique filenames -- nothing missing).
  This is a re-parcellation, NOT new subjects or new conditions -- every number in this file
  that concerns subjects/conditions/coverage (144 subjects, 139 with all 7 conditions, age
  range 20-86, etc.) is identical to data/osf_guha/'s, by construction. Only V (and therefore
  p and the connectome values themselves) differs.
  Schaefer-200 (OSF component e4nph, V=200, p=19,900) is a further sibling of this same
  dataset. It has NOT been fetched as part of this work (explicitly out of scope for this
  session) -- V=200 is close enough to the existing V=229 that it unlocks no new capability
  the way V=100 does for raw-edge ARD/SSVS.

CITE: Rieck et al. 2021 (Neurobiology of Aging, the source data collection); Rodriguez-Acosta,
Guha & Dinov 2024 (BMGC, the paper this project benchmarks against on this data); Schaefer et
al. 2018 (the Schaefer parcellation itself).
