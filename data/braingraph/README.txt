DATASET: braingraph.org OASIS-3 -- STRUCTURAL (diffusion) connectomes, GraphML
=========================================================================================

WHAT THIS FOLDER IS
Structural (diffusion-tractography) brain connectomes for the OASIS-3 cohort, at three
parcellation resolutions, used as a second modality alongside the ABIDE/ACPI functional
connectomes -- structural networks are typically sparser and more low-rank, so this is the
candidate dataset where the graph-spectral method might win on raw accuracy, not just
parsimony. Connectomes are here; participant-level outcome labels are NOT (see caveats).

EXACT SOURCE
  Direct download (7z archives), no HCP-style click-through:
    scale1 (124 nodes, 372 MB): https://braingraph.org/static/oasis3_graphmls_scale1.7z
    scale2 (170 nodes, 624 MB): https://braingraph.org/static/oasis3_graphmls_scale2.7z
    scale3 (272 nodes, 1.2 GB): https://braingraph.org/static/oasis3_graphmls_scale3.7z
    scale4 (502 nodes, 2.5 GB): https://braingraph.org/static/oasis3_graphmls_scale4.7z
    scale5 (1058 nodes, 5.5 GB): https://braingraph.org/static/oasis3_graphmls_scale5.7z
    Download page (also lists HCP structural): https://braingraph.org/cms/download-pit-group-connectomes/
  Extract with p7zip, e.g.:
    7z x oasis3_graphmls_scale2.7z -odata/braingraph/oasis3_graphmls_scale2/

OUTCOME / LABELS (gated, free registration, not yet obtained)
  OASIS-3 age/CDR/diagnosis: register at https://central.xnat.org, accept the OASIS-3 data
  use agreement, then download "ADRC Clinical Data" (age at entry, CDR, dx1-dx5). Save as
  data/braingraph/oasis3_demographics.csv with columns subject_id, age, cdr. Join key is the
  sub-OAS3xxxx token already present in every GraphML filename in this folder.

HOW TO USE (label-free analyses work today)
  Rscript Codes/pipelines/braingraph_oasis_pipeline.R data/braingraph/oasis3_graphmls_scale2
  (Note: the script now lives under Codes/pipelines/, not directly under Codes/ as an
  earlier note in this file said.)

FILES IN THIS FOLDER
  oasis3_graphmls_scale1/, oasis3_graphmls_scale2/, oasis3_graphmls_scale3/  -- one directory
    per resolution, each holding one GraphML file per scan session. Confirmed: 975 files in
    each of scale1 and (by the same extraction) scale2/scale3. Filename pattern:
    sub-OAS3xxxx_ses-dNNNN_atlas-L2018_res-scaleK_conndata-network_connectivity.graphml
    (sub-OAS3xxxx = OASIS-3 subject id; ses-dNNNN = days since a study reference date, i.e.
    session/visit; scaleK = the resolution, K in {1,2,3} present here). No scale4/scale5.
  Each GraphML file is a single subject-session's structural connectivity graph: nodes = the
  atlas's regions (124/170/272 depending on scale), edges carry weights including
  number_of_fibers and normalized_fiber_density. There is no subject.csv or connectomes_full
  in this folder -- these are raw GraphML files, not yet assembled into the project's shared
  tabular format; Codes/pipelines/braingraph_oasis_pipeline.R reads them directly via igraph.
  There is no graph-level attribute carrying age, CDR, or diagnosis in these files -- confirmed
  by inspection when this pipeline was first written.

VARIABLE DICTIONARY
  Not applicable -- no subjects.csv exists for this folder. Any subject-level variable would
  have to come from the separate, gated oasis3_demographics.csv described above, joined on
  the sub-OAS3xxxx token.

COUNTS
  975 GraphML files in each of scale1/scale2/scale3 (one file per session; multiple sessions
  per subject are possible, e.g. sub-OAS30001 and sub-OAS30002 each appear with two distinct
  ses- tokens in scale1). 695 distinct subjects once sessions are collapsed to one-per-subject
  (not independently recomputed here, since doing so needs
  the same session-collapsing logic already implemented in braingraph_oasis_pipeline.R).
  Nodes V = 124 (scale1), 170 (scale2), 272 (scale3).

KNOWN CAVEATS
  Connectomes ARE here, but participant-level outcomes (age, CDR, diagnosis) are NOT -- they
  sit behind the OASIS-3 data use agreement at central.xnat.org. No supervised analysis is
  possible on this folder until oasis3_demographics.csv is obtained and joined.
  The unnormalized graph Laplacian basis is DEGENERATE on these structural connectomes (an
  earlier session found ~54/55 low-frequency features have essentially zero between-subject
  variance) -- use the adjacency or normalized-Laplacian basis instead, never the raw
  unnormalized Laplacian, on this dataset.
  Only scale1/2/3 were downloaded; scale4 (502 nodes) and scale5 (1058 nodes) were not.

CITE: Varga & Grolmusz 2020 (braingraph.org); OASIS-3 (LaMontagne et al. 2019).
