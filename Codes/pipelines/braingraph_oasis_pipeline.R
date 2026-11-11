## =====================================================================================
## braingraph.org OASIS-3 STRUCTURAL connectomes -> graph-spectral pipeline.
##

## distributes OASIS-3 structural connectomes with DIRECT download links,
## as GraphML inside 7z archives:
##
##   scale1: https://braingraph.org/static/oasis3_graphmls_scale1.7z  (372 MB, 696 subj, 124 nodes)
##   scale2: https://braingraph.org/static/oasis3_graphmls_scale2.7z  (624 MB, 696 subj, 170 nodes)
##   scale3: https://braingraph.org/static/oasis3_graphmls_scale3.7z  (1.2 GB, 696 subj, 272 nodes)
##   scale4: .../oasis3_graphmls_scale4.7z (2.5 GB, 502 nodes); scale5 (5.5 GB, 1058 nodes)
##
## ----------------------------------------------------------------------------------------------
## HOW TO USE:
##   1. Download one scale, e.g. scale2, into data/braingraph/
##   2. Extract it (needs p7zip:  brew install p7zip):
##        7z x data/braingraph/oasis3_graphmls_scale2.7z -odata/braingraph/oasis3_scale2/
##   3. OUTCOME: structural connectomes carry no bundled phenotype. OASIS-3 age / clinical
##      diagnosis (CDR) come from the OASIS-3 release (free application at central.xnat.org).
##      Save them as data/braingraph/oasis3_demographics.csv with columns: subject_id, age, cdr
##      (subject_id must match the GraphML file-name subject token). If graph attributes embed
##      age/sex, this script will pick them up automatically and the CSV is optional.
##   4. Run:  Rscript Codes/pipelines/braingraph_oasis_pipeline.R data/braingraph/oasis3_scale2
##
## Structural connectomes are sparser and more low-rank than functional ones, so this is the
## candidate dataset where graph-spectral shrinkage may win on ACCURACY (not just parsimony).
## =====================================================================================

suppressMessages({
  source("Codes/core/abide_pipeline.R")          # build_features, graph_spectral_ard, helpers
  stopifnot(requireNamespace("igraph", quietly = TRUE))
})

## Read one GraphML into a weighted adjacency matrix + any embedded subject attributes.
## braingraph OASIS-3 edge attributes (verified): number_of_fibers, fiber_density,
## normalized_fiber_density, fiber_length_*, shore_* (NODDI/SHORE). We default to
## normalized_fiber_density (size-normalized streamline density), a standard structural weight.
read_graphml_connectome <- function(path, weight_attr = NULL) {
  g <- igraph::read_graph(path, format = "graphml")
  eattrs <- igraph::edge_attr_names(g)
  if (is.null(weight_attr)) {
    pref <- c("normalized_fiber_density", "fiber_density", "number_of_fibers",
              "edge_weight", "weight")
    weight_attr <- pref[pref %in% eattrs][1]
    if (is.na(weight_attr) && length(eattrs)) weight_attr <- eattrs[1]
  }
  A <- as.matrix(igraph::as_adjacency_matrix(
    g, attr = if (!is.na(weight_attr)) weight_attr else NULL, sparse = FALSE))
  A <- (A + t(A)) / 2; diag(A) <- 0
  gattr <- igraph::graph_attr(g)            # may contain age/sex/dx if braingraph embedded them
  list(A = A, gattr = gattr, weight_attr = weight_attr)
}

## subject id from a braingraph filename: sub-OAS30001_ses-d0757_..._scaleN_... -> OAS30001
oasis_subject_id <- function(fname) sub(".*sub-(OAS3[0-9]+)_.*", "\\1", basename(fname))

load_braingraph <- function(dir, demographics_csv = file.path(dirname(dir), "oasis3_demographics.csv"),
                            weight_attr = NULL, one_session_per_subject = TRUE) {
  files <- list.files(dir, pattern = "[.]graphml$", full.names = TRUE)
  if (!length(files)) stop("No .graphml files found in ", dir,
                           " -- extract the 7z first (see header).")
  if (one_session_per_subject) {                 # keep earliest session per subject
    sid_all <- oasis_subject_id(files)
    files <- files[!duplicated(sid_all)]
  }
  cat(sprintf("Reading %d GraphML connectomes from %s ...\n", length(files), dir))
  conn <- list(); ids <- character(); embedded <- list()
  V <- NULL; wattr <- NULL
  for (f in files) {
    r <- read_graphml_connectome(f, weight_attr = weight_attr)
    if (is.null(V)) { V <- nrow(r$A); wattr <- r$weight_attr }
    if (nrow(r$A) != V) next                # keep a single resolution
    conn[[length(conn) + 1]] <- r$A
    ids <- c(ids, oasis_subject_id(f))
    embedded[[length(embedded) + 1]] <- r$gattr
  }
  cat(sprintf("  loaded %d connectomes at V=%d (edge weight: %s)\n", length(conn), V, wattr))

  ## outcome: prefer embedded graph attribute 'age'/'sex'; else join external demographics
  subj <- data.frame(subject_id = ids, stringsAsFactors = FALSE)
  emb_age <- sapply(embedded, function(a) if (!is.null(a$age)) as.numeric(a$age) else NA)
  if (any(is.finite(emb_age))) {
    subj$AGE_AT_SCAN <- emb_age
    cat("  using age embedded in GraphML graph attributes\n")
  } else if (file.exists(demographics_csv)) {
    dem <- read.csv(demographics_csv)
    key <- intersect(c("subject_id", "ID", "MR_ID"), names(dem))[1]
    subj <- merge(subj, dem, by.x = "subject_id", by.y = key, all.x = TRUE, sort = FALSE)
    if (!"AGE_AT_SCAN" %in% names(subj) && "age" %in% names(subj))
      subj$AGE_AT_SCAN <- subj$age
    cat("  joined external demographics from ", demographics_csv, "\n")
  } else {
    cat("  WARNING: no embedded age and no demographics CSV at", demographics_csv,
        "\n           -> connectomes load, but the regression needs an outcome. See header.\n")
  }
  list(conn = conn, subjects = subj, V = V, N = length(conn))
}

run_braingraph <- function(dir, outcome = "AGE_AT_SCAN", Ks = c(10, 20, 30),
                           bases = c("laplacian", "adjacency"), n_splits = 10) {
  d <- load_braingraph(dir)
  if (!outcome %in% names(d$subjects) || all(!is.finite(d$subjects[[outcome]]))) {
    cat("\nNo usable outcome -- stopping after connectome load. Provide demographics (see header).\n")
    return(invisible(d))
  }
  ## reuse the same repeated-split machinery as ABIDE
  y_all <- d$subjects[[outcome]]; keep <- which(is.finite(y_all))
  conn <- d$conn[keep]; y_all <- y_all[keep]; N <- length(conn)
  cat(sprintf("braingraph OASIS-3: N=%d, V=%d, outcome=%s\n", N, d$V, outcome))
  res <- list()
  for (s in seq_len(n_splits)) {
    set.seed(300 + s)
    te <- sample(seq_len(N), round(0.30 * N)); tr <- setdiff(seq_len(N), te)
    for (basis in bases) for (K in Ks) {
      ft <- build_features(conn, tr, K, basis = basis)
      gs <- graph_spectral_ard(ft$X_spec[tr, ], y_all[tr], ft$X_spec[te, ])
      rg <- raw_edge_ridge_woodbury(ft$X_edge[tr, ], y_all[tr], ft$X_edge[te, ])
      res[[length(res) + 1]] <- data.frame(split = s, basis = basis, K = K, dim = ft$q,
        gs_rmse = rmse(y_all[te], gs$pred), ridge_rmse = rmse(y_all[te], rg$pred))
    }
  }
  df <- do.call(rbind, res)
  dir.create("results", showWarnings = FALSE)
  write.csv(df, "results/braingraph_oasis_results.csv", row.names = FALSE)
  cat("\nMean test RMSE by (basis,K):\n")
  print(aggregate(cbind(gs_rmse, ridge_rmse) ~ basis + K, df, mean))
  cat("Saved results/braingraph_oasis_results.csv\n")
  invisible(df)
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) run_braingraph(args[1]) else
    cat("Usage: Rscript Codes/pipelines/braingraph_oasis_pipeline.R <dir-with-extracted-graphmls>\n")
}
