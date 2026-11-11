## =====================================================================================
## Assemble data/osf_master/{master.rds, subjects.csv, meta.csv} from the osf_guha cache.
##
## This is the master-assembly step referenced throughout CLAUDE.md section 9 (the outcome
## explorer, Rieck-replication, and mode-mapping scripts all read data/osf_master/master.rds).
## No saved script for this step could be found under Codes/ when this file was added
## (2026-09-11) -- it was evidently run ad hoc in an earlier session. This version
## reconstructs it from the schema of the existing master.rds/subjects.csv (verified by
## inspecting data/osf_master/subjects.csv and reading master.rds back in R) so re-running it
## after a cache update (e.g. Codes/fetch/osf_guha_fetch.py picking up newly-available files)
## reproduces the same columns:
##   FILE_ID, SUBJECT, COND, DOMAIN, LOAD, AGE, AGE_BIN, AGE_Q, IS_TASK
##
## Verified derivation rules (checked against the existing master.rds before this file existed):
##   DOMAIN:  gng_inhibition -> inhibition; gng_initiation -> initiation;
##            nbk_0back/1back/2back -> wm; tsw_shifting -> shifting; rest -> rest
##   LOAD:    nbk_0back -> 0, nbk_1back -> 1, nbk_2back -> 2, else NA
##   IS_TASK: 0 for rest, 1 for every task condition
##   AGE:     joined from data/osf_guha/participant_ages.csv by SUBJECT (sub-NNN -> integer NNN)
##   AGE_BIN: 1 if AGE > median(AGE) ("older"), else 0 ("younger") -- median taken over the
##            144 subjects in participant_ages.csv (verified: median 47.5, AGE_BIN=1 range 48-86,
##            AGE_BIN=0 range 20-47)
##   AGE_Q:   quartile of AGE (1 = youngest .. 4 = oldest), via quantile(AGE, 0:4/4)
##
## Only subjects/conditions present in the cache are included -- rows are NOT invented for
## still-missing files. SUBJECT is stored as an integer (matches existing master.rds).
## =====================================================================================

CACHE   <- "data/osf_guha/cache"
AGES_CSV <- "data/osf_guha/participant_ages.csv"
OUT_DIR <- "data/osf_master"

domain_of <- function(cond) {
  ifelse(cond == "gng_inhibition", "inhibition",
  ifelse(cond == "gng_initiation", "initiation",
  ifelse(cond %in% c("nbk_0back", "nbk_1back", "nbk_2back"), "wm",
  ifelse(cond == "tsw_shifting", "shifting",
  ifelse(cond == "rest", "rest", NA_character_)))))
}
load_of <- function(cond) {
  ifelse(cond == "nbk_0back", 0L, ifelse(cond == "nbk_1back", 1L, ifelse(cond == "nbk_2back", 2L, NA_integer_)))
}

build_osf_master <- function(cache_dir = CACHE, ages_csv = AGES_CSV, out_dir = OUT_DIR) {
  files <- list.files(cache_dir, pattern = "\\.npy$", full.names = FALSE)
  ## cache filenames are "Power229_10_sub-<id>_<cond>.txt.npy" (the ".txt" is the original OSF
  ## filename extension, preserved by osf_guha_fetch.py before it appends ".npy" for np.save) --
  ## strip BOTH suffixes to recover the bare condition name (verified against a sample file: cond
  ## group must be "gng_inhibition", not "gng_inhibition.txt").
  pat <- regmatches(files, regexec("^Power229_10_sub-([0-9]+)_(.+)\\.txt\\.npy$", files))
  ok <- lengths(pat) == 3
  if (!all(ok)) cat(sprintf("  skipping %d cache files that don't match the expected name pattern\n", sum(!ok)))
  files <- files[ok]; pat <- pat[ok]
  subj_str <- vapply(pat, `[`, character(1), 2)
  cond     <- vapply(pat, `[`, character(1), 3)
  subj_int <- as.integer(subj_str)

  ages <- read.csv(ages_csv, stringsAsFactors = FALSE)
  ages$SUBJECT <- as.integer(sub("^sub-", "", ages$subjectID))

  s <- data.frame(
    FILE_ID = sub("\\.txt\\.npy$", "", files),
    SUBJECT = subj_int,
    COND    = cond,
    stringsAsFactors = FALSE
  )
  s <- merge(s, ages[, c("SUBJECT", "Age")], by = "SUBJECT", all.x = TRUE)
  names(s)[names(s) == "Age"] <- "AGE"

  n_before <- nrow(s)
  s <- s[!is.na(s$AGE), ]
  if (nrow(s) < n_before) cat(sprintf("  dropped %d cached files whose SUBJECT has no entry in %s\n", n_before - nrow(s), ages_csv))

  s$DOMAIN <- domain_of(s$COND)
  s$LOAD   <- load_of(s$COND)
  med <- median(ages$Age)
  s$AGE_BIN <- as.integer(s$AGE > med)
  qs <- quantile(ages$Age, 0:4 / 4)
  s$AGE_Q <- as.integer(cut(s$AGE, qs, include.lowest = TRUE, labels = 1:4))
  s$IS_TASK <- as.integer(s$COND != "rest")

  ## order matches the pre-existing master.rds/subjects.csv column order
  s <- s[, c("FILE_ID", "SUBJECT", "COND", "DOMAIN", "LOAD", "AGE", "AGE_BIN", "AGE_Q", "IS_TASK")]
  s <- s[order(s$SUBJECT, s$COND), ]
  rownames(s) <- NULL

  conn <- lapply(s$FILE_ID, function(fid) {
    m <- readBin  # placeholder to avoid lint; not used
    reticulate_not_used <- NULL
    m <- npy_load(file.path(cache_dir, paste0(fid, ".txt.npy")))
    m
  })
  V <- ncol(conn[[1]])

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  master <- list(conn = conn, subjects = s, V = V)
  saveRDS(master, file.path(out_dir, "master.rds"))
  write.csv(s, file.path(out_dir, "subjects.csv"), row.names = FALSE)
  write.csv(data.frame(N = nrow(s), V = V, atlas = "Power229", pipeline = "osf_master",
                        strategy = "all_conditions"),
            file.path(out_dir, "meta.csv"), row.names = FALSE)

  cat(sprintf("Wrote %s: %d connectomes (V=%d), %d subjects\n",
              file.path(out_dir, "master.rds"), nrow(s), V, length(unique(s$SUBJECT))))
  invisible(master)
}

## Minimal .npy reader (float32/float64, C-order, no fortran flag) -- avoids requiring
## reticulate/RcppCNPy just to read back what osf_guha_fetch.py wrote with np.save().
npy_load <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  magic <- readBin(con, "raw", 6)
  ver <- readBin(con, "integer", 2, size = 1, signed = FALSE)
  hlen <- if (ver[1] == 1) readBin(con, "integer", 1, size = 2, signed = FALSE, endian = "little")
          else readBin(con, "integer", 1, size = 4, endian = "little")
  header <- rawToChar(readBin(con, "raw", hlen))
  dtype <- regmatches(header, regexpr("'descr':\\s*'([^']+)'", header))
  dtype <- sub(".*'descr':\\s*'([^']+)'.*", "\\1", header)
  shape_str <- sub(".*'shape':\\s*\\(([^)]*)\\).*", "\\1", header)
  shape <- as.integer(strsplit(gsub("[ ,]+$", "", shape_str), ",\\s*")[[1]])
  n <- prod(shape)
  size <- if (grepl("f4", dtype)) 4 else 8
  vals <- readBin(con, "double", n, size = size, endian = "little")
  matrix(vals, nrow = shape[1], ncol = shape[2], byrow = TRUE)
}

if (sys.nframe() == 0) {
  build_osf_master()
}
