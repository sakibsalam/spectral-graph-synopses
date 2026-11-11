## =====================================================================================
## Parameterized version of Codes/fetch/build_osf_master.R, generalized to any of the three
## OSF parcellations of the SAME 144-subject / 7-condition dataset (Power-229, Schaefer-100,
## Schaefer-200), rather than being hardcoded to the Power-229 filename pattern. Does NOT
## modify build_osf_master.R itself -- that file, and the master it produces
## (data/osf_master/master.rds), are untouched. This is a NEW file used only to build
## data/osf_schaefer100/master.rds and data/osf_schaefer200/master.rds.
##
## Same schema as the existing master: FILE_ID, SUBJECT, COND, DOMAIN, LOAD, AGE, AGE_BIN,
## AGE_Q, IS_TASK -- derivation rules copied unchanged from build_osf_master.R (domain_of(),
## load_of(), AGE_BIN = median split, AGE_Q = quartile via quantile(0:4/4), IS_TASK = COND !=
## "rest"). Ages are joined from data/osf_guha/participant_ages.csv (the SAME 144 subjects
## across all three parcellations -- there is no separate ages file for Schaefer-100/200).
##
## Usage:
##   source("Codes/fetch/build_osf_master_generic.R")
##   build_osf_master_generic(cache_dir = "data/osf_schaefer100/cache", prefix = "Schaefer100_17",
##                            out_dir = "data/osf_schaefer100")
##   build_osf_master_generic(cache_dir = "data/osf_schaefer200/cache", prefix = "Schaefer200_17",
##                            out_dir = "data/osf_schaefer200")
## =====================================================================================

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

## Minimal .npy reader (float32/float64, C-order, no fortran flag) -- same as build_osf_master.R,
## duplicated here so this file has no dependency on that one (both stay independently runnable).
npy_load_generic <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  magic <- readBin(con, "raw", 6)
  ver <- readBin(con, "integer", 2, size = 1, signed = FALSE)
  hlen <- if (ver[1] == 1) readBin(con, "integer", 1, size = 2, signed = FALSE, endian = "little")
          else readBin(con, "integer", 1, size = 4, endian = "little")
  header <- rawToChar(readBin(con, "raw", hlen))
  dtype <- sub(".*'descr':\\s*'([^']+)'.*", "\\1", header)
  shape_str <- sub(".*'shape':\\s*\\(([^)]*)\\).*", "\\1", header)
  shape <- as.integer(strsplit(gsub("[ ,]+$", "", shape_str), ",\\s*")[[1]])
  n <- prod(shape)
  size <- if (grepl("f4", dtype)) 4 else 8
  vals <- readBin(con, "double", n, size = size, endian = "little")
  matrix(vals, nrow = shape[1], ncol = shape[2], byrow = TRUE)
}

build_osf_master_generic <- function(cache_dir, prefix, out_dir,
                                     ages_csv = "data/osf_guha/participant_ages.csv") {
  files <- list.files(cache_dir, pattern = "\\.npy$", full.names = FALSE)
  ## cache filenames are "<prefix>_sub-<id>_<cond>.txt.npy" -- strip BOTH suffixes.
  rx <- paste0("^", prefix, "_sub-([0-9]+)_(.+)\\.txt\\.npy$")
  pat <- regmatches(files, regexec(rx, files))
  ok <- lengths(pat) == 3
  if (!all(ok)) cat(sprintf("  skipping %d cache files that don't match the expected name pattern (prefix=%s)\n",
                            sum(!ok), prefix))
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
  if (nrow(s) < n_before) cat(sprintf("  dropped %d cached files whose SUBJECT has no entry in %s\n",
                                      n_before - nrow(s), ages_csv))

  s$DOMAIN <- domain_of(s$COND)
  s$LOAD   <- load_of(s$COND)
  med <- median(ages$Age)
  s$AGE_BIN <- as.integer(s$AGE > med)
  qs <- quantile(ages$Age, 0:4 / 4)
  s$AGE_Q <- as.integer(cut(s$AGE, qs, include.lowest = TRUE, labels = 1:4))
  s$IS_TASK <- as.integer(s$COND != "rest")

  s <- s[, c("FILE_ID", "SUBJECT", "COND", "DOMAIN", "LOAD", "AGE", "AGE_BIN", "AGE_Q", "IS_TASK")]
  s <- s[order(s$SUBJECT, s$COND), ]
  rownames(s) <- NULL

  conn <- lapply(s$FILE_ID, function(fid) npy_load_generic(file.path(cache_dir, paste0(fid, ".txt.npy"))))
  V <- ncol(conn[[1]])

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  master <- list(conn = conn, subjects = s, V = V)
  saveRDS(master, file.path(out_dir, "master.rds"))
  write.csv(s, file.path(out_dir, "subjects.csv"), row.names = FALSE)
  write.csv(data.frame(N = nrow(s), V = V, atlas = prefix, pipeline = "osf_master_generic",
                        strategy = "all_conditions"),
            file.path(out_dir, "meta.csv"), row.names = FALSE)

  n_subj <- length(unique(s$SUBJECT))
  cond_counts <- table(table(s$SUBJECT))   # distribution of #conditions-per-subject
  nback_ok <- length(unique(s$SUBJECT[s$COND %in% c("nbk_0back", "nbk_1back", "nbk_2back")]))
  cat(sprintf("Wrote %s: %d connectomes (V=%d), %d subjects\n",
              file.path(out_dir, "master.rds"), nrow(s), V, n_subj))
  cat("  conditions-per-subject distribution (subjects with N conditions):\n")
  print(cond_counts)
  invisible(master)
}
