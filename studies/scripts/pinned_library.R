## A long study runs on a pinned library (studies/tools/pin_library.R), named by
## the environment variable LIMN_STUDY_LIB. Source this at the top of the
## script, before any package is loaded; call lmn_study_workers(cl) on a
## cluster before its workers load anything, and lmn_study_header() at the end
## for the lines that go at the top of the findings.

lmn_study_lib <- Sys.getenv("LIMN_STUDY_LIB")
if (!nzchar(lmn_study_lib) || !file.exists(file.path(lmn_study_lib, "PINNED.txt")))
  stop("set LIMN_STUDY_LIB to a pinned library (studies/tools/pin_library.R) before running ",
       "this study.", call. = FALSE)
.libPaths(c(lmn_study_lib, .Library))

## Every package's version and the modification time of its installed files,
## to compare before and after the run.
lmn_lib_snapshot <- function(lib = lmn_study_lib) {
  pk <- list.dirs(lib, recursive = FALSE, full.names = FALSE)
  data.frame(package = pk, version = vapply(pk, function(p)
    tryCatch(read.dcf(file.path(lib, p, "DESCRIPTION"), "Version")[1, 1], error = function(e) NA_character_), ""),
    mtime = vapply(pk, function(p) {
      f <- list.files(file.path(lib, p), recursive = TRUE, full.names = TRUE)
      if (length(f)) format(max(file.mtime(f)), "%Y-%m-%d %H:%M:%S") else NA_character_ }, ""),
    row.names = NULL)
}
lmn_lib_before <- lmn_lib_snapshot()

lmn_study_workers <- function(cl) {
  lib <- lmn_study_lib
  invisible(parallel::clusterCall(cl, function(l) { .libPaths(c(l, .Library)); invisible(NULL) }, lib))
}

## The header lines: what the library is pinned to, and whether any package in
## it changed during the run (which the pinned library should make impossible).
lmn_study_header <- function() {
  after <- lmn_lib_snapshot()
  m <- merge(lmn_lib_before, after, by = "package", all = TRUE, suffixes = c(".before", ".after"))
  changed <- m$package[is.na(m$version.before) | is.na(m$version.after) |
                         m$version.before != m$version.after | m$mtime.before != m$mtime.after]
  pinned <- readLines(file.path(lmn_study_lib, "PINNED.txt"))
  key <- c("illume", "illumex", "fable", "fabletools", "tsibble")
  v <- after[after$package %in% key, ]
  c(sprintf("Pinned library: %s (%s; %d packages). %s", sub("^label: ", "", pinned[1]),
            sub("^made: ", "made ", pinned[2]), nrow(after),
            paste(sprintf("%s %s", v$package, v$version), collapse = ", ")),
    if (length(changed)) sprintf("**Changed during the run:** %s. These results are not from one build.",
                                 paste(changed, collapse = ", "))
    else "No package in the library changed during the run.")
}
