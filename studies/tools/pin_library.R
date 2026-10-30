## pin_library: copy every package a study loads into a library of its own, so
## that a reinstall into the shared libraries cannot reach a running study
## (a rule set on 26 September). Nothing is
## downloaded or built: installed packages are copied as they are.
##
## Rscript studies/tools/pin_library.R <dest> <label> [extra packages...]
##   dest   the new library, outside the repository, such as
##          "<somewhere>/study-libs/illume-9c29c0e"
##   label  what the library is pinned to, recorded in its manifest, such as
##          "illume main 9c29c0e"
## The packages are limn's own dependencies (Depends, Imports, Suggests),
## pkgload, the study packages below and any extra ones named, the installed
## Suggests of all of those, and the dependencies of everything. Base packages
## come with R and are not copied. Check a new library with
## studies/tools/check_pinned.R before a long run.
args <- commandArgs(TRUE)
if (length(args) < 2L) stop("usage: Rscript studies/tools/pin_library.R <dest> <label> [packages...]")
dest <- args[1]; label <- args[2]
if (dir.exists(dest) && length(list.files(dest)))
  stop("'", dest, "' already exists and is not empty; a pinned library is made once.")
libs <- .libPaths()
ip <- utils::installed.packages(lib.loc = libs)
ip <- ip[!duplicated(ip[, "Package"]), , drop = FALSE]       # the first on the path wins
desc <- read.dcf("DESCRIPTION")
field <- function(f) if (f %in% colnames(desc)) desc[, f] else ""
deps <- unlist(lapply(c("Depends", "Imports", "Suggests"), function(f)
  trimws(sub("[(].*", "", strsplit(field(f), ",")[[1]]))))
want <- unique(c(setdiff(deps, c("R", "")), "pkgload", "tsibble", "fable", "fabletools",
                 "distributional", "dplyr", "sandwich", "scoringRules", args[-(1:2)]))
missing <- setdiff(want, rownames(ip))
if (length(missing)) stop("not installed: ", paste(missing, collapse = ", "))
## packages a study's packages load only when a feature runs, such as fable's
## ARIMA(), which needs feasts and urca: their Suggests, where installed. The
## first pinned library left these out, and every ARIMA fit failed.
sugg <- unlist(tools::package_dependencies(want, db = ip, which = "Suggests"))
sugg <- intersect(setdiff(sugg, want), rownames(ip))
want <- unique(c(want, sugg))
all_pkgs <- unique(c(want, unlist(tools::package_dependencies(
  want, db = ip, recursive = TRUE, which = c("Depends", "Imports", "LinkingTo")))))
base <- rownames(ip)[ip[, "Priority"] %in% "base"]
all_pkgs <- setdiff(intersect(all_pkgs, rownames(ip)), base)
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
for (p in all_pkgs) {
  ok <- file.copy(file.path(ip[p, "LibPath"], p), dest, recursive = TRUE, copy.date = TRUE)
  if (!ok) stop("could not copy ", p)
}
man <- data.frame(package = all_pkgs, version = ip[all_pkgs, "Version"],
                  from = ip[all_pkgs, "LibPath"],
                  built = vapply(all_pkgs, function(p) {
                    d <- read.dcf(file.path(dest, p, "DESCRIPTION"))
                    if ("Built" %in% colnames(d)) d[, "Built"] else NA_character_ }, ""),
                  row.names = NULL)
utils::write.csv(man, file.path(dest, "MANIFEST.csv"), row.names = FALSE)
writeLines(c(paste("label:", label), paste("made:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
             paste("R:", R.version.string), paste("packages:", nrow(man))),
           file.path(dest, "PINNED.txt"))
cat("pinned", nrow(man), "packages to", dest, "\n")
