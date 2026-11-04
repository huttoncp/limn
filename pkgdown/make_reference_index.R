## Generate a grouped pkgdown reference index from the actual Rd topics, so
## that no topic can be left out of the index by hand.
## Run from the package root.

rd <- list.files("man", pattern = "[.]Rd$")
topic <- sub("[.]Rd$", "", rd)
## pkgdown indexes by the Rd \name, which is what the filename encodes except
## where roxygen escaped a character; read the name to be certain
nm <- vapply(file.path("man", rd), function(f) {
  x <- readLines(f, warn = FALSE)
  i <- grep("^\\\\name\\{", x)[1]
  if (is.na(i)) NA_character_ else sub("^\\\\name\\{(.*)\\}.*$", "\\1", x[i])
}, character(1))
topic <- ifelse(is.na(nm), topic, nm)
topic <- setdiff(topic, "limn-package")

## ---- groups, in the order they should appear -------------------------------
## Each entry: title, a one-line description, a predicate over topic names,
## and optionally the topics to list first, in that order -- the front door
## of a group before its parts. Whatever else matches follows alphabetically.
##
## Edit THIS file, not _pkgdown.yml: the yml is regenerated from it, and a
## hand edit there is lost the next time an export is added.
grp <- list(
  list("Predictions and forecasts",
       "Predictive distributions from a fitted illume model, the rows a forecast needs, and covariates whose future is not known.",
       function(x) x %in% c("ilm_predict", "ilm_forecast", "ilm_future", "ilm_covariate_paths",
                            "ilm_as_draws", "as.data.frame.ilm_pred"),
       c("ilm_predict", "ilm_forecast", "ilm_future", "ilm_covariate_paths", "ilm_as_draws",
         "as.data.frame.ilm_pred")),

  list("Totals and averages",
       "Sums and weighted averages of predictions over rows, taken draw by draw.",
       function(x) x %in% c("ilm_aggregate", "ilm_average")),

  list("Scores and calibration",
       "How well predictive distributions match the outcomes: proper scores, scaled errors, probability integral transforms and calibration.",
       function(x) grepl("^ilm_(scores|naive_scale|pit|calibration)", x),
       c("ilm_scores.ilm_pred", "ilm_naive_scale", "ilm_pit", "ilm_calibration.ilm_pred")),

  list("Cross-validation and comparison",
       "Folds that follow the data's structure, cross-validated scores, and model comparisons with standard errors.",
       function(x) x %in% c("ilm_folds", "ilm_cv", "ilm_compare"),
       c("ilm_folds", "ilm_cv", "ilm_compare")),

  list("Plots",
       "Fan charts, scores by horizon, model comparisons, PIT histograms and calibration.",
       function(x) grepl("^ilm_plot_", x),
       c("ilm_plot_fan", "ilm_plot_horizon", "ilm_plot_compare", "ilm_plot_pit",
         "ilm_plot_calibration")),

  list("Extending limn",
       "Methods a package defines so that limn can predict from, and refit, its own models.",
       function(x) x %in% c("ilm_predict_draws", "ilm_cv_refit"))
)

assigned <- character(0)
lines <- c("reference:")
for (g in grp) {
  hit <- Filter(function(x) isTRUE(tryCatch(g[[3]](x), error = function(e) FALSE)),
                setdiff(topic, assigned))
  hit <- sort(hit)
  if (length(g) >= 4L) hit <- c(intersect(g[[4]], hit), setdiff(hit, g[[4]]))
  if (!length(hit)) next
  assigned <- c(assigned, hit)
  lines <- c(lines,
             paste0("- title: \"", g[[1]], "\""),   # quoted: titles may hold a colon
             paste0("  desc: >"),
             paste0("    ", g[[2]]),
             "  contents:",
             paste0("  - ", hit))
}

left <- sort(setdiff(topic, assigned))
if (length(left)) {
  lines <- c(lines, "- title: \"Other topics\"",
             "  desc: >",
             "    Topics not yet placed in a group above.",
             "  contents:",
             paste0("  - ", left))
}

head <- c(
  "url: https://huttoncp.github.io/limn/",
  "template:",
  "  bootstrap: 5",
  "  bslib:",
  "    primary: \"#2a6f97\"",
  "",
  "home:",
  "  title: A Projection Engine for Distributional Forecasting",
  "",
  "navbar:",
  "  structure:",
  "    left:  [intro, reference, articles, news, illume, illumex]",
  "    right: [search, github]",
  "  components:",
  "    illume:",
  "      text: Modelling (illume)",
  "      href: https://huttoncp.github.io/illume/",
  "    illumex:",
  "      text: Exploration (illumex)",
  "      href: https://huttoncp.github.io/illumex/",
  "",
  "articles:",
  "- title: Articles",
  "  navbar: ~",
  "  contents:",
  "  - predictions",
  "  - forecasting",
  "")

writeLines(c(head, lines), "_pkgdown.yml")
cat("topics:", length(topic), " assigned:", length(assigned),
    " unassigned:", length(left), "\n")
if (length(left)) cat("UNASSIGNED:", paste(left, collapse = ", "), "\n")
