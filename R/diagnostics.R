## ---- diagnostics: one source of text ------------------------------------------------
##
## Every check limn makes is one entry here, with illume's check fields plus a
## reference:
##   check      its name
##   status     OK, NOTE, WARN or FAIL
##   detail     what was found
##   cause      what it means, and why the remedy helps
##   suggestion the remedy, and how to apply it
##   reference  where the reasoning comes from, when a source exists
## Warnings, errors, printed verdicts and the help pages all compose their text
## from these entries, so they say the same thing. Braces are placeholders,
## filled from the situation; `doc` gives the words the help pages use for them.
## With the same fields as illume's checks, the entries can move unchanged to a
## mechanism shared across the family, should one be made.

lmn_ref <- list(
  kenward_roger = "Kenward and Roger (1997), Biometrics 53, 983-997",
  maas_hox = "Maas and Hox (2005), Methodology 1, 86-92",
  fpp3 = "Hyndman and Athanasopoulos (2021), Forecasting: Principles and Practice, 3rd edition, section 10.4",
  fpp3_cv = "Hyndman and Athanasopoulos (2021), Forecasting: Principles and Practice, 3rd edition, section 5.10",
  gneiting = "Gneiting, Balabdaoui and Raftery (2007), Journal of the Royal Statistical Society B 69, 243-268",
  van_calster = "Van Calster and others (2019), BMC Medicine 17, 230",
  coroneo_iacone = "Coroneo and Iacone (2020), Journal of Applied Econometrics 35(4), 391-409",
  diebold_mariano = "Diebold and Mariano (1995), Journal of Business and Economic Statistics 13, 253-263")

## A reference with the pages (or section) that support the entry citing it.
## A reference supports an entry's finding and its cause; a remedy cites it only
## where its text says so ("the reference, p. 4"). Anything in an entry that no
## source gives is limn's own, and says so (lmn_txt$ours), or is limn's own study.
lmn_cite <- function(ref, at) paste0(ref, ", ", at)

lmn_txt <- list(
  ours = "(limn's own reasoning)",
  few = "{term} has {k} levels, fewer than about 20.",
  narrow = paste("Too narrow: {coverage} of outcomes fell inside the central {level}% of their",
                 "predictive distributions, against {level}% expected."),
  narrow_cause = "The predictions are more certain than the data bear out.",
  few_new_cause = paste("With so few groups the variance between them is poorly determined: the",
                        "reference found its standard error too small with 50 groups or fewer, for",
                        "linear models fitted by restricted maximum likelihood. In limn's coverage",
                        "study with 8 groups, intervals for a new unit's observation covered 0.87 to",
                        "0.90 at 90%."),
  remedies = paste("illume::ilm_remedies() lists the model changes the fit's own checks support,",
                   "and illume::ilm_apply_remedy() refits the model with the one you choose."))

lmn_diagnostics <- list(
  ## ---- predictions (ilm_predict(), ilm_forecast()) ----
  few_groups_means = list(
    check = "few groups, means", status = "WARN", detail = lmn_txt$few,
    cause = paste("With so few groups the fit's normal approximation understates the uncertainty",
                  "in the fixed effects, since it takes no account of the uncertainty in the",
                  "variances (the reference shows this for linear mixed models), so intervals for",
                  "means cover as z-intervals do: in limn's coverage study with 8 groups, a 90%",
                  "interval for a typical unit's mean covered 0.84 to 0.86 for gaussian, Poisson,",
                  "negative binomial and beta outcomes (binary outcomes covered 0.92)."),
    suggestion = paste("No remedy for these intervals has been measured yet, so read them as somewhat",
                       "too narrow; more groups is the sure cure."),
    reference = lmn_cite(lmn_ref$kenward_roger, "pp. 983-984"),
    doc = list(term = "a random term", k = "k")),
  few_groups_new = list(
    check = "few groups, new units", status = "WARN", detail = lmn_txt$few,
    cause = lmn_txt$few_new_cause,
    ## Browne (1998), whom the reference cites for the REML statement below:
    ## awaiting the primary source; on the literature list.
    suggestion = paste("Refit with reml = TRUE in illume::ilm_model(): restricted maximum likelihood",
                       "is almost always at least as good as maximum likelihood, and sometimes better,",
                       "especially at estimating variances (the reference, p. 89, citing Browne",
                       "1998). In limn's prototype",
                       "study of few groups (studies/prototypes/) it raised 8-group coverage from 0.872",
                       "to 0.892 (Poisson) and from 0.891 to 0.903 (gaussian)."),
    reference = lmn_cite(lmn_ref$maas_hox, "pp. 86 and 89-90"),
    doc = list(term = "a random term", k = "k")),
  ## the same check on a fit that already used REML: the refit is not offered
  ## again (studies/prototypes/p10_few_groups_reml.R, results/p10_log.txt)
  few_groups_new_reml = list(
    check = "few groups, new units, REML fit", status = "WARN", detail = lmn_txt$few,
    cause = lmn_txt$few_new_cause,
    suggestion = paste("The fit already used restricted maximum likelihood (REML), the refit this",
                       "check otherwise suggests. In limn's prototype study of few groups",
                       "(studies/prototypes/), with REML and 8 groups of 6, 90% intervals for a new",
                       "unit's observation covered 0.892 (Poisson) and 0.903 (gaussian), with",
                       "simulation standard errors of 0.004 to 0.005. More groups is the sure cure."),
    reference = lmn_cite(lmn_ref$maas_hox, "pp. 86 and 89-90"),
    doc = list(term = "a random term", k = "k")),
  unseen_levels = list(
    check = "unseen levels", status = "NOTE",
    detail = "{n} of the rows name a level the fit has not seen.",
    cause = paste("Under groups = \"fitted\" a unit with no estimated effects can only be",
                  "predicted as a new one."),
    suggestion = paste("They are predicted as new units, with effects drawn from the fitted",
                       "distribution, and marked in the new_level column. If they were meant to be",
                       "units in the fit, check how their levels are spelt."),
    reference = NULL, doc = list(n = "Some")),
  boundary_draws = list(
    check = "boundary draws", status = "FAIL",
    detail = "{bad} of {S} draws of {what} are not finite.",
    cause = paste("Its estimate (SD {sd}) is at the edge of its range, where the data cannot",
                  "resolve it, and the fit did not hold it there."),
    suggestion = paste("Predict with uncertainty = \"fixed\", which holds the variances at their",
                       "estimates, or refit with boundary = \"avoid\" in illume::ilm_model(), which",
                       "keeps the variance off the edge with a small penalty."),
    reference = NULL, doc = list(bad = "Some", S = "the", what = "a variance", sd = "near zero")),
  nan_average = list(
    check = "draws not computed", status = "FAIL",
    detail = "{bad} of {S} draws of the {what} could not be computed.",
    cause = paste("Averaging over units failed where a draw's variance reaches {v}, as happens",
                  "when the variance is poorly determined, for instance with few groups."),
    suggestion = "Predict with uncertainty = \"fixed\", which holds the variances at their estimates.",
    reference = NULL, doc = list(bad = "Some", S = "the", what = "mean", v = "a large value")),
  nan_draws = list(
    check = "draws not computed", status = "FAIL",
    detail = "{bad} of {S} draws of the {what} could not be computed.",
    cause = "A draw's value overflowed or was undefined.",
    suggestion = "Predict with uncertainty = \"fixed\", which holds the variances at their estimates.",
    reference = NULL, doc = list(bad = "Some", S = "the", what = "mean")),
  dm_few_origins = list(
    check = "few origins", status = "WARN",
    detail = "The comparison has {T} origins over {H} horizons, fewer than about 100.",
    cause = paste("With so few origins the usual Diebold-Mariano test rejects too often when",
                  "forecasts overlap: in limn's simulation, up to 12% of the time at 5% with 50",
                  "origins or fewer (studies/findings/compare_dm_fixed.md)."),
    suggestion = paste("Use method = \"wpe\", the reference's fixed-smoothing test, which held its",
                       "size from 10 origins up in the same simulation; or compare over more origins."),
    reference = lmn_cite(lmn_ref$coroneo_iacone,
                         "section 4.2 (the test) and section 5 (its bandwidth)"),
    doc = list(T = "T", H = "several")),
  dm_negative_variance = list(
    check = "Diebold-Mariano variance", status = "WARN",
    detail = paste("The long-run variance truncated at lag {lag} was negative, so it was taken as 0",
                   "and the null hypothesis of equal accuracy was rejected."),
    cause = paste("A long-run variance truncated with equal weights can come out negative.",
                  "Diebold and Mariano's rule for that case is to take it as 0 and reject the",
                  "null, so the statistic is infinite and the p-value 0. They call a negative",
                  "estimate a rare event. In limn's own check against the simulation design of",
                  "Harvey, Leybourne and Newbold (1997, p. 284), it arose in 35% of samples with 8",
                  "origins at horizon 5 (lag 4), and in 17% with 16, and the rule took the",
                  "modified test's size at nominal 10% to 41.9% there, and to 29.7% with 8 origins",
                  "at horizon 3 (lag 2), against their Table 1's 12.9% and 18.1% (p. 285)",
                  "(studies/findings/compare_dm_negative.md)."),
    suggestion = paste("Read this rejection with care, or use method = \"wpe\", the default, which needs",
                       "no lag and held its size from 10 origins up in limn's simulation."),
    reference = lmn_cite(lmn_ref$diebold_mariano, "section 1.1"), doc = list(lag = "L")),
  nan_family = list(
    check = "draws not computed", status = "FAIL",
    detail = "{bad} of {S} draws of a new observation's distribution could not be computed.",
    cause = paste("In those draws the linear predictor or the dispersion ran so far out that the",
                  "family overflowed. That happens when a parameter sits at the edge of its range",
                  "(a variance or residual SD near zero, a correlation near 1) and the fit did not",
                  "hold it there, so its draws spread without limit."),
    suggestion = paste("Predict with uncertainty = \"fixed\", which holds the variances at their",
                       "estimates, or uncertainty = \"latent\", which holds the dispersion too."),
    reference = NULL, doc = list(bad = "Some", S = "the")),
  cv_interpolation = list(
    check = "folds and time", status = "WARN",
    detail = "The model has a correlation over time, and the folds are {kind}.",
    cause = paste("A test row then sits between training rows of its own series, whose cells",
                  "inform it, so the scores are those of interpolation, not forecasting, and",
                  paste0("flatter the model's forecasts ", lmn_txt$ours, ".")),
    suggestion = paste("To score forecasts, use rolling-origin folds, whose training rows all come",
                       "before the rows they score (the reference):",
                       "ilm_folds(type = \"rolling\", time = , group = , initial = , horizon = )."),
    reference = lmn_ref$fpp3_cv, doc = list(kind = "random or blocked")),
  cv_failed = list(
    check = "cross-validation", status = "FAIL",
    detail = "None of the {K} folds could be refitted and scored. The first fold's error: {error}",
    cause = paste("ilm_cv_refit() refits the model to each fold's training rows; for an illume fit",
                  "it is illume::ilm_refit(), which repeats the model's own call. A model it cannot",
                  "refit, or folds whose training rows cannot support it, leave nothing to score."),
    suggestion = paste("Check that every fold's training rows can support the model, and fit it",
                       "with a direct call to illume::ilm_model(): a fit made inside a function that",
                       "passes its arguments on with ... can leave a call ilm_refit() cannot repeat."),
    reference = NULL, doc = list(K = "K", error = "(the first error)")),
  ## ---- extension points (ilm_predict_draws(), ilm_cv_refit()) ----
  draws_shape = list(
    check = "parameter draws", status = "FAIL",
    detail = "The parameter draws for this fit are not in the shape limn needs: {what}.",
    cause = paste("Every prediction is computed from the draws ilm_predict_draws() returns, laid",
                  "out as illume::ilm_draws() lays them out. The method for this class of fit",
                  "returned something else."),
    suggestion = paste("If the fit's class comes from another package, its ilm_predict_draws()",
                       "method needs fixing there: ?ilm_predict_draws gives the shape. For an",
                       "illume fit this is a fault in limn; please report it."),
    reference = NULL, doc = list(what = "what is wrong")),
  trend = list(
    check = "deterministic trend", status = "WARN",
    detail = "The fixed part trends in `{time}` ({terms}).",
    cause = paste("Forecasts past the data assume that the trend's slope does not change (the",
                  "reference), so their intervals leave out any change in it",
                  paste0(lmn_txt$ours, ".")),
    suggestion = paste("For longer horizons a stochastic trend is safer (the reference): a random walk,",
                       "ar = ilm_rw1(~ {time} | {group}), lets the level move, and the intervals",
                       "widen with the horizon{instead}."),
    reference = lmn_ref$fpp3,
    doc = list(time = "time", terms = "its trend terms", group = "series", instead = "")),

  ## ---- totals (ilm_aggregate()) ----
  missing_totals = list(
    check = "totals not simulated", status = "FAIL",
    detail = "{bad} of {S} simulated totals are missing.",
    cause = paste("The family could not draw outcomes from some parameter draws, as happens when",
                  "a dispersion parameter sits at its limit and the fit did not hold it there (a",
                  "negative binomial whose size runs to infinity, say)."),
    suggestion = paste("If the data show no overdispersion beyond the model's random terms, refit",
                       "with family = \"poisson\" in illume::ilm_model(), whose totals are exact."),
    reference = NULL, doc = list(bad = "Some", S = "the")),

  ## ---- calibration (ilm_calibration() on a prediction) ----
  calib_narrow_counts = list(
    check = "coverage", status = "WARN", detail = lmn_txt$narrow, cause = lmn_txt$narrow_cause,
    suggestion = paste(paste0("Overdispersion is the usual cause for counts ", lmn_txt$ours, ":"),
                       "refit with family = \"nbinom\" in illume::ilm_model(), which lets the",
                       "variance exceed the mean."),
    reference = lmn_cite(lmn_ref$gneiting, "pp. 252-253"),
    doc = list(coverage = "a share", level = "L")),
  calib_narrow_dispersion = list(
    check = "coverage", status = "WARN", detail = lmn_txt$narrow, cause = lmn_txt$narrow_cause,
    suggestion = paste("Let the dispersion vary with the predictors through a dispformula in",
                       "illume::ilm_model(); with few groups and predictions for new units, refit",
                       "with reml = TRUE."),
    reference = lmn_cite(lmn_ref$gneiting, "pp. 252-253"),
    doc = list(coverage = "a share", level = "L")),
  calib_narrow = list(
    check = "coverage", status = "WARN", detail = lmn_txt$narrow, cause = lmn_txt$narrow_cause,
    suggestion = lmn_txt$remedies, reference = lmn_cite(lmn_ref$gneiting, "pp. 252-253"),
    doc = list(coverage = "a share", level = "L")),
  calib_wide = list(
    check = "coverage", status = "WARN",
    detail = paste("Too wide: {coverage} of outcomes fell inside the central {level}% of their",
                   "predictive distributions, against {level}% expected."),
    cause = "The predictions are less certain than the data show.",
    suggestion = paste("An overdispersed family, or a dispersion model the data do not need, widens",
                       paste0("every interval ", lmn_txt$ours, ":"), lmn_txt$remedies),
    reference = lmn_cite(lmn_ref$gneiting, "pp. 252-253"),
    doc = list(coverage = "a share", level = "L")),
  calib_bias = list(
    check = "bias", status = "WARN",
    detail = "Biased: outcomes are on average {size} {direction} their predicted means (SE {se}).",
    cause = "The mean structure is missing something the outcomes follow.",
    suggestion = "illume::ilm_remedies() and illume's residual checks point at which term.",
    reference = NULL, doc = list(size = "some amount", direction = "above or below", se = "s")),
  ## a shortfall in coverage alongside a biased mean: the bias alone can cause
  ## it, so the dispersion remedy waits until the mean is fixed
  calib_narrow_bias = list(
    check = "coverage", status = "WARN", detail = lmn_txt$narrow,
    cause = "A biased mean alone can cause this: the intervals are centred in the wrong place.",
    suggestion = paste("Fix the mean first, as above, then check the coverage again; widen the",
                       "intervals only if it still falls short."),
    reference = NULL, doc = list(coverage = "a share", level = "L")),
  calib_ok = list(
    check = "calibration", status = "OK",
    detail = paste("Consistent with calibration: coverage within sampling error at every level,",
                   "no bias in the mean, and the PITs consistent with uniform ({p})."),
    cause = "", suggestion = "", reference = NULL, doc = list(p = "p")),
  binary_intercept = list(
    check = "calibration intercept", status = "WARN",
    detail = paste("Off on average: {observed} of outcomes were 1, against {predicted} predicted;",
                   "the calibration intercept is {est} (95% interval {lo} to {hi})."),
    cause = "The probabilities are too high or too low overall.",
    suggestion = paste("A difference between the data the model was fitted to and the data scored,",
                       "such as a different rate of the outcome, can cause this (the reference, p. 2):",
                       lmn_txt$remedies),
    reference = lmn_cite(lmn_ref$van_calster, "pp. 2 and 4"),
    doc = list(observed = "a share", predicted = "another", est = "b", lo = "l", hi = "u")),
  binary_overconfident = list(
    check = "calibration slope", status = "WARN",
    detail = paste("Overconfident: the calibration slope is {est} (95% interval {lo} to {hi});",
                   "predicted probabilities are too extreme."),
    cause = "The model has fitted noise, so its probabilities stray too far from the average.",
    suggestion = paste("Use fewer terms, or penalised ones (the reference, p. 4): smooths in",
                       "illume::ilm_model() are penalised."),
    reference = lmn_cite(lmn_ref$van_calster, "pp. 2 and 4"), doc = list(est = "b", lo = "l", hi = "u")),
  ## a slope whose interval lies below 0 is more than too extreme
  binary_reversed = list(
    check = "calibration slope", status = "WARN",
    detail = paste("Reversed: the calibration slope is {est} (95% interval {lo} to {hi});",
                   "the predictions order the outcomes the wrong way."),
    cause = paste("This verdict is limn's own: the reference defines the calibration slope, with a",
                  "target of 1, and reads a slope below 1 as too extreme, but does not discuss one",
                  "below 0. Outcomes the model rates more likely turn out less likely; it may have",
                  "fitted noise, or the scored data may differ from the data it was fitted to",
                  paste0(lmn_txt$ours, ".")),
    suggestion = paste("Use fewer terms, or penalised ones: smooths in illume::ilm_model() are",
                       "penalised; and check that the scored data come from the same population",
                       "and coding."),
    reference = lmn_cite(lmn_ref$van_calster, "p. 4"), doc = list(est = "b", lo = "l", hi = "u")),
  binary_underconfident = list(
    check = "calibration slope", status = "WARN",
    detail = paste("Underconfident: the calibration slope is {est} (95% interval {lo} to {hi});",
                   "predicted probabilities are too moderate."),
    cause = "The probabilities stay too close to the average.",
    suggestion = paste("A penalty stronger than the data warrant can cause this",
                       paste0(lmn_txt$ours, ":"), lmn_txt$remedies),
    reference = lmn_cite(lmn_ref$van_calster, "p. 4"), doc = list(est = "b", lo = "l", hi = "u")),
  binary_ok = list(
    check = "calibration", status = "OK",
    detail = "Consistent with calibration: the calibration intercept includes 0 and the slope includes 1.",
    cause = "", suggestion = "", reference = NULL, doc = list()))

## Fill a template's {placeholders} from `values`.
#' @keywords internal
#' @noRd
lmn_fill <- function(txt, values) {
  ## numbers passed raw are counts and levels (rows, draws, origins, folds,
  ## groups; 90 for a 90% level): whole, grouped from 1,000. Estimates come
  ## already written, by lmn_fmt_sig() and the others.
  shown <- function(v) if (is.numeric(v))
    ifelse(is.finite(v) & v == round(v), lmn_fmt_count(v), lmn_fmt_sig(v)) else format(v)
  for (nm in names(values))
    txt <- gsub(paste0("{", nm, "}"), paste(shown(values[[nm]]), collapse = ", "), txt, fixed = TRUE)
  txt
}

## One diagnostic, filled in: a one-row data frame of its fields, with the
## values its words were made from kept in the list column `values`, so the
## check can be read as data as well as words.
#' @keywords internal
#' @noRd
lmn_diag <- function(id, ...) {
  e <- lmn_diagnostics[[id]]
  if (is.null(e)) stop("internal: no diagnostic '", id, "'.", call. = FALSE)
  v <- list(...)
  d <- data.frame(id = id, check = e$check, status = e$status, detail = lmn_fill(e$detail, v),
                  cause = lmn_fill(e$cause, v), suggestion = lmn_fill(e$suggestion, v),
                  reference = if (is.null(e$reference)) NA_character_ else e$reference,
                  stringsAsFactors = FALSE)
  d$values <- list(v)
  d
}

## No diagnostics: the zero-row frame of lmn_diag()'s columns.
#' @keywords internal
#' @noRd
lmn_no_checks <- function() lmn_diag("calib_ok", p = NA)[0, , drop = FALSE]

## Evaluate `expr`, keeping every limn diagnostic it signals (warnings pass
## on to the caller as before; an error still stops). Returns the value and
## the diagnostics as lmn_diag() rows, duplicates dropped.
#' @keywords internal
#' @noRd
lmn_collect <- function(expr) {
  got <- list()
  v <- withCallingHandlers(expr, lmn_diagnostic = function(cnd) {
    if (!inherits(cnd, "error")) got[[length(got) + 1L]] <<- cnd$diagnostic
  })
  list(value = v, checks = lmn_bind_checks(got))
}

## Diagnostics from several places as one frame: rows bound, and a
## diagnostic raised more than once with the same text kept once.
#' @keywords internal
#' @noRd
lmn_bind_checks <- function(x) {
  x <- Filter(function(d) !is.null(d) && nrow(d), x)
  if (!length(x)) return(lmn_no_checks())
  d <- do.call(rbind, x)
  d <- d[!duplicated(paste(d$id, d$detail, sep = "\r")), , drop = FALSE]
  rownames(d) <- NULL
  d
}

## What was found and why, with the reference that supports them: a reference
## follows the cause, not the remedy, since a remedy draws on it only where its
## own text says so.
#' @keywords internal
#' @noRd
lmn_found_text <- function(d) {
  parts <- c(d$detail, d$cause)
  txt <- paste(parts[nzchar(parts)], collapse = " ")
  if (!is.na(d$reference)) txt <- paste0(txt, " Reference: ", d$reference, ".")
  txt
}

## A diagnostic's text, as a warning or an error says it.
#' @keywords internal
#' @noRd
lmn_diag_text <- function(d) {
  parts <- c(lmn_found_text(d), d$suggestion)
  paste(parts[nzchar(parts)], collapse = " ")
}

## A diagnostic's text as a printed verdict: what was found and why (with its
## reference) on one line, the remedy on the next.
#' @keywords internal
#' @noRd
lmn_verdict_lines <- function(d) {
  unlist(lapply(seq_len(nrow(d)), function(i) {
    r <- d[i, ]
    c(lmn_found_text(r),
      if (nzchar(r$suggestion)) paste0("Remedy: ", r$suggestion))
  }))
}

## Signal a diagnostic as a warning or an error. The condition carries the
## filled fields, and its class names it (lmn_<id>, lmn_diagnostic).
#' @keywords internal
#' @noRd
lmn_signal <- function(id, ..., type = c("warning", "error")) {
  type <- match.arg(type)
  d <- lmn_diag(id, ...)
  cond <- structure(class = c(paste0("lmn_", id), "lmn_diagnostic", type, "condition"),
                    list(message = lmn_diag_text(d), call = NULL, diagnostic = d))
  if (type == "warning") warning(cond) else stop(cond)
}

## The help pages' text for a set of diagnostics, from the same entries: a
## roxygen section, for @eval.
#' @keywords internal
#' @noRd
lmn_diag_rd <- function(ids, title = "Checks") {
  items <- vapply(ids, function(id) {
    e <- lmn_diagnostics[[id]]
    d <- do.call(lmn_diag, c(list(id), e$doc))
    paste0("* **", e$check, "** (", tolower(e$status), "): ", lmn_diag_text(d))
  }, "", USE.NAMES = FALSE)
  c(paste0("@section ", title, ":"),
    paste("What each check means and what to do about it, in the words it is printed with.",
          "A reference supports what was found and why; a remedy draws on it only where it says",
          "so, and anything marked as limn's own reasoning is not the source's:"),
    "", items)
}
