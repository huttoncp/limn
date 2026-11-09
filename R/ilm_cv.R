#' Cross-validated scores of a fitted illume model
#'
#' Scores a model on data it was not fitted to, fold by fold: the model is
#' refitted to each fold's training rows ([ilm_cv_refit()], which for an
#' illume fit is [illume::ilm_refit()]), predicts the
#' fold's test rows, and each test row's predictive distribution is scored
#' against its outcome. Every row is scored by a model
#' that did not see it.
#'
#' @section Which unit a test row is predicted for:
#' `groups` says how the random effects and the correlation over time enter,
#' as in [ilm_predict()]. By default it follows the folds:
#' * **grouped folds** hold out whole units, so a test row's unit is new to the
#'   training fit, and it is predicted as a new unit (`"new"`);
#' * **random, blocked and rolling folds** keep a unit's other rows in the
#'   training data, so it is predicted as itself (`"fitted"`).
#'
#' Random or blocked folds on a model with a correlation over time score
#' interpolation, not forecasting: a test row sits between training rows of its
#' own series, so its neighbours' cells inform it. `ilm_cv()` warns, and names
#' rolling-origin folds, which score forecasts.
#'
#' @section Ex-ante and ex-post:
#' Rolling-origin folds score forecasts. A covariate whose future value would
#' not have been known at the origin, such as the weather on the day forecast,
#' makes the evaluation *ex-post*: it flatters the model. The folds record which
#' covariates are known in advance (`ilm_folds(known = )`); the result says
#' which evaluation it is.
#'
#' @param fit A fitted `"ilm_model"`.
#' @param data The data the folds index: the rows the model is refitted to and
#'   scored on.
#' @param folds Folds from [ilm_folds()] built on `data`.
#' @param rules The scores, as losses (lower is better): any of `"log"`, `"crps"`,
#'   `"brier"` (binary outcomes), `"pinball"` (at `probs`), `"interval"` (for the
#'   central `level` interval), `"wis"` (the weighted interval score, over the
#'   central intervals at levels 1 - `alpha`), `"squared"` (about the predictive mean) and
#'   `"absolute"` (about the median). Each is computed from the predictive
#'   mixture. [ilm_scores.ilm_pred()] defines each rule and gives its source.
#' @param groups How the random effects enter the predictions, as in
#'   [ilm_predict()]; `NULL`, the default, follows the folds (see *Which unit a
#'   test row is predicted for*).
#' @param uncertainty As in [ilm_predict()].
#' @param refit `"each"`: refit the model to every fold's training rows.
#' @param ncores Number of worker processes. The folds are shared out as
#'   workers become free; the results are the same whatever the number.
#' @param nsim Number of parameter draws for each fold's predictions.
#' @param seed Seed; each fold's draws use their own seed derived from it, so
#'   the results do not depend on `ncores`.
#' @param ... For the scores: `probs`, the pinball score's probabilities
#'   (default 0.1, 0.5, 0.9); `level`, the interval score's coverage
#'   (default 0.9); and `alpha`, the weighted interval score's intervals
#'   (default 0.02, 0.05, 0.1, ..., 0.9).
#' @return A data frame of class `"ilm_cv"`: one row per test row, with its
#'   fold, its row in `data`, for rolling folds its origin and horizon, and one
#'   column per score. A fold whose refit or prediction failed has `NA` scores
#'   and its error in `failed`. `print()` gives the mean of each score, by
#'   horizon for rolling folds; to compare two models' scores with standard
#'   errors that respect the folds, use [ilm_compare()].
#' @eval lmn_diag_rd(c("cv_interpolation", "cv_failed"))
#' @references Hyndman, R. J. and Athanasopoulos, G. (2021). *Forecasting:
#'   Principles and Practice*, 3rd edition. OTexts. Section 5.10, "Time series
#'   cross-validation" (online edition,
#'   <https://otexts.com/fpp3/toolbox.html#sec-tscv>, read 27 October 2026).
#' @seealso [ilm_folds()], [ilm_compare()], [ilm_plot_horizon()],
#'   [ilm_scores.ilm_pred()] for the scoring rules, and [ilm_cv_refit()] for
#'   how a fit is refitted to each fold.
#' @examples
#' set.seed(1)
#' d <- data.frame(site = factor(rep(1:12, each = 10)), t = rep(1:10, 12), x = rnorm(120))
#' d$y <- rpois(120, exp(0.5 + 0.3 * d$x + rnorm(12, 0, 0.4)[d$site]))
#' fit <- illume::ilm_model(y ~ x + (1 | site), data = d, family = "poisson",
#'                          verbose = FALSE)
#' folds <- ilm_folds(d, type = "grouped", group = "site", K = 4)
#' cv <- ilm_cv(fit, d, folds, nsim = 200)
#' cv
#' @export
ilm_cv <- function(fit, data, folds, rules = c("log", "crps"), groups = NULL,
                   uncertainty = c("joint", "fixed", "latent", "none"), refit = "each",
                   ncores = 1L, nsim = 1000L, seed = 1L, ...) {
  if (!inherits(fit, "ilm_model"))
    stop("`fit` must be a fitted ilm_model; it is ", class(fit)[1], ".", call. = FALSE)
  if (!inherits(folds, "ilm_folds"))
    stop("`folds` must come from ilm_folds().", call. = FALSE)
  if (!is.data.frame(data) || nrow(data) != attr(folds, "n"))
    stop("`data` must be the data frame the folds were built on (", attr(folds, "n"), " rows).",
         call. = FALSE)
  uncertainty <- match.arg(uncertainty)
  refit <- match.arg(refit, "each")
  ncores <- lmn_check_count(ncores, "ncores")
  nsim <- lmn_check_count(nsim, "nsim")
  type <- attr(folds, "type")
  if (is.null(groups)) groups <- if (identical(type, "grouped")) "new" else "fitted"
  vc <- illume::ilm_varcorr(fit)
  ## the checks cross-validation raises are warned as before, and kept on it,
  ## with those of each fold's prediction (silenced there, since they repeat)
  interp <- lmn_collect(if (!is.null(vc$latent) && type %in% c("random", "blocked"))
    lmn_signal("cv_interpolation", kind = type))$checks
  response <- lmn_response_name(fit)
  y <- data[[response]]
  score_args <- list(...)
  one <- function(k) {
    fo <- folds[[k]]
    base <- data.frame(fold = k, row = fo$test)
    if (!is.null(fo$origin)) { base$origin <- fo$origin; base$horizon <- fo$horizon }
    res <- tryCatch({
      fk <- suppressMessages(suppressWarnings(ilm_cv_refit(fit, data[fo$train, , drop = FALSE])))
      pk <- suppressWarnings(ilm_predict(fk, data[fo$test, , drop = FALSE], target = "response",
                                         groups = groups, uncertainty = uncertainty, nsim = nsim,
                                         seed = seed + k))
      sc <- do.call(lmn_pred_scores, c(list(pk, y[fo$test], rules = rules), score_args))
      structure(cbind(base, sc[attr(sc, "rules")], failed = NA_character_), checks = pk$checks)
    }, error = function(e) cbind(base, failed = conditionMessage(e)))
    res
  }
  ks <- seq_along(folds)
  out <- if (ncores == 1L) lapply(ks, one) else
    lmn_pool_lapply(ks, one, ncores, packages = lmn_method_packages(fit))
  checks <- lmn_bind_checks(c(list(interp), lapply(out, attr, "checks")))
  cols <- unique(unlist(lapply(out, names)))
  out <- do.call(rbind, lapply(out, function(x) { x[setdiff(cols, names(x))] <- NA; x[cols] }))
  score_cols <- setdiff(cols, c("fold", "row", "origin", "horizon", "failed"))
  if (all(!is.na(out$failed)))
    lmn_signal("cv_failed", K = length(folds), error = out$failed[1L], type = "error")
  evaluation <- NA_character_; covs <- character(0)
  if (identical(type, "rolling")) {
    covs <- lmn_unknown_covariates(fit, data, folds)
    evaluation <- if (length(covs)) "ex-post" else "ex-ante"
  }
  structure(out, class = c("ilm_cv", "data.frame"), rules = score_cols, type = type,
            groups = groups, uncertainty = uncertainty, refit = refit,
            evaluation = evaluation, unknown_covariates = if (length(covs)) covs,
            orientation = "lower is better", checks = checks)
}

## The fixed part's covariates whose values change over time within a series
## and are not named as known in advance: what makes a rolling evaluation
## ex-post. A covariate constant within each series (a clinic's size) is known.
#' @keywords internal
#' @noRd
lmn_unknown_covariates <- function(fit, data, folds) {
  ff <- if (!is.null(fit$fixed_formula)) fit$fixed_formula else stats::formula(fit)
  vars <- setdiff(all.vars(stats::delete.response(stats::terms(ff))),
                  c(attr(folds, "time"), attr(folds, "group"), attr(folds, "known")))
  vars <- intersect(vars, names(data))
  g <- attr(folds, "group")
  key <- if (is.null(g)) rep("1", nrow(data)) else as.character(data[[g]])
  vars[vapply(vars, function(v) any(tapply(data[[v]], key, function(z) length(unique(z)) > 1L)), TRUE)]
}

## The package's source directory, if it was loaded from source (pkgload's
## load_all()), or NULL if it is installed. An installed package also has a
## DESCRIPTION and an R/ directory; only it has Meta/package.rds.
#' @keywords internal
#' @noRd
lmn_source_path <- function(path) {
  if (file.exists(file.path(path, "Meta", "package.rds"))) return(NULL)
  if (file.exists(file.path(path, "DESCRIPTION")) && dir.exists(file.path(path, "R"))) path
}

## Folds shared out to a pool of worker processes, as each becomes free (B5):
## one TMB thread a worker, and each fold's own seed, so the results are the
## same whatever the number of workers.
#' @keywords internal
#' @noRd
lmn_pool_lapply <- function(x, fun, ncores, packages = character(0)) {
  cl <- parallel::makePSOCKcluster(min(ncores, length(x)))
  on.exit(parallel::stopCluster(cl), add = TRUE)
  lib <- .libPaths()
  ## an installed limn is loaded as such; one loaded from its source (while
  ## developing it) is loaded from the same source in each worker
  src <- lmn_source_path(getNamespaceInfo("limn", "path"))
  ## the set-up must not carry limn's namespace with it: a worker would then load
  ## limn to receive it, before the library paths are set
  ## then the packages whose methods the fit dispatches to (ilm_predict_draws(),
  ## ilm_cv_refit()), after limn, so that they register on limn's generics
  setup <- function(lib, src, packages) {
    .libPaths(lib)
    Sys.setenv(OMP_NUM_THREADS = "1")
    if (is.null(src)) suppressPackageStartupMessages(requireNamespace("limn", quietly = TRUE))
    else get("load_all", envir = asNamespace("pkgload"))(src, quiet = TRUE)
    for (p in packages)
      if (!suppressPackageStartupMessages(requireNamespace(p, quietly = TRUE)))
        stop("a worker could not load '", p, "', whose methods the fit needs.", call. = FALSE)
    invisible(NULL)
  }
  environment(setup) <- globalenv()
  parallel::clusterCall(cl, setup, lib, src, packages)
  parallel::clusterApplyLB(cl, x, fun)
}

#' @export
print.ilm_cv <- function(x, ...) {
  rules <- attr(x, "rules")
  type <- attr(x, "type")
  nf <- length(unique(x$fold)); failed <- unique(x$fold[!is.na(x$failed)])
  cat("<ilm_cv> ", lmn_fmt_count(nrow(x)), " scored rows over ", lmn_fmt_count(nf), " ", type,
      " folds; refit ", attr(x, "refit"), ", groups = \"", attr(x, "groups"), "\"\n", sep = "")
  ok <- is.na(x$failed)
  ## scores are estimates: three significant figures, trailing zeros kept
  cat("Mean scores (lower is better):\n")
  means <- colMeans(as.data.frame(x)[ok, rules, drop = FALSE], na.rm = TRUE)
  print(stats::setNames(lmn_fmt_sig(means), names(means)), quote = FALSE)
  if (!is.null(x$horizon)) {
    cat("By horizon:\n")
    by_h <- stats::aggregate(as.data.frame(x)[ok, rules, drop = FALSE], list(horizon = x$horizon[ok]),
                             mean, na.rm = TRUE)
    by_h$horizon <- lmn_fmt_data(by_h$horizon)
    for (r in rules) by_h[[r]] <- lmn_fmt_sig(by_h[[r]])
    print(by_h, row.names = FALSE, right = TRUE)
  }
  if (identical(attr(x, "evaluation"), "ex-post"))
    cat("Ex-post: ", paste(attr(x, "unknown_covariates"), collapse = ", "),
        " would not be known at the origin, which flatters the forecasts (see ilm_folds(known = )).\n",
        sep = "")
  if (length(failed))
    cat(lmn_fmt_count(length(failed)), " fold(s) failed and are left out: ", paste(failed, collapse = ", "),
        " (see the failed column).\n", sep = "")
  cat("For standard errors, compare models with ilm_compare().\n")
  invisible(x)
}
