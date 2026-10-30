#' Predictive distributions from a fitted illume model
#'
#' The engine every other prediction in `limn` goes through. For each row of
#' `newdata` it returns the distribution of the predicted quantity: the mean, a
#' new observation, or the linear predictor. That distribution carries the
#' fit's own uncertainty, drawn jointly from the Laplace approximation that
#' gives the fit its standard errors.
#'
#' @section How the random effects enter (`groups`):
#' A mixed model predicts different things for the same row, depending on the
#' question asked. So `groups` has no default for a model with random effects
#' or a correlation over time; it must be named. The choices:
#'
#' * `"fitted"`: the unit the row names, with that unit's own estimated effects.
#'   A level the fit has not seen is predicted as a new unit, and marked in the
#'   `new_level` column.
#' * `"new"`: a new, unseen unit. Its effects are drawn from the fitted
#'   distribution of effects, afresh in every parameter draw. Rows naming the
#'   same unit share its effects within a draw.
#' * `"typical"`: a typical unit, with its effects at zero. For a nonlinear
#'   link this is not the average unit.
#' * `"population"`: the average over units. For `target = "mean"`, the mean is
#'   averaged over the distribution of effects, as illume's
#'   `predict(groups = "population")` does. For a new observation, a randomly chosen unit's outcome is the same
#'   thing as a new unit's, so it equals `"new"`.
#'
#' Give one choice for every term, or a list per term, as in
#' `list(site = "fitted", ar = "new")`. The correlation over time is the term
#' `"ar"`. Smooth terms are part of the mean structure and always enter as
#' fitted.
#'
#' @section Over time:
#' With a correlation over time (`ilm_ar1()`, `ilm_car1()` or `ilm_rw1()`), a
#' seen series under `"fitted"` is continued from its fitted cells:
#' * **after its last cell**, a forecast carried forward from that cell;
#' * **between two cells**, drawn from the bridge between them;
#' * **before its first cell**, carried backward from it.
#'
#' A new series, or any series under `"new"`, is a fresh path. AR(1) and CAR(1)
#' paths start from their stationary distribution, and a random walk starts
#' from zero at the series' first time. Each series' new times are drawn in
#' sequence, each given the one before, so a draw is one path through them,
#' and rows naming the same series and time share one value, across scenarios
#' too. The model's own parameters are used throughout; nothing is refitted.
#'
#' For a mean or the linear predictor, the estimate at a seen series' new time
#' is the point forecast: the fitted cells carried forward with no new shocks,
#' at the estimates, as a Kalman filter's point forecast is. For a new series,
#' as for a new unit, it is the predictive median. `as.data.frame()` gives the
#' expected value beside it, as `mean`; for counts it can be well above the
#' estimate at long horizons (see [as.data.frame.ilm_pred()]).
#'
#' @section Covariates given as draws:
#' A covariate that is not known in advance, such as next month's
#' temperature, can be given as draws from its own forecast. `newdata` then
#' carries a `.draw` column: one block of rows per covariate draw, every block
#' listing the same rows in the same order, with only the drawn covariates
#' differing. Parameter draw s uses covariate draw ((s - 1) mod C) + 1 of the
#' C, so the prediction mixes over both. `nsim` is rounded up to a multiple of
#' C, with a message, so that each covariate draw is used equally often.
#'
#' The prediction has one row per row of a block. A drawn covariate shows as
#' `NA` in the rows where it varies, `meta` records C and the drawn columns,
#' and a drawn row's estimate is its predictive median, as a new unit's is.
#' Intervals conditional on a covariate's point forecast leave out the
#' uncertainty in forecasting it (fpp3, section 10.3). In limn's own study
#' (`studies/findings/covariate_draws.md`) they covered too little: 0.75 to 0.83
#' for 90% gaussian intervals, against about 0.90 with the covariate given as
#' draws.
#' [ilm_covariate_paths()] builds the blocks from an illume model of the
#' covariate, and [ilm_as_draws()] from paths simulated by another tool.
#'
#' @section Uncertainty:
#' * `"joint"` (the default): every parameter is drawn jointly, the variances
#'   included. In limn's prototype study of few groups (studies/prototypes/), it
#'   came closest to nominal coverage for a new unit's observation, though still
#'   below it: 0.88 for 90% intervals with 8 groups, against 0.85-0.87 for the
#'   other methods it tried, each with a simulation standard error of about
#'   0.01.
#' * `"fixed"`: the variance parameters are held at their estimates, and
#'   everything else is drawn.
#' * `"latent"`: every parameter is held at its estimate, the dispersion
#'   included, and only the random effects and the cells of a correlation over
#'   time are drawn, from their distribution given the parameters. For a
#'   gaussian model with a correlation over time this is a Kalman filter's
#'   forecast, which treats the parameters as known; limn's tests check it
#'   against `stats::KalmanForecast()`. ARIMA prediction intervals likewise
#'   leave out the variation in the parameter estimates (fpp3, section 9.8).
#' * `"none"`: the plug-in prediction at the estimates, with no parameter
#'   uncertainty. A new unit's effects, and a series' cells at times it was not
#'   fitted at (a forecast, a bridge or a backcast), are still drawn, and
#'   covariates given as draws are still mixed over; otherwise a new
#'   observation's distribution is the family's alone.
#'
#' @eval lmn_diag_rd(c("few_groups_means", "few_groups_new", "few_groups_new_reml", "unseen_levels", "trend", "boundary_draws", "nan_average", "nan_draws", "nan_family"))
#'
#' @section What is not here yet:
#' * **Multinomial outcomes** come later, with forecasts of their latent
#'   state.
#' * **Binomial outcomes with more than one trial** are predicted per trial.
#'
#' @param fit A fitted `"ilm_model"`.
#' @param newdata A data frame of rows to predict. Or a named list of data
#'   frames, one per scenario, stacked with a `.scenario` column.
#'   A `.draw` column gives covariates as draws; see *Covariates given as
#'   draws*.
#' @param target `"mean"`, the expected response; `"response"`, a new
#'   observation; or `"link"`, the linear predictor.
#' @param groups How the random effects enter; see *How the random effects
#'   enter*. Required when the model has random-effect or correlation-over-time
#'   terms.
#' @param uncertainty `"joint"`, `"fixed"`, `"latent"` or `"none"`; see
#'   *Uncertainty*.
#' @param set A named list of values to give every row, such as
#'   `list(dose = c(0, 10))`. Each row is repeated once per combination.
#' @param nsim Number of parameter draws.
#' @param seed Seed for the draws. The caller's random-number state is
#'   restored afterwards.
#' @param time,group For a correlation over time built from vectors, the new
#'   rows' times and groups, passed on to `illume::ilm_matrices()`.
#' @section The prediction object:
#' Read a prediction through [as.data.frame()][as.data.frame.ilm_pred] or
#' [summary()][summary.ilm_pred] where you can. Its elements are listed here so
#' that nothing in it is hidden. Every prediction, class `"ilm_pred"`, has:
#' * `kind`: `"response"` for a new observation's distribution, held as the
#'   family's parameters in each draw; `"mean"` or `"link"` for draws of the
#'   mean or the linear predictor; or `"sample"` for simulated outcomes, which
#'   is how [ilm_aggregate()] holds some totals.
#' * `rows`: a data frame with one row per predicted row: `newdata`'s columns
#'   (one block's, with drawn covariates `NA`; see *Covariates given as draws*),
#'   and `new_level`, `TRUE` where a row names a level the fit has not seen.
#'   A total or an average has its `by` columns and `n_rows` instead, and an
#'   average also `weight`.
#' * `weights`: the draws' weights, which sum to 1.
#' * `meta`: a list recording how the prediction was made. It has `target`,
#'   `groups` (one `"term = choice"` per term), `uncertainty`, `draws` (the
#'   number of parameter draws), `nsim`, `family`, `seed`, and `response`, the
#'   name of the response column. With covariates given as draws it also has
#'   `covariate_draws` and `drawn`. With a correlation over time it has `time`
#'   and `series`, the columns holding them.
#' * `checks`: a data frame of the checks the prediction raised (see *Checks*),
#'   one row each, with columns `id`, `check`, `status`, `detail`, `cause`,
#'   `suggestion`, `reference` and `values`. It has no rows when nothing was
#'   found.
#' * `origin`: each row's forecast origin, its series' last observed time, or
#'   `NA` for a row that is not ahead of its series' data.
#'
#' A prediction straight from `ilm_predict()` or [ilm_forecast()] also has
#' `new_unit`, `TRUE` where a row's effects are drawn as a new unit's or a new
#' series'. A total or an average has `what` instead, saying how it was taken,
#' and `meta` records the same in `aggregate`.
#'
#' What else it holds depends on `kind`:
#' * `"response"`: `dist`, the family's distribution functions as
#'   `illume::ilm_dist()` gives them; `lp`, a named list of rows-by-draws
#'   matrices of linear predictors; `par`, the family's other parameters in
#'   each draw; and `pvec`, any parameter that is a vector in each draw, such as
#'   an ordered outcome's thresholds.
#' * `"mean"` and `"link"`: `draws`, a rows-by-draws matrix of the predicted
#'   values, and `estimate`, the value at the estimates (a new unit's or a
#'   drawn row's predictive median; see *Value*).
#' * `"sample"`: `draws`, a rows-by-samples matrix of simulated outcomes;
#'   `estimate`, their median; and `dist`, the family and whether it is
#'   discrete.
#'
#' @return A predictive distribution, class `"ilm_pred"`: see *The prediction
#'   object*. [summary()][summary.ilm_pred] and
#'   [as.data.frame()][as.data.frame.ilm_pred] give point estimates and
#'   intervals.
#'   - For `"mean"` and `"link"`, the point estimate is the value at the
#'     estimates, the number `predict()` gives.
#'   - For `"response"`, it is the predictive median, with the mean beside it.
#' @references Hyndman, R. J. and Athanasopoulos, G. (2021). *Forecasting:
#'   Principles and Practice*, 3rd edition. OTexts. Section 9.8, "Forecasting"
#'   (ARIMA), and section 10.3, "Forecasting" (dynamic regression) (online
#'   edition, <https://otexts.com/fpp3/>, read 17 October 2026).
#' @seealso [ilm_forecast()] for forecasts; [ilm_aggregate()] and
#'   [ilm_average()] for totals and averages over rows; [ilm_pit()] and
#'   [ilm_scores.ilm_pred()] to check it against outcomes; [ilm_plot_fan()].
#' @examples
#' set.seed(1)
#' d <- data.frame(site = factor(rep(1:24, each = 5)), x = rnorm(120))
#' d$y <- rpois(120, exp(0.5 + 0.4 * d$x + rnorm(24, 0, 0.5)[d$site]))
#' fit <- illume::ilm_model(y ~ x + (1 | site), data = d, family = "poisson",
#'                          verbose = FALSE)
#' nd <- data.frame(site = c("1", "2", "new"), x = 0)
#' as.data.frame(ilm_predict(fit, nd, groups = "fitted", nsim = 200))
#' as.data.frame(ilm_predict(fit, nd, target = "response", groups = "new",
#'                           nsim = 200))
#' @export
ilm_predict <- function(fit, newdata, target = c("mean", "response", "link"), groups,
                        uncertainty = c("joint", "fixed", "latent", "none"),
                        set = NULL,
                        nsim = 1000L, seed = 1L, time = NULL, group = NULL) {
  if (!inherits(fit, "ilm_model"))
    stop("`fit` must be a fitted ilm_model; it is ", class(fit)[1], ".", call. = FALSE)
  target <- match.arg(target)
  uncertainty <- match.arg(uncertainty)
  nsim <- lmn_check_count(nsim, "nsim")
  nd <- lmn_scenarios(newdata)
  if (!is.null(set)) nd <- lmn_expand_set(nd, set)

  dist <- illume::ilm_dist(fit)
  if (identical(dist$family, "multinomial"))
    stop("predictions for a multinomial outcome are not in limn yet.", call. = FALSE)
  if (target == "mean" && !is.null(dist$categories))
    stop("an ordered outcome's mean is a probability for each category, not one ",
         "number. Predict with target = \"response\": its distribution is those ",
         "probabilities, and its median is the most central category.", call. = FALSE)
  vc <- illume::ilm_varcorr(fit)
  re_terms <- names(Filter(function(m) !identical(attr(m, "kind"), "smooth"), vc$re))
  has_ar <- !is.null(vc$latent)
  spec <- lmn_groups_spec(if (missing(groups)) NULL else groups, re_terms, has_ar)

  m <- illume::ilm_matrices(fit, nd, time = time, group = group)
  ## covariates given as draws: blocks of the same rows, one block
  ## per draw, paired with the parameter draws in turn
  cb <- lmn_covariate_blocks(nd, m)
  if (!is.null(cb)) nsim <- lmn_round_nsim(nsim, cb$C)
  ## the variances' draws matter only when a unit's effects are drawn afresh or
  ## integrated over, which a population mean also does
  natural <- lmn_needs_sampling(spec, m, "response")
  dr <- ilm_predict_draws(fit, nsim = nsim, seed = seed, uncertainty = uncertainty,
                          natural = natural)
  lmn_check_draws(dr, nsim, uncertainty, natural)
  if (natural) lmn_check_natural(dr, vc, names(spec))
  sample_needed <- lmn_needs_sampling(spec, m, target) || !is.null(cb)
  if (uncertainty == "none" && sample_needed) dr <- lmn_repeat_draw(dr, nsim)

  lpd <- lmn_pick_blocks(lmn_assemble(fit, m, dr, spec, dist, target, seed), cb)
  ## the checks the prediction raises are warned as before, and kept on it
  ck <- list(lmn_collect(lmn_warn_few_groups(spec, m, dr, target,
                                             reml = isTRUE(fit$reml)))$checks,
             lmn_collect(lmn_warn_trend(fit, newdata, nd, m, vc))$checks)
  rows <- if (is.null(cb)) nd else cb$rows
  rows$new_level <- lpd$new_level
  if (any(rows$new_level)) ck <- c(ck, list(lmn_diag("unseen_levels", n = sum(rows$new_level))))
  meta <- list(target = target, groups = paste(names(spec), spec, sep = " = "),
               uncertainty = uncertainty, draws = ncol(dr$draws), nsim = nsim,
               family = dist$family,
               seed = lmn_seed_record(seed, "draws"))
  if (!is.null(cb)) { meta$covariate_draws <- cb$C; meta$drawn <- cb$drawn }
  checks <- lmn_bind_checks(ck)
  fcm <- lmn_forecast_meta(fit, m, nd, cb)
  meta$time <- fcm$time; meta$series <- fcm$series
  ## the response column, so a score's scale can be read from the history
  meta$response <- all.vars(stats::formula(fit)[[2]])[1]

  if (target == "response") {
    lmn_check_values(lpd$lp$mu, lpd, "linear predictor")
    pred <- lmn_new_pred(dist, lpd$lp, par = lpd$par, pvec = lpd$pvec, rows = rows, meta = meta)
    lmn_check_family(pred)
    pred$checks <- checks; pred$origin <- fcm$origin; pred$new_unit <- lpd$new_unit
    return(pred)
  }
  ## the mean or the linear predictor: draws of it, and its value at the estimates
  ## at the estimates, a seen series' new times carry its cells forward with no
  ## new shocks: the point forecast
  lp0 <- lmn_assemble(fit, m, lmn_mode_draw(dr, vc), spec, dist, target, seed, noise = FALSE)
  est <- if (target == "link") lp0$lp$mu[, 1] else lmn_mean_values(lp0, dist)[, 1]
  if (!is.null(cb)) est <- est[cb$idx[, 1]]
  vals <- if (target == "link") lpd$lp$mu else lmn_mean_values(lpd, dist)
  lmn_check_values(vals, lpd, if (target == "link") "linear predictor" else "mean")
  ## a new unit or series, or a row whose covariates are drawn, has no value at
  ## the estimates: its estimate is the predictive median
  med <- lpd$new_unit | (if (is.null(cb)) FALSE else cb$drawn_row)
  if (any(med)) est[med] <- apply(vals[med, , drop = FALSE], 1L, stats::median)
  pred <- lmn_new_draw_pred(vals, est, kind = target, rows = rows, meta = meta)
  pred$checks <- checks; pred$origin <- fcm$origin; pred$new_unit <- lpd$new_unit
  pred
}

## ---- pieces of the engine --------------------------------------------------------

## For a correlation over time: which columns of the new rows hold its time and
## its series, and each row's origin, the series' last observed time, for rows
## after it. The fit's cells run on one grid of times for every series, so a
## series that ends early has cells past its data: the observed times are read
## from the fit's own rows, whose response was seen. Rows up to a series' last
## observed time, in a series the fit has not seen, or with no correlation over
## time have no origin (NA).
#' @keywords internal
#' @noRd
lmn_forecast_meta <- function(fit, m, nd, cb = NULL) {
  ar <- m$ar
  n <- if (is.null(cb)) nrow(nd) else nrow(cb$rows)
  none <- list(time = NULL, series = NULL, origin = rep(NA_real_, n))
  if (is.null(ar) || !nrow(ar)) return(none)
  same_num <- function(z) {
    v <- suppressWarnings(as.numeric(if (is.factor(z)) as.character(z) else z))
    length(v) == nrow(ar) && !anyNA(v) && isTRUE(all(v == as.numeric(ar$time)))
  }
  same_chr <- function(z) length(z) == nrow(ar) && isTRUE(all(as.character(z) == as.character(ar$group)))
  tcol <- names(nd)[vapply(nd, same_num, TRUE)][1]
  gcol <- if (all(is.na(ar$group))) NA else names(nd)[vapply(nd, same_chr, TRUE)][1]
  mf <- stats::model.frame(fit)
  if (!is.na(tcol) && tcol %in% names(mf) && (is.na(gcol) || gcol %in% names(mf))) {
    key <- function(g) if (is.na(gcol)) rep("", length(g)) else as.character(g)
    last <- tapply(as.numeric(mf[[tcol]]), key(mf[[if (is.na(gcol)) tcol else gcol]]), max)
    lt <- as.vector(last[key(ar$group)])
    origin <- ifelse(!ar$new_group & !is.na(lt) & as.numeric(ar$time) > lt, lt, NA_real_)
  } else {                       # the fit's columns not found: the cells alone
    ahead <- is.na(ar$cell) & !is.na(ar$prev_cell) & is.na(ar$next_cell) & !ar$new_group
    origin <- ifelse(ahead, as.numeric(ar$time) - ar$dt_prev, NA_real_)
  }
  if (!is.null(cb)) origin <- origin[cb$idx[, 1]]
  list(time = if (is.na(tcol)) NULL else tcol, series = if (is.na(gcol)) NULL else gcol, origin = origin)
}


## Covariates given as draws: `newdata`
## carries a `.draw` column, and every block of it lists the same rows in the
## same order, with only the covariates differing. The rows are checked by what
## the model makes of them: the scenario, each random term's unit, and the
## series and time of a correlation over time. NULL without a `.draw` column.
#' @keywords internal
#' @noRd
lmn_covariate_blocks <- function(nd, m) {
  if (is.null(nd$.draw)) return(NULL)
  if (anyNA(nd$.draw)) stop("`.draw` has missing values; every row needs its draw.", call. = FALSE)
  labs <- unique(nd$.draw)
  pos <- split(seq_len(nrow(nd)), factor(nd$.draw, levels = labs))
  n <- lengths(pos)
  if (any(n != n[1]))
    stop("every `.draw` block must list the same rows: block ", format(labs[which(n != n[1])[1]]),
         " has ", n[n != n[1]][1], " rows and block ", format(labs[1]), " has ", n[1], ".",
         call. = FALSE)
  idx <- matrix(unlist(pos, use.names = FALSE), n[1])
  ids <- list(.scenario = nd$.scenario)
  for (t in names(m$re)) ids[[paste0("the units of '", t, "'")]] <- m$re[[t]]$level
  if (!is.null(m$ar)) {
    ids[["the series"]] <- m$ar$group
    ids[["the times"]] <- m$ar$time
  }
  for (k in names(ids)) {
    v <- ids[[k]]
    if (is.null(v)) next
    same <- apply(idx, 2L, function(i) identical(v[i], v[idx[, 1]]))
    if (!all(same))
      stop("every `.draw` block must list the same rows in the same order, with only the ",
           "covariates differing; block ", format(labs[which(!same)[1]]), " differs from block ",
           format(labs[1]), " in ", if (k == ".scenario") "`.scenario`" else k, ".", call. = FALSE)
  }
  rows <- nd[idx[, 1], setdiff(names(nd), ".draw"), drop = FALSE]
  drawn_row <- rep(FALSE, nrow(rows)); drawn <- character(0)
  for (v in names(rows)) {
    col <- nd[[v]]; varies <- rep(FALSE, nrow(rows))
    for (c in seq_len(ncol(idx))[-1L])
      varies <- varies | !lmn_same_values(col[idx[, c]], col[idx[, 1]])
    if (!any(varies)) next
    drawn <- c(drawn, v); drawn_row <- drawn_row | varies
    rows[[v]][varies] <- NA
  }
  rownames(rows) <- NULL
  list(C = length(labs), idx = idx, rows = rows, drawn = drawn, drawn_row = drawn_row)
}

## Element by element: are two vectors equal, missing values included?
#' @keywords internal
#' @noRd
lmn_same_values <- function(a, b) {
  a <- as.character(a); b <- as.character(b)
  (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
}

## Each covariate draw is used equally often, so the number of parameter draws
## is rounded up to a multiple of the number of covariate draws, and it says so.
#' @keywords internal
#' @noRd
lmn_round_nsim <- function(nsim, C) {
  S <- as.integer(ceiling(nsim / C) * C)
  if (S != nsim)
    message("nsim is rounded up from ", nsim, " to ", S, ", a multiple of the ", C,
            " covariate draws, so that each is used equally often.")
  S
}

## Parameter draw s takes its rows from covariate draw ((s - 1) mod C) + 1.
#' @keywords internal
#' @noRd
lmn_pick_blocks <- function(lpd, cb) {
  if (is.null(cb)) return(lpd)
  S <- ncol(lpd$lp[[1]]); cs <- (seq_len(S) - 1L) %% cb$C + 1L
  pick <- function(M) {
    out <- matrix(0, nrow(cb$idx), S)
    for (c in seq_len(cb$C)) {
      s <- which(cs == c)
      out[, s] <- M[cb$idx[, c], s, drop = FALSE]
    }
    out
  }
  lpd$lp <- lapply(lpd$lp, pick)
  lpd$v_int <- pick(lpd$v_int)
  for (k in c("new_level", "new_unit", "drawn")) lpd[[k]] <- lpd[[k]][cb$idx[, 1]]
  lpd
}

## A seed derived from the user's, for one step's random draws; NULL stays NULL.
#' @keywords internal
#' @noRd
lmn_sub_seed <- function(seed, k) if (is.null(seed)) NULL else seed + k

## A named list of data frames is a set of scenarios, stacked and labelled.
#' @keywords internal
#' @noRd
lmn_scenarios <- function(newdata) {
  if (is.data.frame(newdata)) return(newdata)
  if (!is.list(newdata) || is.null(names(newdata)) || any(!nzchar(names(newdata))) ||
      !all(vapply(newdata, is.data.frame, TRUE)))
    stop("`newdata` must be a data frame, or a named list of data frames, one ",
         "per scenario.", call. = FALSE)
  out <- do.call(rbind, lapply(names(newdata), function(k)
    cbind(.scenario = k, newdata[[k]], stringsAsFactors = FALSE)))
  rownames(out) <- NULL
  out
}

## Every row repeated once per combination of the values in `set`.
#' @keywords internal
#' @noRd
lmn_expand_set <- function(nd, set) {
  if (!is.list(set) || is.null(names(set)) || any(!nzchar(names(set))))
    stop("`set` must be a named list of values, such as list(dose = c(0, 10)).",
         call. = FALSE)
  grid <- expand.grid(set, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  out <- nd[rep(seq_len(nrow(nd)), each = nrow(grid)), , drop = FALSE]
  for (v in names(grid)) out[[v]] <- rep(grid[[v]], times = nrow(nd))
  rownames(out) <- NULL
  out
}

## The choice for every term, from one string or a list per term.
#' @keywords internal
#' @noRd
lmn_groups_spec <- function(groups, re_terms, has_ar) {
  terms <- c(re_terms, if (has_ar) "ar")
  choices <- c("fitted", "new", "typical", "population")
  if (!length(terms)) return(stats::setNames(character(0), character(0)))
  if (is.null(groups))
    stop("`groups` is needed: this model has ",
         paste(c(if (length(re_terms)) paste0("random effects (", paste(re_terms, collapse = ", "), ")"),
                 if (has_ar) "a correlation over time (\"ar\")"), collapse = " and "),
         ", and a row's prediction depends on which unit it is for. Name one of\n",
         "  \"fitted\":     the unit the row names, with its estimated effects;\n",
         "  \"new\":        a new, unseen unit, its effects drawn from their fitted distribution;\n",
         "  \"typical\":    a typical unit, with its effects at zero;\n",
         "  \"population\": the average over units (for a new observation, the same as \"new\").\n",
         "Give one for every term, or one per term: list(", terms[1], " = \"fitted\", ...).",
         call. = FALSE)
  if (is.character(groups) && length(groups) == 1L && is.null(names(groups)))
    groups <- stats::setNames(rep(groups, length(terms)), terms)
  groups <- unlist(groups)
  bad <- setdiff(unique(groups), choices)
  if (length(bad))
    stop("`groups` must be one of ", paste0("\"", choices, "\"", collapse = ", "),
         "; not ", paste0("\"", bad, "\"", collapse = ", "), ".", call. = FALSE)
  unknown <- setdiff(names(groups), terms)
  if (length(unknown))
    stop("the model has no term ", paste0("'", unknown, "'", collapse = ", "),
         "; its terms are ", paste0("'", terms, "'", collapse = ", "), ".", call. = FALSE)
  missing_terms <- setdiff(terms, names(groups))
  if (length(missing_terms))
    stop("say how ", paste0("'", missing_terms, "'", collapse = ", "),
         " enter too: `groups` must cover every term.", call. = FALSE)
  groups[terms]
}

## Parameter draws: joint, with the variances held, with every parameter held
## (the effects and cells alone drawn), or the estimates alone.
## Transforming each draw's variances to their natural scale is the slow part
## of drawing, so it is asked for only when `natural`; variances held at their
## estimates are ilm_varcorr()'s in every draw.
#' @keywords internal
#' @noRd
lmn_param_draws <- function(fit, uncertainty, nsim, seed, vc, natural = TRUE) {
  if (uncertainty == "none")
    return(lmn_mode_draw(illume::ilm_draws(fit, nsim = 1L, seed = seed, natural = FALSE), vc))
  given <- c(joint = "none", fixed = "theta", latent = "parameters")[[uncertainty]]
  dr <- illume::ilm_draws(fit, nsim = nsim, seed = seed, given = given,
                          natural = natural && uncertainty == "joint")
  if (natural && uncertainty != "joint")
    dr$natural <- lmn_rep_natural(lmn_vc_natural(vc), ncol(dr$draws))
  dr
}

## The natural-scale variances ilm_varcorr() reports (the same transform as the
## draws'), as arrays of one draw.
#' @keywords internal
#' @noRd
lmn_vc_natural <- function(vc) {
  nat <- list(re = lapply(vc$re, function(M) array(unclass(M), c(dim(M), 1L))))
  if (!is.null(vc$latent)) {
    lt <- vc$latent
    for (k in intersect(c("Sigma", "var_per_time"), names(lt)))
      lt[[k]] <- array(lt[[k]], c(dim(as.matrix(lt[[k]])), 1L))
    nat$latent <- lt
  }
  nat
}

## One draw's natural-scale variances, repeated for S draws.
#' @keywords internal
#' @noRd
lmn_rep_natural <- function(nat, S) {
  rep_arr <- function(a) array(rep(a, S), c(dim(a)[1:2], S))
  nat$re <- lapply(nat$re, rep_arr)
  if (!is.null(nat$latent)) {
    for (k in intersect(c("Sigma", "var_per_time"), names(nat$latent)))
      nat$latent[[k]] <- rep_arr(nat$latent[[k]])
    if (!is.null(nat$latent$rho)) nat$latent$rho <- rep(nat$latent$rho[1], S)
  }
  nat
}

## One "draw" at the estimates: the mode, with ilm_varcorr()'s variances.
#' @keywords internal
#' @noRd
lmn_mode_draw <- function(dr, vc) {
  dr$draws <- matrix(dr$mode, ncol = 1L, dimnames = list(names(dr$mode), NULL))
  nat <- lmn_vc_natural(vc)
  dr$natural$re <- nat$re
  dr$natural$latent <- nat$latent
  dr
}

## The estimates repeated, for plug-in predictions that still draw new units.
#' @keywords internal
#' @noRd
lmn_repeat_draw <- function(dr, S) {
  dr$draws <- dr$draws[, rep(1L, S), drop = FALSE]
  dr$natural <- lmn_rep_natural(dr$natural, S)
  dr
}

## With few groups, for means as well as new observations,
## intervals cover below their nominal level (studies/findings/engine_coverage.md).
## For means and linear predictors, whatever the unit, they cover as z-intervals
## do; no remedy is measured yet. For a new observation the shortfall is for
## units whose effects are drawn or averaged over, and REML narrows it
## (prototype P10, studies/prototypes/p10_few_groups_reml.R). On a fit that
## already used REML (`reml`, from fit$reml) that refit is not suggested again.
#' @keywords internal
#' @noRd
lmn_warn_few_groups <- function(spec, m, dr, target, reml = FALSE) {
  for (t in setdiff(names(spec), "ar")) {
    ch <- spec[[t]]
    uses_var <- (ch == "new") || (ch == "population" && target != "link") ||
      (ch == "fitted" && any(m$re[[t]]$new_group))
    if (target == "response" && !uses_var) next
    k <- length(unique(dr$map$level[dr$map$block == "bvec" & dr$map$term == t]))
    if (k >= 20L) next
    id <- if (target != "response") "few_groups_means" else
      if (reml) "few_groups_new_reml" else "few_groups_new"
    lmn_signal(id, term = paste0("'", t, "'"), k = k)
  }
  invisible(TRUE)
}

## fpp3 section 10.4: a trend in time in the fixed part assumes its slope never
## changes, so forecasts' intervals leave out a change in it. A stochastic trend
## (a random walk, illume's ilm_rw1()) lets the level move, and they widen with
## the horizon. Warned when rows lie past the data and there is no random walk.
## The time variable comes from the correlation over time, or from
## ilm_future()'s rows.
#' @keywords internal
#' @noRd
lmn_warn_trend <- function(fit, newdata, nd, m, vc) {
  if (identical(vc$latent$type, "rw1")) return(invisible(FALSE))
  vars <- if (!is.null(m$ar)) attr(illume::ilm_cells(fit), "vars") else NULL
  tv <- if (!is.null(vars)) vars[["time"]] else attr(newdata, "time")
  gv <- if (!is.null(vars)) vars[["group"]] else attr(newdata, "group")
  if (is.null(tv)) return(invisible(FALSE))
  ahead <- if (!is.null(m$ar)) any(is.na(m$ar$cell) & is.na(m$ar$next_cell) & !is.na(m$ar$prev_cell))
    else !is.null(nd$horizon) && any(nd$horizon > 0, na.rm = TRUE)
  if (!ahead) return(invisible(FALSE))
  labs <- attr(stats::terms(stats::formula(fit)), "term.labels")
  labs <- labs[!grepl("|", labs, fixed = TRUE)]
  trend <- labs[vapply(labs, function(l) lmn_uses_outside(str2lang(l), tv), TRUE)]
  if (!length(trend)) return(invisible(FALSE))
  lmn_signal("trend", time = tv, terms = paste(trend, collapse = ", "),
             group = if (is.null(gv)) "series" else gv,
             instead = if (is.null(vc$latent)) "" else
               paste0(" (in place of the model's ", vc$latent$type,
                      " term: a model has one correlation over time)"))
  invisible(TRUE)
}

## Does `expr` use the variable `v` anywhere but inside a periodic term?
#' @keywords internal
#' @noRd
lmn_uses_outside <- function(expr, v) {
  if (is.name(expr)) return(identical(as.character(expr), v))
  if (!is.call(expr)) return(FALSE)
  fn <- expr[[1]]
  fname <- if (is.name(fn)) as.character(fn)
    else if (is.call(fn) && identical(fn[[1]], as.name("::"))) as.character(fn[[3]]) else ""
  if (fname %in% c("ilm_fourier", "ilm_cyclic", "iml_fourier", "iml_cyclic")) return(FALSE)
  any(vapply(as.list(expr)[-1], lmn_uses_outside, TRUE, v = v))
}

## A draw whose value is NaN is a numerical failure, not a value, so the
## prediction stops and says where. Inf still sorts correctly into an interval.
#' @keywords internal
#' @noRd
lmn_check_values <- function(vals, lpd, what) {
  bad <- colSums(is.nan(vals)) > 0
  if (!any(bad)) return(invisible(TRUE))
  v <- lpd$v_int[, bad, drop = FALSE]
  if (any(v > 0))
    lmn_signal("nan_average", bad = sum(bad), S = length(bad), what = what,
               v = lmn_fmt_sig(max(v)), type = "error")
  lmn_signal("nan_draws", bad = sum(bad), S = length(bad), what = what, type = "error")
}

## A new observation's distribution must be computable in every draw. A draw
## whose linear predictor or dispersion runs far out (from a parameter at the
## edge of its range that the fit did not hold) can overflow the family, which
## then has no quantiles; the prediction stops and says so.
#' @keywords internal
#' @noRd
lmn_check_family <- function(pred) {
  i <- seq_len(nrow(pred$lp[[1]]))
  q <- tryCatch(suppressWarnings(lmn_draw_quantiles(pred, i, 0.5)), error = function(e) NULL)
  ## is.na(), not is.nan(): an ordered outcome's quantile is a category's
  ## integer index, so a draw that cannot be computed gives NA_integer_
  bad <- if (!is.null(q)) colSums(is.na(q)) > 0 else {
    ## some families' quantile functions stop, rather than return NaN, on a
    ## parameter that is not a number (a zero part's, for one): then draw by
    ## draw, and a draw that stops counts as one that cannot be computed
    vapply(seq_along(pred$weights), function(s) {
      a <- if (length(pred$pvec)) lmn_draw_args(pred, i, s) else lmn_block_args(pred, i, s)
      qs <- tryCatch(suppressWarnings(pred$dist$q(rep(0.5, length(i)), a$lp, a$par)),
                     error = function(e) NA)
      any(is.na(qs))
    }, TRUE)
  }
  if (any(bad)) lmn_signal("nan_family", bad = sum(bad), S = length(bad), type = "error")
  invisible(TRUE)
}

## A variance whose draws are not finite has its estimate at the edge of its
## range, in a direction the fit did not hold; it cannot give a new unit's
## effect. Checked for the terms a prediction draws or integrates over.
#' @keywords internal
#' @noRd
lmn_check_natural <- function(dr, vc, terms) {
  bad_draws <- function(a) {                      # per draw: any value not finite
    d <- dim(a)
    if (is.null(d)) !is.finite(a) else !apply(array(is.finite(a), d), length(d), all)
  }
  S <- ncol(dr$draws)
  stop_at <- function(what, bad, sds)
    lmn_signal("boundary_draws", bad = bad, S = S, what = what,
               sd = paste(lmn_fmt_sig(sds), collapse = ", "), type = "error")
  for (t in setdiff(terms, "ar")) {
    bad <- sum(bad_draws(dr$natural$re[[t]]))
    if (bad) stop_at(paste0("the covariance of '", t, "'"), bad, sqrt(diag(as.matrix(vc$re[[t]]))))
  }
  if ("ar" %in% terms) {
    lt <- dr$natural$latent
    v <- if (!is.null(lt$var_per_time)) "var_per_time" else "Sigma"
    bad <- sum(bad_draws(lt[[v]]) | (if (is.null(lt$rho)) FALSE else bad_draws(lt$rho)))
    if (bad) stop_at("the correlation over time", bad, sqrt(as.numeric(vc$latent[[v]])))
  }
  invisible(TRUE)
}

## Will any row need a unit's effects, or a latent value, drawn afresh?
#' @keywords internal
#' @noRd
lmn_needs_sampling <- function(spec, m, target) {
  if (!length(spec)) return(FALSE)
  for (t in names(spec)) {
    ch <- spec[[t]]
    if (ch == "new" || (ch == "population" && target == "response")) return(TRUE)
    newlev <- if (t == "ar") m$ar$new_group else m$re[[t]]$new_group
    if (ch == "fitted" && any(newlev)) return(TRUE)
    ## a seen series at a time off its cells: a forecast, a bridge or a backcast
    if (t == "ar" && ch == "fitted" && any(is.na(m$ar$cell))) return(TRUE)
  }
  FALSE
}

## The linear predictors for every row and draw, with the family's parameters.
## Effects enter by `spec`: fitted draws, zero, fresh draws for a new unit
## (shared by rows naming the same unit), or, for a population mean, a
## variance for ilm_normal_expect() to integrate over.
#' @keywords internal
#' @noRd
lmn_assemble <- function(fit, m, dr, spec, dist, target, seed, noise = TRUE) {
  D <- dr$draws; map <- dr$map; S <- ncol(D); n <- nrow(m$X)
  key <- paste(map$block, map$term, map$level, map$dim, sep = "\r")
  ## a design column's draw, matched by name; the zero part's and the dispersion
  ## model's draws may carry their part's prefix ("zi:x", "disp:x") where the
  ## design's columns do not
  by_name <- function(block, names, prefix = NULL) {
    terms <- map$term[map$block == block]
    k <- match(names, terms)
    if (!is.null(prefix)) k[is.na(k)] <- match(paste0(prefix, names[is.na(k)]), terms)
    r <- map$row[map$block == block][k]
    if (anyNA(r))
      stop("internal: no '", block, "' draw for ", paste(names[is.na(r)], collapse = ", "),
           call. = FALSE)
    r
  }
  eta <- m$X %*% D[by_name("beta", colnames(m$X)), , drop = FALSE]
  for (sm in names(m$smooth)) {
    rows <- map$row[map$block == "bvec" & map$term == sm]
    eta <- eta + m$smooth[[sm]] %*% D[rows, , drop = FALSE]
  }
  new_level <- rep(FALSE, n); new_unit <- rep(FALSE, n); drawn <- rep(FALSE, n)
  v_int <- matrix(0, n, S)                       # variance to integrate (population mean)
  z_draws <- 0                                   # counter for the seed of fresh effects
  for (t in setdiff(names(spec), "ar")) {
    r <- m$re[[t]]; Z <- as.matrix(r$Z); ch <- spec[[t]]
    known <- !r$new_group
    if (ch == "fitted") {
      ## a random intercept's design can come without a column name
      dims <- colnames(Z)
      if (is.null(dims)) dims <- unique(map$dim[map$block == "bvec" & map$term == t])
      if (length(dims) != ncol(Z))
        stop("internal: the design of '", t, "' does not match its draws", call. = FALSE)
      for (j in seq_len(ncol(Z))) {
        k <- match(paste("bvec", t, as.character(r$level), dims[j], sep = "\r"), key)
        if (any(known & is.na(k)))
          stop("internal: no draw for a fitted level of '", t, "'", call. = FALSE)
        if (any(known))
          eta[known, ] <- eta[known, ] + Z[known, j] * D[map$row[k[known]], , drop = FALSE]
      }
      new_level <- new_level | !known
    }
    fresh <- if (ch == "fitted") !known else if (ch == "typical") rep(FALSE, n)
      else rep(TRUE, n)
    if (!any(fresh)) next
    Sig <- dr$natural$re[[t]]
    if (ch == "population" && target == "mean") {
      for (s in seq_len(S)) v_int[fresh, s] <- v_int[fresh, s] +
          rowSums((Z[fresh, , drop = FALSE] %*% Sig[, , s]) * Z[fresh, , drop = FALSE])
      next
    }
    if (ch == "population" && target == "link") next   # the average of eta is eta
    new_unit <- new_unit | fresh
    z_draws <- z_draws + 1L
    eta[fresh, ] <- eta[fresh, ] + lmn_fresh_effects(Z[fresh, , drop = FALSE],
                                                     as.character(r$level)[fresh], Sig,
                                                     lmn_sub_seed(seed, 1000L * z_draws))
  }
  if ("ar" %in% names(spec)) {
    res <- lmn_ar_contribution(fit, m, dr, spec[["ar"]], target, lmn_sub_seed(seed, 7919L), noise)
    eta <- eta + res$eta
    v_int <- v_int + res$v_int
    new_level <- new_level | res$new_level
    new_unit <- new_unit | res$new_unit
    drawn <- res$drawn
  }
  lp <- list(mu = eta)
  if ("zi" %in% dist$lp_names)
    lp$zi <- m$zi %*% D[by_name("gzi", colnames(m$zi), "zi:"), , drop = FALSE]
  if ("disp" %in% dist$lp_names)
    lp$disp <- m$disp %*% D[by_name("gamma", colnames(m$disp), "disp:"), , drop = FALSE]
  par <- list(); pvec <- list()
  for (p in dist$par_names) {
    rows <- map$row[map$block == p]
    if (length(rows) == 1L) par[[p]] <- D[rows, ] else pvec[[p]] <- D[rows, , drop = FALSE]
  }
  list(lp = lp, par = par, pvec = pvec, v_int = v_int, new_level = new_level,
       new_unit = new_unit, drawn = drawn)
}

## A new unit's effects: in each draw, one draw from N(0, Sigma_s) per unit,
## shared by that unit's rows, times the rows' designs.
#' @keywords internal
#' @noRd
lmn_fresh_effects <- function(Z, unit, Sig, seed) {
  S <- dim(Sig)[3]; d <- ncol(Z)
  u_of <- match(unit, unique(unit)); U <- length(unique(unit))
  lmn_with_seed(seed, {
    out <- matrix(0, nrow(Z), S)
    for (s in seq_len(S)) {
      e <- eigen(matrix(Sig[, , s], d, d), symmetric = TRUE)
      root <- e$vectors %*% diag(sqrt(pmax(e$values, 0)), d)   # boundaries allowed
      u <- matrix(stats::rnorm(U * d), U, d) %*% t(root)
      out[, s] <- rowSums(Z * u[u_of, , drop = FALSE])
    }
    out
  })
}

## The correlation over time's part of each row (see R/latent.R). Under
## "fitted", a row at a fitted cell takes that cell's draw (zero at a random
## walk's anchor), and any other time of a seen series continues its cells: a
## forecast, a bridge or a backcast. A new series, or any under "new", and any
## under "population" for a new observation, is a fresh path. "typical" is zero,
## and a population mean integrates over the latent's variance at each time.
#' @keywords internal
#' @noRd
lmn_ar_contribution <- function(fit, m, dr, ch, target, seed, noise = TRUE) {
  ar <- m$ar; D <- dr$draws; map <- dr$map; S <- ncol(D); n <- nrow(ar)
  out <- list(eta = matrix(0, n, S), v_int = matrix(0, n, S),
              new_level = if (ch == "fitted") ar$new_group else rep(FALSE, n),
              new_unit = rep(FALSE, n), drawn = rep(FALSE, n))
  if (ch == "typical" || (ch == "population" && target == "link")) return(out)
  cells <- illume::ilm_cells(fit)
  par <- lmn_latent_par(dr$natural$latent, S, attr(cells, "step"))
  if (ch == "population" && target == "mean") {
    out$v_int <- lmn_latent_var(ar, cells, par)
    return(out)
  }
  a_cell <- matrix(0, nrow(cells), S)
  br <- map$block == "B_ar"
  a_cell[map$cell[br], ] <- D[map$row[br], , drop = FALSE]
  fresh <- if (ch == "fitted") ar$new_group else rep(TRUE, n)
  lp <- lmn_latent_paths(ar, cells, a_cell, par, fresh, seed, noise)
  out$eta <- lp$values
  out$new_unit <- fresh
  out$drawn <- lp$drawn
  out
}

## Each draw's mean of the response. For a population mean the latent part is
## integrated by illume's own quadrature, as predict(groups = "population") does.
#' @keywords internal
#' @noRd
lmn_mean_values <- function(lpd, dist) {
  lp <- lpd$lp; n <- nrow(lp$mu); S <- ncol(lp$mu)
  out <- matrix(NA_real_, n, S)
  for (s in seq_len(S)) {
    lps <- lapply(lp, function(M) M[, s])
    par <- lapply(lpd$par, function(v) v[s])
    for (p in names(lpd$pvec)) par[[p]] <- lpd$pvec[[p]][, s]
    sd <- sqrt(lpd$v_int[, s])
    if (any(sd > 0)) {
      out[, s] <- illume::ilm_normal_expect(function(e) {
        l2 <- lps; l2$mu <- e; dist$mean(l2, par) }, eta = lps$mu, sd = sd)
    } else out[, s] <- dist$mean(lps, par)
  }
  out
}

## ---- methods ------------------------------------------------------------------------

#' Point estimates and intervals from a predictive distribution
#'
#' @param x,object A predictive distribution from [ilm_predict()].
#' @param level Coverage of the central intervals.
#' @param row.names,optional,... Unused; for the generics.
#' @return A data frame: the rows predicted, the point estimate, the mean, and
#'   for each level the interval's lower and upper limits.
#'   - For a new observation, the estimate is the predictive median.
#'   - For a mean or a linear predictor, the estimate is the value at the
#'     estimates: for a seen series' new time, its fitted cells carried forward
#'     with no new shocks. For a gaussian model that is a Kalman filter's point
#'     forecast; limn's tests check it against `stats::KalmanForecast()`.
#'   - `mean` is the expected value, averaged over everything the prediction is
#'     uncertain about: the parameters, a new unit's effects, and a series'
#'     future shocks. On a log or logit scale it can sit well above the
#'     estimate, most for forecasts at long horizons and for new units. Use it
#'     for an expected count, for totals (expected values add up, medians do
#'     not) and for squared error; use the estimate for a central value. For a
#'     forecast, `mean` is the point forecast as fpp3 and its fable package give
#'     it: the mean of the forecast distribution, in fable's `.mean` column
#'     (sections 5.1 and 5.5).
#'   - An ordered outcome is summarised in its categories, with no mean.
#'
#'   `summary()` is `as.data.frame()` with the 50% and 90% intervals.
#' @references Hyndman, R. J. and Athanasopoulos, G. (2021). *Forecasting:
#'   Principles and Practice*, 3rd edition. OTexts. Section 5.1, the point
#'   forecast in fable's `.mean` column, and section 5.5, "Distributional
#'   forecasts and prediction intervals" (online edition,
#'   <https://otexts.com/fpp3/>, read 17 October 2026).
#' @seealso [ilm_predict()], which makes the prediction; [ilm_aggregate()] for
#'   totals.
#' @examples
#' set.seed(1)
#' d <- data.frame(site = factor(rep(1:12, each = 6)), x = rnorm(72))
#' d$y <- rpois(72, exp(0.5 + 0.3 * d$x + rnorm(12, 0, 0.4)[d$site]))
#' fit <- illume::ilm_model(y ~ x + (1 | site), data = d, family = "poisson",
#'                          verbose = FALSE)
#' p <- ilm_predict(fit, data.frame(site = c("1", "new"), x = 0),
#'                  target = "response", groups = "fitted", nsim = 200)
#' as.data.frame(p, level = c(0.5, 0.9))
#' summary(p)
#' @export
as.data.frame.ilm_pred <- function(x, row.names = NULL, optional = FALSE, ...,
                                   level = 0.9) {
  out <- x$rows
  ## an ordered outcome is summarised in its categories: the median category and
  ## the interval's end categories, with no numeric mean
  cats <- if (lmn_is_empirical(x)) NULL else x$dist$categories
  as_cat <- function(v) if (is.null(cats)) v else factor(cats[v], levels = cats, ordered = TRUE)
  if (lmn_is_empirical(x)) {
    out$estimate <- x$estimate
    out$mean <- lmn_pred_mean(x)
  } else {
    out$estimate <- as_cat(lmn_pred_quantile(x, 0.5))
    if (is.null(cats)) out$mean <- lmn_pred_mean(x)
  }
  for (L in level) {
    out[[paste0("lower_", L)]] <- as_cat(lmn_pred_quantile(x, (1 - L) / 2))
    out[[paste0("upper_", L)]] <- as_cat(lmn_pred_quantile(x, (1 + L) / 2))
  }
  out
}

#' @rdname as.data.frame.ilm_pred
#' @export
summary.ilm_pred <- function(object, level = c(0.5, 0.9), ...) {
  as.data.frame(object, level = level)
}
