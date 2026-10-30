## Ferro (2014), cited below for the samples' fair CRPS: awaiting the primary
## source. It has not been read; the paper is on the literature list.

#' Totals and averages of predictions over rows
#'
#' Sums or averages a prediction over groups of its rows: a series' forecasts
#' over a horizon, every series at each time, or all rows. The sum is taken
#' draw by draw. So whatever the rows share within a draw (the parameters, a
#' unit's effects, a series' path over time) carries into the total, and totals
#' at different levels add up.
#'
#' @section How a total is computed:
#' * **A prediction of the mean:** each draw's means are summed, and the
#'   estimate is the sum of the rows' estimates.
#' * **New observations from a Poisson or gaussian model with no zero part:**
#'   within a draw, a sum of such outcomes is again Poisson (with the summed
#'   mean) or gaussian (with the summed mean and variance). So the total is an
#'   exact mixture, and it is scored like any other prediction.
#' * **New observations from any other family:** outcomes are simulated from
#'   each draw, four per draw, and summed. The total is then held as those
#'   samples:
#'   - its CDF, quantiles, mean and PIT are the samples';
#'   - its CRPS is the fair estimator from them (Ferro 2014);
#'   - it has no log score, since samples have no density.
#'
#'   A fixed internal seed means totals over different groupings of the same
#'   prediction add up.
#' * **The linear predictor** has no total.
#'
#' The method used is printed, and kept in `meta$aggregate`. An average
#' (`fun = "mean"`) is the total over the number of rows. Where the family does
#' not stay closed under that division, it is scored through the total, so an
#' average of counts keeps its whole-number steps. Rows of different scenarios
#' are never added together.
#'
#' @eval lmn_diag_rd("missing_totals")
#' @param pred A prediction from [ilm_predict()] or [ilm_forecast()].
#' @param by Names of columns in the prediction's rows, as `as.data.frame()`
#'   shows them. There is one total for each combination of their values.
#'   `NULL`, the default, totals every row.
#' @param fun `"sum"`, the default, or `"mean"`.
#' @return A prediction for the totals (see *The prediction object* in
#'   [ilm_predict()]), with one row per combination of `by` (and of scenario),
#'   and `n_rows`, the number of rows each one combines. Its `what` and
#'   `meta$aggregate` say how the total was taken.
#' @references Ferro, C. A. T. (2014). Fair scores for ensemble forecasts.
#'   *Quarterly Journal of the Royal Meteorological Society* 140, 1917-1923.
#' @seealso [ilm_forecast()], [ilm_future()]; [ilm_average()] for weighted
#'   averages of means.
#' @examples
#' set.seed(1)
#' d <- data.frame(clinic = factor(rep(sprintf("c%02d", 1:20), each = 24)),
#'                 month = rep(1:24, 20))
#' d$visits <- rpois(nrow(d), exp(1.5 + ave(rnorm(nrow(d), 0, 0.15), d$clinic,
#'                                          FUN = cumsum)))
#' fit <- illume::ilm_model(visits ~ 1, data = d, family = "poisson",
#'                          ar = illume::ilm_rw1(~ month | clinic), verbose = FALSE)
#' fut <- ilm_future(d, time = "month", group = "clinic", h = 3, response = "visits")
#' fc <- ilm_forecast(fit, fut, nsim = 500)
#' ## each clinic's total over the quarter, and every clinic's total by month
#' head(as.data.frame(ilm_aggregate(fc, by = "clinic")))
#' as.data.frame(ilm_aggregate(fc, by = "month"))
#' @export
ilm_aggregate <- function(pred, by = NULL, fun = c("sum", "mean")) {
  if (!inherits(pred, "ilm_pred"))
    stop("`pred` must be a prediction from ilm_predict() or ilm_forecast(); it is ",
         class(pred)[1], ".", call. = FALSE)
  fun <- match.arg(fun)
  g <- lmn_agg_groups(pred$rows, by)
  out <- lmn_aggregate(pred, g$groups, g$out, fun = fun)
  out$checks <- pred$checks                      # a total keeps the prediction's checks
  out$origin <- lmn_agg_origin(pred$origin, g$groups, nrow(g$out))
  out
}

## A total's origin: the one its rows were all forecast from, or NA when they
## were forecast from different times, or not forecast at all.
#' @keywords internal
#' @noRd
lmn_agg_origin <- function(origin, groups, K) {
  if (is.null(origin)) return(rep(NA_real_, K))
  vapply(split(origin, factor(groups, levels = seq_len(K))), function(o)
    if (length(o) && !anyNA(o) && all(o == o[1])) o[1] else NA_real_, 0, USE.NAMES = FALSE)
}

## The groups of rows a total or an average combines: one per combination of
## `by`'s values, and never across scenarios. `out` has one row per group, with
## `n_rows`.
#' @keywords internal
#' @noRd
lmn_agg_groups <- function(rows, by) {
  if (!is.null(by)) {
    if (!is.character(by) || !length(by))
      stop("`by` must name columns of the prediction's rows.", call. = FALSE)
    miss <- setdiff(by, names(rows))
    if (length(miss))
      stop("`by` names column(s) the prediction's rows do not have: ",
           paste(miss, collapse = ", "), ". They have: ", paste(names(rows), collapse = ", "),
           ".", call. = FALSE)
  }
  ## never across scenarios
  keys <- unique(c(if (".scenario" %in% names(rows)) ".scenario", by))
  key <- if (length(keys)) do.call(paste, c(unname(as.list(rows[keys])), sep = "\r"))
    else rep("", nrow(rows))
  groups <- match(key, unique(key))
  out <- if (length(keys)) rows[!duplicated(groups), keys, drop = FALSE]
    else data.frame(row.names = 1L)
  rownames(out) <- NULL
  out$n_rows <- tabulate(groups, nrow(out))
  list(groups = groups, out = out)
}
