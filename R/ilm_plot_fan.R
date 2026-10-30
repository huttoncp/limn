#' Fan chart of predictive distributions over time
#'
#' Draws each row's predictive median as a line, with nested bands for the
#' central intervals at `levels`, against time and optionally one panel per
#' series. The observed history can be drawn before them. The bands are read
#' from each row's own predictive distribution, so for a forecast they widen
#' with the horizon as far as the model says the future is uncertain.
#'
#' @section Reading the plot:
#' Time runs along the horizontal axis, with one panel per series when `group`
#' is given.
#' * The grey points and line are the observed history.
#' * The dark line is each row's predictive median: half the predicted outcomes
#'   fall above it, half below.
#' * The shaded bands are the central intervals, darkest for the narrowest (by
#'   default 50%, then 80% and 95%). If the model is right, about that share of
#'   future outcomes fall inside each band.
#' * Bands that widen with the horizon show the model's uncertainty about the
#'   future growing. A band that is too narrow or too wide is something only
#'   outcomes can show: score the forecasts against them with
#'   [ilm_scores.ilm_pred()], or check their coverage with
#'   [ilm_calibration()][ilm_calibration.ilm_pred].
#'
#' @param pred A predictive distribution, class `"ilm_pred"`, whose rows have a
#'   time column (and a series column, for panels).
#' @param time The name of the time column, in the prediction's rows and in
#'   `history`.
#' @param group The name of the series column, for one panel per series, or
#'   `NULL` for a single panel.
#' @param history Optional data frame of observed values to draw before the
#'   predictions, with the `time`, `group` and `response` columns.
#' @param response The name of the response column in `history`.
#' @param levels Coverages of the central intervals drawn as bands.
#' @param ... Passed to [tinyplot::tinyplot()]: for example
#'   `facet.args = list(free = TRUE)` to give each series' panel its own y
#'   scale, or `ylim` to set the scale.
#' @return The data behind the plot, invisibly: one row per history row and
#'   predicted row, with the observed value, the median, and each band's
#'   limits.
#' @seealso [ilm_forecast()], which makes the forecast; [ilm_future()] for the
#'   rows a forecast needs.
#' @examples
#' set.seed(1)
#' d <- data.frame(id = factor(rep(c("a", "b"), each = 30)), t = rep(1:30, 2))
#' d$y <- 10 + ave(rnorm(60, 0, 0.4), d$id, FUN = cumsum) + rnorm(60, 0, 0.5)
#' fit <- illume::ilm_model(y ~ 1, data = d, family = "gaussian",
#'                          ar = illume::ilm_rw1(~ t | id), verbose = FALSE)
#' fc <- ilm_forecast(fit, ilm_future(d, time = "t", group = "id", h = 6,
#'                                    response = "y"), nsim = 300)
#' ilm_plot_fan(fc, time = "t", group = "id", history = d, response = "y")
#' @export
ilm_plot_fan <- function(pred, time, group = NULL, history = NULL, response = NULL,
                         levels = c(0.5, 0.8, 0.95), ...) {
  if (!inherits(pred, "ilm_pred"))
    stop("`pred` must be a predictive distribution (class \"ilm_pred\").",
         call. = FALSE)
  rows <- pred$rows
  lmn_check_column(rows, time, "time")
  if (!is.null(group)) lmn_check_column(rows, group, "group")
  if (!is.numeric(levels) || any(levels <= 0 | levels >= 1))
    stop("`levels` must be coverages strictly between 0 and 1.", call. = FALSE)
  levels <- sort(unique(levels), decreasing = TRUE)      # widest band first

  fc <- data.frame(time = rows[[time]],
                   group = if (is.null(group)) "all" else as.character(rows[[group]]),
                   observed = NA_real_, median = lmn_pred_quantile(pred, 0.5))
  for (L in levels) {
    fc[[paste0("lower_", L)]] <- lmn_pred_quantile(pred, (1 - L) / 2)
    fc[[paste0("upper_", L)]] <- lmn_pred_quantile(pred, (1 + L) / 2)
  }
  d <- fc
  if (!is.null(history)) {
    if (is.null(response))
      stop("give `response`, the history's column of observed values.", call. = FALSE)
    lmn_check_column(history, time, "time")
    y <- lmn_check_column(history, response, "response")
    hs <- fc[rep(NA_integer_, nrow(history)), , drop = FALSE]
    hs$time <- history[[time]]
    hs$group <- if (is.null(group)) "all" else
      as.character(lmn_check_column(history, group, "group"))
    hs$observed <- as.numeric(y)
    d <- rbind(hs, fc)
  }
  d <- d[order(d$group, d$time), , drop = FALSE]
  rownames(d) <- NULL

  facet <- if (is.null(group)) NULL else d$group
  alpha <- seq(0.18, 0.45, length.out = length(levels))
  band <- function(L) list(d[[paste0("lower_", L)]], d[[paste0("upper_", L)]])
  ## an empty first layer over every value drawn (history, medians, every band's
  ## limits), so that each panel's scale, shared or free, holds all of them
  vals <- d[setdiff(names(d), c("time", "group"))]
  all_y <- unlist(vals, use.names = FALSE)
  keep <- !is.na(all_y)
  tinyplot::tinyplot(x = rep(d$time, ncol(vals))[keep], y = all_y[keep],
                     facet = if (is.null(facet)) NULL else rep(facet, ncol(vals))[keep],
                     ymin = NULL, ymax = NULL,   # named, so the layers below can set or drop them
                     type = "n", xlab = time,
                     ylab = if (is.null(response)) "prediction" else response, ...)
  for (k in seq_along(levels)) {
    bk <- band(levels[k])
    tinyplot::tinyplot_add(x = d$time, y = NULL, ymin = bk[[1]], ymax = bk[[2]], facet = facet,
                           type = tinyplot::type_ribbon(alpha = alpha[k]), col = "#2a6f97")
  }
  tinyplot::tinyplot_add(x = d$time, y = d$median, ymin = NULL, ymax = NULL, facet = facet,
                         type = "l", col = "#1b4965", lwd = 2)
  if (!is.null(history))
    tinyplot::tinyplot_add(x = d$time, y = d$observed, ymin = NULL, ymax = NULL, facet = facet,
                           type = "o", pch = 16, cex = 0.7, col = "grey20")
  invisible(d)
}
