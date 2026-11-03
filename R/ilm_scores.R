## Hyndman and Koehler (2006), cited below for the scaled errors: awaiting the
## primary source. It has not been read; the paper is on the literature list.
## fpp3 section 5.8 is cited for the formulas meanwhile.

#' Proper scores of a prediction against observed outcomes
#'
#' Scores each row of a prediction against its observed outcome, by proper
#' scoring rules computed from the predictive distribution itself, and by the
#' errors of the point forecasts it gives. Every score is a loss: lower is
#' better.
#'
#' @section The scoring rules:
#' For a row with outcome `y`, predictive CDF `F`, density or probability mass
#' `f`, mean `m` and median `med`, written as R expressions (a comparison such
#' as `(y < q)` counts as 1 when true and 0 when false):
#' * `"log"`: `-log(f(y))`, the negative of the logarithmic score (Gneiting and
#'   Raftery 2007, section 4.1, equation (19), p. 365; for a discrete outcome,
#'   the probability mass, section 3.1, Example 3, p. 363).
#' * `"crps"`: the integral over `z` of `(F(z) - (z >= y))^2`, the continuous
#'   ranked probability score in its negative orientation (section 4.2,
#'   equation (20), p. 367). For a count or an ordered outcome, `F` steps only
#'   at whole numbers, so the integral is the sum over them of
#'   `(F(k) - (k >= y))^2`.
#' * `"brier"`: `(P(Y = 1) - y)^2`, for a binary outcome coded 0/1: the Brier
#'   score for a dichotomous event (Table 1, p. 365), as a loss.
#' * `"pinball"`: `(y - q) * (p - (y < q))`, with `q` the predictive
#'   `p`-quantile, for each `p` in `probs`: the negative of the quantile score
#'   (section 6.1, equation (41), p. 370), the "tick" or "check" loss.
#' * `"interval"`: `(u - l) + (2 / alpha) * (l - y) * (y < l) + (2 / alpha) *
#'   (y - u) * (y > u)`, for the central interval from `l` to `u` at `level` =
#'   1 - `alpha`: the interval score (section 6.2, equation (43), p. 370).
#' * `"squared"`: `(y - m)^2`. The energy score with exponent 2 reduces to
#'   this, in the opposite orientation (section 4.3, equation (23), p. 367), a
#'   proper but not strictly proper score.
#' * `"absolute"`: `abs(y - med)`, twice the pinball loss at `p` = 1/2.
#' * `"mase"`: `abs(y - med) / s`, with `s` the series' in-sample mean absolute
#'   naive error at lag `period`. The mean over rows is the mean absolute
#'   scaled error.
#' * `"rmsse"`: `(y - m)^2 / s2`, with `s2` the series' in-sample mean squared
#'   naive error. The square root of the mean over rows is the root mean
#'   squared scaled error.
#'
#' The scaled errors and their in-sample naive scales are as fpp3 defines them
#' (section 5.8, "Scaled errors"), where they are credited to Hyndman and
#' Koehler (2006). fpp3 scales the error of a point forecast; `limn` takes the
#' median as the point forecast for the absolute error and the mean for the
#' squared error, the point forecasts each error rewards (its own choice).
#'
#' @section Giving MASE and RMSSE their scale:
#' The scaled errors (`mase`, `rmsse`) divide each row's error by its series'
#' in-sample naive-forecast error, so they need each series' history. There are
#' two ways to give it, and they give identical scores wherever both apply:
#'
#' - **`history = train`**, the common case: the data the model was trained
#'   on. The prediction already knows its time, series and response columns
#'   (it records them), so nothing else is needed. `period` sets the naive
#'   forecast's lag: 1, or the season's length for the seasonal naive.
#' - **`naive = ilm_naive_scale(train, time =, group =, response =)`**: the
#'   scale computed once and passed in. Use it when the prediction doesn't know
#'   its time column (see *Time and series columns* in [ilm_predict()]), or when several
#'   models' predictions must share one scale, so their scaled errors compare
#'   on the same footing.
#'
#' When `history =` can't work, because the prediction has no time column, the
#' error says so and gives the `naive =` call to use instead. Give one of the
#' two, not both.
#'
#' @param object A prediction from [ilm_predict()] or [ilm_forecast()] with
#'   `target = "response"`.
#' @param y The observed outcomes, one per row of the prediction.
#' @param rules The scoring rules: any of `"log"`, `"crps"`, `"brier"`,
#'   `"pinball"`, `"interval"`, `"squared"`, `"absolute"`, `"mase"` and
#'   `"rmsse"`.
#' @param probs The probabilities of the pinball loss.
#' @param level The coverage of the interval score's central interval.
#' @param naive For `mase` and `rmsse`: each series' scale, from
#'   [ilm_naive_scale()]. See *Giving MASE and RMSSE their scale*.
#' @param history For `mase` and `rmsse`: the data the model was trained on,
#'   from which each series' scale is computed. See *Giving MASE and RMSSE
#'   their scale*.
#' @param period With `history`, the naive forecast's lag: 1 for the naive
#'   forecast, the season's length for the seasonal naive.
#' @param ... Unused.
#' @return A data frame: the prediction's rows, with one column per score:
#'   `log`, `crps`, `brier`, `pinball_<p>` for each p, `interval_<level>`,
#'   `squared`, `absolute`, `mase` and `rmsse`, as asked for. Attribute `rules`
#'   names the score columns, and `orientation` is "lower is better". For a
#'   discrete outcome whose CRPS was summed approximately, attributes
#'   `crps_bound` and `crps_bracket` give each row's error estimate and its
#'   rigorous bound.
#' @references Gneiting, T. and Raftery, A. E. (2007). Strictly proper scoring
#'   rules, prediction, and estimation. *Journal of the American Statistical
#'   Association* 102(477), 359-378. Section 3.1, Example 3 (p. 363); Table 1
#'   (p. 365); section 4.1, equation (19) (p. 365); section 4.2, equation (20)
#'   and section 4.3, equation (23) (p. 367); section 6.1, equation (41), and
#'   section 6.2, equation (43) (p. 370).
#'
#'   Hyndman, R. J. and Athanasopoulos, G. (2021). *Forecasting: Principles
#'   and Practice*, 3rd edition. OTexts. Section 5.8, "Evaluating point forecast
#'   accuracy" (online edition,
#'   <https://otexts.com/fpp3/toolbox.html#sec-accuracy>, read 27 October 2026).
#'
#'   Hyndman, R. J. and Koehler, A. B. (2006). Another look at measures of
#'   forecast accuracy. *International Journal of Forecasting* 22, 679-688.
#' @seealso [ilm_naive_scale()] for the scaled errors' scale; [ilm_cv()] for
#'   scores out of sample; [ilm_compare()] to compare two models' scores;
#'   [ilm_pit()] and [ilm_calibration()][ilm_calibration.ilm_pred] for
#'   calibration.
#' @examples
#' set.seed(1)
#' d <- data.frame(site = factor(rep(c("a", "b", "c"), each = 16)), t = rep(1:16, 3))
#' d$y <- rpois(nrow(d), exp(1.5 + ave(rnorm(nrow(d), 0, 0.15), d$site, FUN = cumsum)))
#' train <- d[d$t <= 12, ]; test <- d[d$t > 12, ]
#' fit <- illume::ilm_model(y ~ 1, data = train, family = "poisson",
#'                          ar = illume::ilm_rw1(~ t | site), verbose = FALSE)
#' fc <- ilm_forecast(fit, test[c("site", "t")], nsim = 200)
#' ## the scale read from the training data, through the prediction's own columns
#' s <- ilm_scores(fc, test$y, rules = c("crps", "mase"), history = train)
#' colMeans(s[attr(s, "rules")])
#' ## the same scale computed once, as several models would share it
#' sc <- ilm_naive_scale(train, time = "t", group = "site", response = "y")
#' all.equal(s$mase, ilm_scores(fc, test$y, rules = "mase", naive = sc)$mase)
#' @exportS3Method illume::ilm_scores
ilm_scores.ilm_pred <- function(object, y, rules = c("log", "crps"), probs = c(0.1, 0.5, 0.9),
                                level = 0.9, naive = NULL, history = NULL, period = 1L, ...) {
  if (!is.null(history)) {
    if (!is.null(naive))
      stop("give the scaled errors' scale one way: `naive` or `history`, not both.", call. = FALSE)
    naive <- lmn_naive_from_history(object, history, period)
  }
  lmn_pred_scores(object, y, rules = rules, probs = probs, level = level, naive = naive)
}

## The scale from the training data, through the columns the prediction
## recorded; when it has no time column, an error that gives the naive = call.
#' @keywords internal
#' @noRd
lmn_naive_from_history <- function(pred, history, period) {
  if (!is.data.frame(history))
    stop("`history` must be the data frame the model was trained on.", call. = FALSE)
  m <- pred$meta
  if (is.null(m$time))
    stop("`history =` needs the prediction's time column, and this prediction has none (neither ",
         "the model nor the rows said which it is; ilm_predict() and ilm_forecast() take it as ",
         "`time =`), so the series' order in the history is unknown. Give the scale directly ",
         "instead:\n  ilm_scores(pred, y, rules = \"mase\", naive = ilm_naive_scale(",
         "history, time = \"<time column>\", group = \"<series column>\", response = \"",
         m$response %||% "<response column>", "\"))", call. = FALSE)
  lmn_naive_scale(history, time = m$time, group = m$series, response = m$response, period = period)
}

#' Each series' in-sample naive-forecast error, the scale of MASE and RMSSE
#'
#' The mean absolute and mean squared error of the naive forecast (the value
#' `period` steps back) over each series' history: the scale that makes MASE
#' and RMSSE comparable across series. These are the denominators of fpp3's
#' scaled errors, the mean over t = `period` + 1 to T of |y_t - y_(t -
#' `period`)| and of its square (section 5.8, "Scaled errors"), which fpp3
#' credits to Hyndman and Koehler (2006). Compute it from the data the model was
#' trained on, and pass it to [ilm_scores.ilm_pred()] as `naive`.
#'
#' @inheritSection ilm_scores.ilm_pred Giving MASE and RMSSE their scale
#' @param history A data frame: the observed history, as the model was trained
#'   on it.
#' @param time The name of its time column.
#' @param group The name of the column identifying each series, or `NULL` for a
#'   single series. The prediction's rows need the same column.
#' @param response The name of the outcome column.
#' @param period The lag of the naive forecast: 1 for the naive forecast, the
#'   season's length for the seasonal naive.
#' @return A data frame with one row per series: `key`, and the naive error's
#'   mean absolute value `abs` and mean square `sq`.
#' @references Hyndman, R. J. and Athanasopoulos, G. (2021). *Forecasting:
#'   Principles and Practice*, 3rd edition. OTexts. Section 5.8, "Evaluating
#'   point forecast accuracy" (online edition,
#'   <https://otexts.com/fpp3/toolbox.html#sec-accuracy>, read 27 October 2026).
#'
#'   Hyndman, R. J. and Koehler, A. B. (2006). Another look at measures of
#'   forecast accuracy. *International Journal of Forecasting* 22, 679-688.
#' @seealso [ilm_scores.ilm_pred()], which uses it for `mase` and `rmsse`.
#' @examples
#' d <- data.frame(id = rep(c("a", "b"), each = 5), t = rep(1:5, 2),
#'                 y = c(1, 3, 2, 5, 4, 10, 12, 11, 15, 14))
#' ilm_naive_scale(d, time = "t", group = "id", response = "y")
#' @export
ilm_naive_scale <- function(history, time, group = NULL, response, period = 1L) {
  lmn_naive_scale(history, time, group = group, response = response, period = period)
}
