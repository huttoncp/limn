#' Forecasts from a fitted illume model
#'
#' `ilm_predict()` with the settings a forecast wants: the distribution of a
#' new observation (`target = "response"`), for each series continued from
#' its own fitted state (`groups = "fitted"`). Every argument of
#' `ilm_predict()` can be given, and overrides these.
#'
#' With a correlation over time, each series is carried forward from its last
#' fitted cell, one path per draw, so a series' forecasts at successive times
#' are drawn jointly (see *Over time* in [ilm_predict()]). A series the fit
#' has not seen is a new series: a fresh path, marked in `new_level`. Build
#' the rows with [ilm_future()], and draw the fan chart with [ilm_plot_fan()].
#'
#' @section Checks:
#' A forecast makes [ilm_predict()]'s checks, and they are explained there
#' (*Checks* in [ilm_predict()]): few groups, units the fit has not seen, a
#' deterministic trend in time, and draws that are not finite or could not be
#' computed, each with what to do about it. The checks are warned, and kept on the forecast in
#' `checks`.
#'
#' @inheritParams ilm_predict
#' @param ... Other arguments to [ilm_predict()]: `uncertainty`, `set`, `nsim`,
#'   `seed`, `time` and `group`.
#' @return A predictive distribution, as from [ilm_predict()]; *The prediction
#'   object* there lists its elements.
#' @seealso [ilm_future()] for the rows, [ilm_covariate_paths()] for covariates
#'   not known in advance, [ilm_plot_fan()] for the fan chart, and
#'   [ilm_aggregate()] for totals over series or horizons.
#' @examples
#' set.seed(1)
#' d <- data.frame(id = factor(rep(c("a", "b", "c"), each = 30)), t = rep(1:30, 3))
#' d$y <- 10 + ave(rnorm(90, 0, 0.4), d$id, FUN = cumsum) + rnorm(90, 0, 0.5)
#' fit <- illume::ilm_model(y ~ 1, data = d, family = "gaussian",
#'                          ar = illume::ilm_rw1(~ t | id), verbose = FALSE)
#' fut <- ilm_future(d, time = "t", group = "id", h = 6, response = "y")
#' fc <- ilm_forecast(fit, fut, nsim = 500)
#' head(as.data.frame(fc))
#' ilm_plot_fan(fc, time = "t", group = "id", history = d, response = "y")
#' @export
ilm_forecast <- function(fit, newdata, target = c("response", "mean", "link"),
                         groups = "fitted", ...) {
  ilm_predict(fit, newdata, target = match.arg(target), groups = groups, ...)
}
