#' limn: A Projection Engine for Distributional Forecasting
#'
#' Predictive distributions, forecasts, proper scoring rules, calibration,
#' cross-validation and model comparison for models fitted by
#' [illume::ilm_model()]. `limn` is the third package in a family, after
#' `illumex` for exploration and `illume` for inference, and shares their
#' `ilm_` prefix, so one prefix finds a function in any of the three.
#'
#' @section A prediction is a distribution:
#' Every prediction `limn` makes is a predictive distribution: a mixture, over
#' draws of the fitted parameters, of the model family's distribution for the
#' response. Intervals, quantiles and scores are all read from that mixture.
#' Scoring it directly, rather than scoring responses simulated from it, removes
#' the simulation's noise from the score.
#'
#' @section Scores are losses:
#' Every score `limn` reports is a loss, so lower is better: the log score is
#' the negative log predictive density, and the CRPS, Brier, pinball, interval,
#' squared and absolute scores are penalties.
#'
#' @section The workflow:
#' * **Predict.** [ilm_predict()] gives the distribution of a mean, a new
#'   observation or the linear predictor, for a unit the fit has seen, a new
#'   one, a typical one or the average over units. [ilm_forecast()] carries
#'   series forward in time. [ilm_future()] builds the rows a forecast needs,
#'   and [ilm_covariate_paths()] and [ilm_as_draws()] give covariates that are
#'   not known in advance as draws.
#' * **Combine.** [ilm_aggregate()] totals a prediction over rows, and
#'   [ilm_average()] takes weighted averages of means (standardisation and
#'   post-stratification).
#' * **Check against outcomes.** [ilm_scores()][ilm_scores.ilm_pred] scores
#'   predictions by proper scoring rules; [ilm_pit()] and [ilm_plot_pit()] show
#'   where outcomes fell in their distributions; and
#'   [ilm_calibration()][ilm_calibration.ilm_pred] gives a verdict on
#'   calibration, with [ilm_plot_calibration()].
#' * **Validate out of sample.** [ilm_folds()] builds folds that follow the
#'   data's structure, and [ilm_cv()] scores a model on them.
#' * **Compare models.** [ilm_compare()] compares scores with standard errors
#'   that respect how the data hang together; [ilm_plot_compare()] and
#'   [ilm_plot_horizon()] draw the comparison.
#' * **Draw.** [ilm_plot_fan()] draws forecasts as fan charts.
#' * **Extend.** [ilm_predict_draws()] and [ilm_cv_refit()] let a package that
#'   fits illume's models another way supply the draws and the refits, so that
#'   every function above works on its fits.
#'
#' @section Checks, and what to do about them:
#' Every check `limn` makes takes its words from one registry. Each says what
#' was found, what it means, and what to do about it, and gives the source it
#' rests on, or says that the reasoning is `limn`'s own. A check is warned when
#' it is raised, with class `lmn_<id>`, and kept on the result in `checks`, so
#' it can be read as data as well as words. Where the fix is a change to the
#' model, the check names the call: `illume::ilm_remedies()` lists the changes
#' the fit's own checks support, and `illume::ilm_apply_remedy()` refits the
#' model with the one you choose.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(site = factor(rep(1:20, each = 8)), x = rnorm(160))
#' d$y <- rpois(160, exp(0.5 + 0.4 * d$x + rnorm(20, 0, 0.4)[d$site]))
#' train <- d[as.integer(d$site) <= 15, ]; test <- d[as.integer(d$site) > 15, ]
#' fit <- illume::ilm_model(y ~ x + (1 | site), data = train, family = "poisson",
#'                          verbose = FALSE)
#' ## the held-out sites are new units: predict their outcomes, then score
#' ## and check the predictions against what happened
#' p <- ilm_predict(fit, test, target = "response", groups = "new", nsim = 300)
#' head(as.data.frame(p))
#' s <- illume::ilm_scores(p, test$y)
#' colMeans(s[attr(s, "rules")])
#' illume::ilm_calibration(p, test$y)
#'
#' @keywords internal
#' @importFrom illume ilm_model
"_PACKAGE"
