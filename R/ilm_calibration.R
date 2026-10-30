## The reliability table's "Wilson intervals" are Wilson's (1927) score
## intervals: awaiting the primary source. It has not been read; the paper is
## on the literature list: Wilson, E. B. (1927). Probable inference, the law of
## succession, and statistical inference. Journal of the American Statistical
## Association 22(158), 209-212. doi:10.1080/01621459.1927.10502953.

#' Calibration of a prediction against observed outcomes
#'
#' Whether outcomes fell where their predictive distributions said they
#' would: coverage and width level by level, bias in the mean, and the
#' uniformity of the probability integral transforms (PITs; Gneiting,
#' Balabdaoui and Raftery 2007, section 3.1, pp. 251-253). A binary outcome is
#' judged instead by its calibration intercept and slope, which together assess
#' what Van Calster et al. (2019, p. 4) call weak calibration. The result ends
#' with a verdict in plain words and, where something is wrong, a remedy that
#' exists, most often in `illume`.
#'
#' This is limn's method for `illume`'s [illume::ilm_calibration()] generic,
#' for a prediction. `illume`'s own method checks a fitted model's category
#' probabilities.
#'
#' @section What is checked:
#' - **Coverage and width.** At each level in `levels`, the share of outcomes
#'   inside the central interval of their predictive distributions, with its
#'   standard error, and the intervals' mean width. Coverage is read from the
#'   PITs, not from whether the outcome lies inside a quantile interval. A
#'   discrete outcome's quantile interval covers more than its level, because
#'   of ties, while the randomised PIT's coverage is exactly nominal when the
#'   model is right.
#' - **Bias**, on the outcome's own scale: the mean of the outcome minus its
#'   predictive mean, with its standard error. The mean PIT is not used, since
#'   a distribution with the right mean and the wrong shape moves it too.
#' - **PIT uniformity**: a chi-square test over `bins` equal bins, with at
#'   least 5 outcomes expected in each (limn's choice of minimum).
#' - **A binary outcome**: the calibration intercept (are the probabilities
#'   right on average?) and slope (are they too extreme, or too moderate?),
#'   with Wald 95% intervals, and a reliability table of equal-count bins with
#'   Wilson intervals. Interval coverage adds nothing for a 0/1 outcome.
#'
#' A coverage or bias departure is called when it is more than 2.58 standard
#' errors from its ideal (the two-sided 99% point of the normal distribution,
#' limn's own choice), and a
#' binary measure when its 95% interval excludes the ideal. The intercept's
#' ideal is 0 and the slope's is 1 (Van Calster et al. 2019, p. 4).
#'
#' A biased mean is reported first. It alone can make the intervals cover too
#' little, since they are centred in the wrong place. So a shortfall in coverage
#' beside it is reported as one the bias can cause, and the remedies that widen
#' intervals wait until the mean is fixed. A calibration slope whose interval
#' lies wholly below 0 is reported as reversed, a verdict of limn's own: the
#' predictions order the outcomes the wrong way, which is more than being too
#' extreme.
#'
#' @eval lmn_diag_rd(c("calib_bias", "calib_narrow_bias", "calib_narrow_counts", "calib_narrow_dispersion", "calib_narrow", "calib_wide", "calib_ok", "binary_intercept", "binary_reversed", "binary_overconfident", "binary_underconfident", "binary_ok"), "Verdicts")
#'
#' @section The rows are treated as independent:
#' Rows that share a group, a series or a forecast origin are not
#' independent, and then the standard errors and tests here are too
#' confident. For forecasts, check calibration horizon by horizon, or over
#' rows from different series and origins.
#'
#' @param object A prediction from [ilm_predict()] or [ilm_forecast()] with
#'   `target = "response"`.
#' @param y The observed outcomes, one per row of the prediction. At least 20
#'   must be scoreable.
#' @param bins The number of bins for the PIT's uniformity test, and for a
#'   binary outcome's reliability table. There are fewer when the outcomes are
#'   few.
#' @param levels The central interval levels at which coverage is checked.
#' @param seed The seed for the randomised PITs of a discrete outcome. The
#'   caller's random-number state is left as it was.
#' @param ... Unused.
#' @return An object of class `"ilm_pred_calibration"`, a list:
#'   - `pit`, each row's PIT;
#'   - `n`, the outcomes used;
#'   - `histogram` and `uniformity`;
#'   - `bias`;
#'   - `coverage`, a data frame of level, coverage, its standard error and
#'     mean width;
#'   - `family`;
#'   - for a binary outcome, `binary` (the intercept and slope) and
#'     `reliability`;
#'   - `checks`, the diagnostics as data;
#'   - `verdict`, their words as printed.
#'
#'   [ilm_plot_calibration()] draws it.
#' @references Gneiting, T., Balabdaoui, F. and Raftery, A. E. (2007).
#'   Probabilistic forecasts, calibration and sharpness. *Journal of the Royal
#'   Statistical Society B* 69, 243-268. Section 3.1, pp. 251-253: the PIT, its
#'   uniformity, coverage read from it, and the histogram's shapes.
#'
#'   Van Calster, B., McLernon, D. J., van Smeden, M., Wynants, L. and
#'   Steyerberg, E. W. (2019). Calibration: the Achilles heel of predictive
#'   analytics. *BMC Medicine* 17, 230. p. 2: causes of poor calibration; p. 4:
#'   weak calibration, the intercept and slope and their targets, and remedies.
#' @seealso [ilm_pit()], [ilm_plot_calibration()], [ilm_scores.ilm_pred()]
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rnorm(400))
#' d$y <- rnbinom(400, size = 1.5, mu = exp(1 + 0.5 * d$x))
#' train <- d[1:200, ]; test <- d[201:400, ]
#' ## Poisson on overdispersed counts: intervals too narrow, and the remedy
#' fp <- illume::ilm_model(y ~ x, data = train, family = "poisson", verbose = FALSE)
#' illume::ilm_calibration(ilm_predict(fp, test, target = "response", nsim = 200), test$y)
#' ## the negative binomial the verdict names
#' fn <- illume::ilm_model(y ~ x, data = train, family = "nbinom", verbose = FALSE)
#' illume::ilm_calibration(ilm_predict(fn, test, target = "response", nsim = 200), test$y)
#' @exportS3Method illume::ilm_calibration
ilm_calibration.ilm_pred <- function(object, y, bins = 10L, levels = c(0.5, 0.8, 0.9, 0.95),
                                     seed = 1L, ...) {
  if (!is.numeric(levels) || !length(levels) || anyNA(levels) || any(levels <= 0 | levels >= 1))
    stop("`levels` must be interval levels between 0 and 1.", call. = FALSE)
  lmn_pred_calibration(object, y, bins = bins, levels = levels, seed = seed)
}
