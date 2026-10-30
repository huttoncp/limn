## ---- diagnostic plots -------------------------------------------------------------
##
## Each draws with tinyplot, as illumex does, and returns the data behind the
## plot invisibly, so the numbers can be read or tested without the picture.

#' Scores by forecast horizon
#'
#' Plots each model's mean score at each forecast horizon, with confidence
#' intervals, as in the reference's plot of forecast accuracy against the
#' horizon (section 5.10, Figure 5.25). Forecasts usually get worse with the
#' horizon. The question is how fast, and whether one model's advantage holds
#' further out.
#'
#' @section Reading the plot:
#' The horizontal axis is the forecast horizon; the vertical axis is the mean
#' score there, and lower is better. Each model's points are joined by a line.
#' * A line that rises is a model whose forecasts get worse further ahead, as
#'   most do. A steeper rise than another model's means its advantage, if any,
#'   fades with the horizon.
#' * Where two models' intervals at a horizon do not overlap, their scores
#'   there differ. Where they overlap, the plot cannot say: [ilm_compare()]
#'   tests the difference directly, with a standard error of the difference
#'   itself, which allows for the two models' scores on the same rows being
#'   correlated.
#'
#' @section Standard errors:
#' Rows forecast from the same origin share that origin's surprises, so when the
#' scores have an `origin` column, each horizon's standard error is clustered by
#' origin. This is the same cluster-robust estimator [ilm_compare()] uses.
#' Otherwise the rows are treated as independent.
#'
#' @param ... One or more tables of scores with a `horizon` column, such as
#'   rolling-origin cross-validation scores. Name them to label the models.
#' @param rule The score to plot; by default the first score in the first
#'   table, as in [ilm_compare()].
#' @param level Confidence level of the intervals.
#' @return The data behind the plot, invisibly: one row per model and horizon,
#'   with the mean score, its standard error and interval, and the numbers of
#'   rows and origins.
#' @references Hyndman, R. J. and Athanasopoulos, G. (2021). *Forecasting:
#'   Principles and Practice*, 3rd edition. OTexts. Section 5.10, "Time series
#'   cross-validation", Figure 5.25 (online edition,
#'   <https://otexts.com/fpp3/toolbox.html#sec-tscv>, read 27 October 2026).
#' @seealso [ilm_compare()] to test whether two models' scores differ;
#'   [ilm_folds()] for rolling-origin folds; [ilm_cv()] for the scores.
#' @examples
#' set.seed(1)
#' s <- expand.grid(horizon = 1:4, origin = 1:12, site = 1:5)
#' ar <- data.frame(s, log = 1 + 0.08 * s$horizon + rnorm(nrow(s), 0, 0.3))
#' iid <- data.frame(s, log = ar$log + 0.03 * s$horizon + rnorm(nrow(s), 0, 0.1))
#' ilm_plot_horizon(ar = ar, iid = iid)
#' @export
ilm_plot_horizon <- function(..., rule = NULL, level = 0.95) {
  tabs <- list(...)
  if (!length(tabs)) stop("give at least one table of scores.", call. = FALSE)
  labels <- vapply(as.list(substitute(list(...)))[-1L], function(e)
    paste(deparse(e), collapse = ""), "")
  if (!is.null(names(tabs))) labels[nzchar(names(tabs))] <- names(tabs)[nzchar(names(tabs))]
  rule <- lmn_default_rule(tabs[[1]], rule)
  rows <- list()
  for (k in seq_along(tabs)) {
    t <- tabs[[k]]
    if (!is.data.frame(t) || !"horizon" %in% names(t) || !rule %in% names(t))
      stop("'", labels[k], "' must be a table of scores with a 'horizon' column ",
           "and a '", rule, "' column.", call. = FALSE)
    for (h in sort(unique(t$horizon))) {
      i <- which(t$horizon == h & is.finite(t[[rule]]))
      v <- t[[rule]][i]
      r <- if ("origin" %in% names(t) && length(unique(t$origin[i])) >= 2L)
        lmn_se_clustered(v, t$origin[i]) else lmn_se_rows(v)
      tq <- stats::qt((1 + level) / 2, r$df)
      rows[[length(rows) + 1L]] <- data.frame(model = labels[k], horizon = h,
        mean = r$estimate, se = r$se, lower = r$estimate - tq * r$se,
        upper = r$estimate + tq * r$se, n = length(i),
        origins = if ("origin" %in% names(t)) length(unique(t$origin[i])) else NA_integer_)
    }
  }
  d <- do.call(rbind, rows)
  m <- length(tabs)
  offset <- if (m > 1L) (match(d$model, labels) - (m + 1) / 2) * 0.25 / m else 0
  ## tinyplot titles the legend with the `by` expression itself, so it is a
  ## plain variable named for what it is
  model <- if (m > 1L) factor(d$model, levels = labels) else NULL
  tinyplot::tinyplot(x = d$horizon + offset, y = d$mean, ymin = d$lower, ymax = d$upper,
                     by = model, type = "pointrange",
                     xaxb = sort(unique(d$horizon)), xlab = "horizon",
                     ylab = paste("mean", rule, "score (lower is better)"))
  tinyplot::tinyplot_add(type = "l")
  invisible(d)
}

#' Plot model comparisons with their intervals
#'
#' Draws each comparison from [ilm_compare()] as its mean score difference with a
#' confidence interval, against a line at zero. An interval that excludes zero
#' is a difference the data can tell apart. It is built from the same standard
#' error and degrees of freedom as the comparison's own test: clustered,
#' row-wise or Diebold-Mariano.
#'
#' @section Reading the plot:
#' Each row is one model compared with the reference, labelled with the score,
#' how the standard error was computed and, for forecasts, the horizon.
#' * A point left of the dashed line at zero means the model scored better than
#'   the reference (scores are losses); right of it, worse.
#' * If the row's interval lies wholly on one side of the line, the data tell
#'   the two models apart at the plot's `level`. If it crosses the line, they
#'   do not: the comparison cannot say which model predicts better.
#' * A wide interval means the comparison has little to go on: with a clustered
#'   or Diebold-Mariano standard error, that is usually too few groups, folds or
#'   origins.
#'
#' @param x A comparison from [ilm_compare()].
#' @param level Confidence level of the intervals.
#' @return The data behind the plot, invisibly: a data frame with one row per
#'   comparison, giving its `label`, the mean score `difference`, the interval's
#'   `lower` and `upper` limits, and `by`, how its standard error was computed.
#' @seealso [ilm_compare()], whose result this plots; [ilm_plot_horizon()] for
#'   scores by horizon.
#' @examples
#' set.seed(1)
#' g <- rep(1:10, each = 5)
#' a <- data.frame(group = g, log = rnorm(50, 1.2) + rnorm(10)[g])
#' b <- data.frame(group = g, log = a$log - 0.1 + rnorm(50, 0, 0.2))
#' ilm_plot_compare(ilm_compare(reference = a, candidate = b, by = "group"))
#' @export
ilm_plot_compare <- function(x, level = 0.95) {
  if (!inherits(x, "ilm_compare"))
    stop("`x` must be a comparison from ilm_compare().", call. = FALSE)
  tq <- stats::qt((1 + level) / 2, x$df)
  ## the same pair compared with different standard errors gets its own row
  how <- ifelse(x$by == "row", "rows independent", paste("by", x$by))
  if (!is.null(x$horizon))
    how <- ifelse(is.na(x$horizon), how, paste0(how, ", horizon ", x$horizon))
  d <- data.frame(label = paste0(x$model, " vs ", x$reference, " (", x$rule, ", ", how, ")"),
                  difference = x$difference, lower = x$difference - tq * x$se,
                  upper = x$difference + tq * x$se, by = x$by)
  ## rows by position, first comparison on top: tinyplot sorts a factor's
  ## levels, which would put the comparisons in alphabetical order instead
  K <- nrow(d); pos <- rev(seq_len(K))
  op <- tinyplot::tpar(las = 1, mar = c(4.5, max(4, 0.55 * max(nchar(d$label))), 2, 1))
  on.exit(tinyplot::tpar(op))
  ## zero is always on the axis, so the reference line is always drawn
  tinyplot::tinyplot(x = pos, y = d$difference, ymin = d$lower, ymax = d$upper,
                     type = "pointrange", flip = TRUE, xlab = "",
                     xaxb = seq_len(K), xaxl = function(v) d$label[K + 1L - v],
                     ylim = range(c(0, d$lower, d$upper)),
                     ylab = "difference in mean score (< 0: model better)")
  tinyplot::tinyplot_add(x = pos, y = d$difference, type = tinyplot::type_vline(v = 0),
                         lty = 2, col = "grey40")
  invisible(d)
}

## The reliability diagram's "Wilson intervals" are Wilson's (1927) score
## intervals: awaiting the primary source (see R/ilm_calibration.R).

#' Plot the calibration of predictive distributions
#'
#' For a continuous or count outcome, plots the observed coverage of the central
#' intervals against their nominal levels. A well-calibrated model's points sit
#' on the diagonal; points below it are intervals too narrow. For a binary
#' outcome it draws the reliability diagram instead: the observed rate against
#' the mean predicted probability, bin by bin, with Wilson intervals. The
#' calibration intercept and slope are given under the title.
#'
#' @section Reading the plot:
#' * **Coverage plot** (a continuous or count outcome): each point is one
#'   interval level, with the share of outcomes inside it and a confidence
#'   interval for that share. Points on the dashed diagonal mean the intervals
#'   cover what they claim. Below it, they cover too little: the predictions
#'   are too certain. Above it, they cover too much: too uncertain.
#' * **Reliability diagram** (a binary outcome): each point is a bin of rows
#'   with similar predicted probabilities, placed at their mean prediction and
#'   their observed rate. On the diagonal, the probabilities are right. A
#'   pattern flatter than the diagonal means probabilities too extreme, and a
#'   steeper one probabilities too moderate.
#'
#' The verdict that goes with the plot, with what to do about it, is printed by
#' [ilm_calibration()][ilm_calibration.ilm_pred].
#'
#' @param x The calibration of a predictive distribution, class
#'   `"ilm_pred_calibration"`.
#' @param level Confidence level of the intervals.
#' @return The data behind the plot, invisibly: for coverage, the calibration's
#'   coverage table (each level's coverage, standard error and mean width) with
#'   the plotted interval limits; for a binary outcome, its reliability table
#'   (each bin's mean prediction, observed rate and Wilson interval).
#' @seealso [ilm_calibration()][ilm_calibration.ilm_pred], whose result it
#'   draws; [ilm_plot_pit()] for the PIT histogram of the same outcomes.
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rnorm(400))
#' d$y <- rnbinom(400, size = 1.5, mu = exp(1 + 0.5 * d$x))
#' train <- d[1:200, ]; test <- d[201:400, ]
#' fit <- illume::ilm_model(y ~ x, data = train, family = "poisson", verbose = FALSE)
#' cal <- illume::ilm_calibration(ilm_predict(fit, test, target = "response",
#'                                            nsim = 200), test$y)
#' ilm_plot_calibration(cal)
#' @export
ilm_plot_calibration <- function(x, level = 0.95) {
  if (!inherits(x, "ilm_pred_calibration"))
    stop("`x` must be the calibration of a predictive distribution.", call. = FALSE)
  z <- stats::qnorm((1 + level) / 2)
  if (is.null(x$reliability)) {
    d <- x$coverage
    d$lower <- pmax(0, d$coverage - z * d$se); d$upper <- pmin(1, d$coverage + z * d$se)
    rng <- range(c(d$level, d$lower, d$upper))
    tinyplot::tinyplot(x = d$level, y = d$coverage, ymin = d$lower, ymax = d$upper,
                       type = "pointrange", xlim = rng, ylim = rng,
                       xlab = "nominal coverage of the central interval",
                       ylab = "observed coverage")
  } else {
    d <- x$reliability
    b <- x$binary
    sub <- sprintf("calibration intercept %s, slope %s",
                   lmn_fmt_sig(b$estimate[b$measure == "intercept"]),
                   lmn_fmt_sig(b$estimate[b$measure == "slope"]))
    tinyplot::tinyplot(x = d$predicted, y = d$observed, ymin = d$lower, ymax = d$upper,
                       type = "pointrange", xlim = c(0, 1), ylim = c(0, 1),
                       xlab = "mean predicted probability", ylab = "observed rate",
                       sub = sub)
  }
  tinyplot::tinyplot_add(x = 0.5, y = 0.5, type = tinyplot::type_abline(a = 0, b = 1),
                         lty = 2, col = "grey40")
  invisible(d)
}

## The score column to use: the one named, or the first score in the table.
#' @keywords internal
#' @noRd
lmn_default_rule <- function(tab, rule) {
  if (!is.null(rule)) return(rule)
  rule <- attr(tab, "rules")[1]
  if (is.null(rule))
    rule <- grep("^(log|crps|brier|pinball_.+|interval_.+|squared|absolute|mase|rmsse)$",
                 names(tab), value = TRUE)[1]
  if (is.null(rule) || is.na(rule))
    stop("name the score with `rule`.", call. = FALSE)
  rule
}
