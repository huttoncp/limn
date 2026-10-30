#' Weighted averages of predicted means: standardisation and post-stratification
#'
#' Averages a prediction of the mean over groups of its rows, with weights,
#' draw by draw. With the rows a reference population and equal weights, it is
#' standardisation: the population's average mean under the model, for
#' example with a covariate set to each of several values. With the rows
#' population cells and their sizes as weights, it is post-stratification (as
#' in multilevel regression and post-stratification, MRP).
#'
#' The average is taken within each draw, so its uncertainty is the joint
#' uncertainty of the rows' means: whatever they share within a draw (the
#' parameters, a unit's effect) carries into the average.
#'
#' Post-stratification weights each cell's predicted mean by the cell's
#' population count, as Gelman and Hill (2007) do for state opinion from
#' national polls (section 14.1, equation (14.1)); their code in the same
#' section computes the weighted average in each simulation draw, as here.
#'
#' @section Means, not new observations:
#' Standardisation and post-stratification average what the model expects,
#' so `pred` must be a prediction of the mean (`target = "mean"`). A weighted
#' average of new observations is a different question, whose answer depends
#' on the sizes of the groups being averaged. For an unweighted average of new
#' observations, use [ilm_aggregate()] with `fun = "mean"`.
#'
#' @param pred A prediction of the mean, from [ilm_predict()] or
#'   [ilm_forecast()] with `target = "mean"`.
#' @param by Names of columns in the prediction's rows, as `as.data.frame()`
#'   shows them. There is one average for each combination of their values.
#'   `NULL`, the default, averages every row. Rows of different scenarios are
#'   never averaged together.
#' @param weights The rows' weights: a non-negative number per row, or the
#'   name of a column in the prediction's rows. They are normalised within
#'   each group, so only their ratios matter. `NULL`, the default, weights the
#'   rows equally.
#' @return A prediction of the mean (see *The prediction object* in
#'   [ilm_predict()]), with one row per combination of `by` (and of scenario),
#'   `n_rows`, the number of rows each average combines, and `weight`, their
#'   total weight.
#' @references Gelman, A. and Hill, J. (2007). *Data Analysis Using Regression
#'   and Multilevel/Hierarchical Models*. Cambridge University Press. Section
#'   14.1, "State-level opinions from national polls": equation (14.1), and
#'   the weighted average computed in each simulation draw.
#' @seealso [ilm_predict()] with `target = "mean"` for the rows' means;
#'   [ilm_aggregate()] for totals, and averages of new observations.
#' @examples
#' set.seed(1)
#' d <- data.frame(clinic = factor(rep(sprintf("c%02d", 1:20), each = 10)),
#'                 age = rep(c("young", "old"), 100), x = rnorm(200))
#' d$visits <- rpois(200, exp(0.5 + 0.4 * d$x + 0.3 * (d$age == "old") +
#'                            rnorm(20, 0, 0.4)[d$clinic]))
#' fit <- illume::ilm_model(visits ~ x + age + (1 | clinic), data = d,
#'                          family = "poisson", verbose = FALSE)
#' ## post-stratification: the mean for a population that is 70% old (the
#' ## population average does not use the clinic, but the rows need its column)
#' cells <- data.frame(age = c("young", "old"), x = 0, size = c(300, 700), clinic = NA)
#' p <- ilm_predict(fit, cells, target = "mean", groups = "population", nsim = 500)
#' as.data.frame(ilm_average(p, weights = "size"))
#' @export
ilm_average <- function(pred, by = NULL, weights = NULL) {
  if (!inherits(pred, "ilm_pred"))
    stop("`pred` must be a prediction from ilm_predict() or ilm_forecast(); it is ",
         class(pred)[1], ".", call. = FALSE)
  if (!identical(pred$kind, "mean"))
    stop("`pred` must be a prediction of the mean (target = \"mean\"); this one is of ",
         if (identical(pred$kind, "link")) "the linear predictor" else "new observations",
         ". Standardisation and post-stratification average what the model expects. ",
         "For an unweighted average of new observations, use ilm_aggregate(fun = \"mean\").",
         call. = FALSE)
  rows <- pred$rows; n <- nrow(rows)
  w <- lmn_average_weights(weights, rows, n)
  g <- lmn_agg_groups(rows, by)
  tot_w <- as.numeric(rowsum(w, g$groups, reorder = TRUE))
  if (any(tot_w <= 0))
    stop("the weights of ", sum(tot_w <= 0), " group(s) sum to zero, so their average ",
         "is undefined; give each group a positive weight.", call. = FALSE)
  d <- 1 / tot_w
  out <- g$out; out$weight <- tot_w
  meta <- pred$meta
  meta$aggregate <- if (is.null(weights)) "each draw's means, averaged within the draw"
    else "each draw's means, a weighted average within the draw"
  meta$aggregate_method <- "draw_means"
  meta$aggregate_fun <- if (is.null(weights)) "mean" else "weighted_mean"
  res <- lmn_new_draw_pred(rowsum(pred$draws * w, g$groups, reorder = TRUE) * d,
                           as.numeric(rowsum(pred$estimate * w, g$groups, reorder = TRUE)) * d,
                           kind = "mean", rows = out, meta = meta)
  res$what <- if (is.null(weights)) "an average of means" else "a weighted average of means"
  res$checks <- pred$checks                      # an average keeps the prediction's checks
  res$origin <- lmn_agg_origin(pred$origin, g$groups, nrow(out))
  res
}

## The rows' weights: equal, a numeric vector, or a column of the rows.
#' @keywords internal
#' @noRd
lmn_average_weights <- function(weights, rows, n) {
  if (is.null(weights)) return(rep(1, n))
  if (is.character(weights) && length(weights) == 1L) {
    if (!weights %in% names(rows))
      stop("`weights` names a column the prediction's rows do not have: ", weights,
           ". They have: ", paste(names(rows), collapse = ", "), ".", call. = FALSE)
    weights <- rows[[weights]]
  }
  if (!is.numeric(weights) || length(weights) != n)
    stop("`weights` must be one number per row of the prediction (", n,
         "), or the name of a column of its rows.", call. = FALSE)
  if (anyNA(weights) || any(weights < 0))
    stop("`weights` must be non-negative and not missing.", call. = FALSE)
  as.numeric(weights)
}
