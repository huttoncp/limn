## ---- totals over rows, from joint draws ---------------------------------------------
##
## A total over rows is summed draw by draw. So whatever the rows share within a
## draw (the parameters, a unit's effects, a series' path) carries into the
## total, and totals at different levels add up. A total of new observations
## stays an exact mixture where a sum stays in its family (Poisson, gaussian);
## otherwise it is summed from outcomes simulated from each draw.

## The family whose sums stay in the family, or NULL.
#' @keywords internal
#' @noRd
lmn_sum_closed <- function(dist) {
  if (!is.null(dist$zero) || !is.null(dist$categories)) return(NULL)
  if (identical(dist$family, "poisson")) return("poisson")
  if (identical(dist$family, "gaussian")) return("gaussian")
  NULL
}

## Each row's family mean, draw by draw: rows x draws.
#' @keywords internal
#' @noRd
lmn_draw_means <- function(pred) {
  n <- lmn_pred_n(pred); S <- length(pred$weights); i <- seq_len(n)
  out <- matrix(0, n, S)
  for (s in seq_len(S)) {
    a <- lmn_draw_args(pred, i, s)
    out[, s] <- pred$dist$mean(a$lp, a$par)
  }
  out
}

## Each row's log dispersion, draw by draw: the dispersion model's linear
## predictor where there is one, the family's own parameter otherwise.
#' @keywords internal
#' @noRd
lmn_draw_logdisp <- function(pred) {
  if (!is.null(pred$lp$disp)) return(pred$lp$disp)
  v <- pred$par$logdisp
  if (is.matrix(v)) v else matrix(v, lmn_pred_n(pred), length(v), byrow = TRUE)
}

## The total, or with `fun = "mean"` the average, of `pred`'s rows within each
## group. `groups` gives each row's group, 1 to K, and `rows` describes the K
## results in that order. An average is the total over the number of rows: done
## directly where the family allows (draws of means, gaussian), and otherwise
## kept as the total's distribution with a `scale` of 1 / rows, so that an
## average of counts is scored through its whole-number total.
#' @keywords internal
#' @noRd
lmn_aggregate <- function(pred, groups, rows, fun = "sum", outcomes = 4L, seed = 1L) {
  meta <- pred$meta
  add <- function(M) rowsum(M, groups, reorder = TRUE)
  m <- tabulate(groups, nrow(rows))
  avg <- identical(fun, "mean")
  how <- function(txt) paste0(txt, if (avg) "; the average is that total over the number of rows")
  meta$aggregate_fun <- fun
  scaled <- function(p) { if (avg) { p$scale <- 1 / m; p$estimate <- p$estimate / m }; p }
  if (lmn_is_draws(pred)) {
    if (pred$kind == "link")
      stop("a total of linear predictors has no meaning on the outcome's scale. ",
           "Aggregate a prediction of the mean (target = \"mean\") or of new ",
           "observations (target = \"response\").", call. = FALSE)
    d <- if (avg) 1 / m else rep(1, length(m))
    meta$aggregate <- if (avg) "each draw's means, averaged within the draw"
      else "each draw's means, summed within the draw"
    meta$aggregate_method <- "draw_means"
    out <- lmn_new_draw_pred(add(pred$draws) * d, as.numeric(add(pred$estimate)) * d,
                             kind = "mean", rows = rows, meta = meta)
    out$what <- if (avg) "an average of means" else "a total of means"
    return(out)
  }
  what <- if (avg) "an average of new observations" else "a total of new observations"
  if (identical(pred$kind, "sample")) {
    meta$aggregate <- how("each sample's values, summed within the sample")
    meta$aggregate_method <- "samples"
    return(scaled(lmn_new_sample_pred(add(pred$draws), pred$dist$discrete, rows, meta,
                                      pred$dist$family, what = what)))
  }
  dist <- pred$dist
  closed <- lmn_sum_closed(dist)
  if (!is.null(closed)) {
    tot <- add(lmn_draw_means(pred))
    if (closed == "poisson") {
      meta$aggregate <- how("exact: within a draw, a sum of Poisson outcomes is Poisson")
      meta$aggregate_method <- "exact_poisson"
      out <- lmn_new_pred(dist, list(mu = log(tot)), weights = pred$weights, rows = rows,
                          meta = meta)
      out$what <- what
      if (avg) out$scale <- 1 / m
      return(out)
    }
    ## gaussian: the means and the variances add; an average divides both scales
    ld <- 0.5 * log(add(exp(2 * lmn_draw_logdisp(pred))))
    if (avg) { tot <- tot / m; ld <- ld - log(m) }
    meta$aggregate <- if (avg) "exact: within a draw, an average of gaussian outcomes is gaussian"
      else "exact: within a draw, a sum of gaussian outcomes is gaussian"
    meta$aggregate_method <- "exact_gaussian"
    lp <- list(mu = tot); par <- list()
    if ("disp" %in% dist$lp_names) lp$disp <- ld else par$logdisp <- ld
    out <- lmn_new_pred(dist, lp, par = par, weights = pred$weights, rows = rows, meta = meta)
    out$what <- what
    return(out)
  }
  ## otherwise: outcomes simulated from each draw, summed within it; unequal
  ## draw weights are honoured by resampling the draws in proportion to them
  S <- length(pred$weights); i <- seq_len(lmn_pred_n(pred))
  sims <- lmn_with_seed(seed, {
    idx <- if (diff(range(pred$weights)) > 1e-12)
      sample.int(S, S, replace = TRUE, prob = pred$weights) else seq_len(S)
    out <- matrix(0, nrow(rows), S * outcomes)
    for (k in seq_len(S)) {
      a <- lmn_draw_args(pred, i, idx[k])
      for (j in seq_len(outcomes))
        out[, (k - 1L) * outcomes + j] <- add(dist$r(a$lp, a$par))
    }
    out
  })
  ## a family that cannot draw from some parameter draws: typically a dispersion
  ## at its limit that the fit did not hold
  if (anyNA(sims))
    lmn_signal("missing_totals", bad = sum(colSums(is.na(sims)) > 0), S = ncol(sims),
               type = "error")
  meta$aggregate <- how(sprintf("simulated: %d outcomes from each draw, summed within it",
                                outcomes))
  meta$aggregate_method <- "simulated"; meta$aggregate_outcomes <- as.integer(outcomes)
  scaled(lmn_new_sample_pred(sims, dist$discrete, rows, meta, dist$family, what = what))
}
