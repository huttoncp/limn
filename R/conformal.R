## ---- conformal recalibration ----------------------------------------------------
##
## Recalibration keeps the model's own predictive distribution and changes only
## which of its central intervals is reported. From out-of-sample PITs u, the
## score s = |u - 0.5| says how far out an outcome fell: the central L interval
## covers exactly when s <= L / 2. So the (1 - alpha) conformal quantile q of
## the calibration scores gives a recalibrated level L' = 2q, and the model's
## central L' interval covers 1 - alpha. The PITs are randomised for discrete
## outcomes, so the scores are continuous, as the theory below needs.
##
## Which quantile depends on what the prediction is for (prototype P11 in
## studies/prototypes/, and Dunn, Wasserman & Ramdas 2023, JASA 118, 2491-2502):
##
##   seen units   the split-conformal quantile over rows. Rows from different
##                groups are not exchangeable, so there is no finite-sample
##                guarantee here, and no method with one borrows strength
##                across groups (Dunn et al. call it an open problem). It is
##                calibrated within bins of the predicted mean ("Mondrian"),
##                because one level for all rows over-covers some and
##                under-covers others (P11a: 0.956 against 0.825).
##   new units    Dunn et al.'s CDF pooling: each group's empirical CDF of the
##                scores averaged with equal weight, valid as the number of
##                groups grows. Or repeated subsampling, one row per group B
##                times with the p-values averaged, which guarantees 1 - 2 alpha
##                and is close to 1 - alpha in practice.
##   forecasts    a level per horizon, from rolling-origin scores.
##
## Internal: ilm_cv() supplies the out-of-sample PITs.

## The split-conformal quantile: the ceiling((n + 1)(1 - alpha))-th smallest
## score, infinite when that is more than n (too few scores for the level).
#' @keywords internal
#' @noRd
lmn_conformal_quantile <- function(s, alpha) {
  s <- sort(s[is.finite(s)])
  n <- length(s); k <- ceiling((n + 1) * (1 - alpha))
  if (k > n) Inf else s[k]
}

## Dunn et al.'s pooled CDF: the average over groups of each group's empirical
## CDF of the scores; the level is its (1 - alpha) quantile.
#' @keywords internal
#' @noRd
lmn_pooled_quantile <- function(s, group, alpha) {
  ok <- is.finite(s)
  s <- s[ok]; group <- as.character(group[ok])
  w <- 1 / as.numeric(table(group)[group])            # each group weighs the same
  w <- w / sum(w)
  o <- order(s)
  cw <- cumsum(w[o])
  s[o][which(cw >= 1 - alpha - 1e-12)[1]]
}

## Repeated subsampling (Dunn et al., Method 3): B subsamples of one score per
## group. For a threshold t, subsample b's conformal p-value is
## (1 + #{scores >= t}) / (k + 1); the set is every t whose average p-value
## exceeds alpha, and the level is its largest t.
#' @keywords internal
#' @noRd
lmn_subsample_quantile <- function(s, group, alpha, B = 100L, seed = 1L) {
  ok <- is.finite(s)
  s <- s[ok]; group <- as.character(group[ok])
  idx <- split(seq_along(s), group)
  k <- length(idx)
  if (k <= 1 / alpha - 1) return(Inf)                  # too few groups: unbounded
  sub <- lmn_with_seed(seed, vapply(seq_len(B), function(b)
    vapply(idx, function(ix) s[ix[sample.int(length(ix), 1L)]], 0), numeric(k)))
  sub <- matrix(sub, nrow = k)
  ## the average p-value is a step function of t, falling at the scores; the
  ## largest t with average p-value > alpha is one of them
  cand <- sort(unique(as.vector(sub)))
  pbar <- vapply(cand, function(t) mean((1 + colSums(sub >= t)) / (k + 1)), 0)
  keep <- which(pbar > alpha)
  if (!length(keep)) return(0)
  cand[max(keep)]
}

## Recalibrated levels, by bin: the level at which the model's own central
## intervals cover `level` out of sample.
##   method "rows"      split conformal over rows (seen units)
##          "pooled"    Dunn et al.'s CDF pooling over `group` (new units)
##          "subsample" repeated subsampling over `group` (new units, guaranteed)
##   by     NULL for one level; a factor of bins (Mondrian), such as predicted-
##          mean thirds or the forecast horizon
#' @keywords internal
#' @noRd
lmn_conformal_levels <- function(pit, level = 0.9, method = c("rows", "pooled", "subsample"),
                                 group = NULL, by = NULL, B = 100L, seed = 1L) {
  method <- match.arg(method)
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1)
    stop("`level` must be one number between 0 and 1.", call. = FALSE)
  if (method != "rows" && is.null(group))
    stop("method \"", method, "\" calibrates over groups, so it needs `group`.",
         call. = FALSE)
  alpha <- 1 - level
  s <- abs(pit - 0.5)
  bins <- if (is.null(by)) factor(rep("all", length(s))) else as.factor(by)
  res <- lapply(levels(bins), function(b) {
    i <- which(bins == b)
    q <- switch(method,
      rows = lmn_conformal_quantile(s[i], alpha),
      pooled = lmn_pooled_quantile(s[i], group[i], alpha),
      subsample = lmn_subsample_quantile(s[i], group[i], alpha, B, seed))
    data.frame(bin = b, n = length(i),
               groups = if (is.null(group)) NA_integer_ else length(unique(group[i])),
               quantile = q, level = min(2 * q, 1), unbounded = is.infinite(q))
  })
  structure(do.call(rbind, res), class = c("lmn_conformal", "data.frame"),
            target = level, method = method)
}

## Mondrian bins for the predicted means: cut points at the calibration rows'
## quantiles, so new rows fall into the same bins. Thirds by default, fewer
## when a bin would hold too few scores for the level.
#' @keywords internal
#' @noRd
lmn_mean_bins <- function(m_cal, level = 0.9, nbins = 3L) {
  need <- ceiling(1 / (1 - level)) - 1L                # fewer scores: unbounded
  nbins <- max(1L, min(nbins, floor(length(m_cal) / (2 * need))))
  cuts <- if (nbins > 1L)
    stats::quantile(m_cal, seq_len(nbins - 1L) / nbins, names = FALSE) else numeric(0)
  list(cuts = cuts, bin = function(m) factor(findInterval(m, cuts) + 1L,
                                             levels = seq_len(nbins)))
}
