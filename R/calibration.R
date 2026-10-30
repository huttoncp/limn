## ---- calibration of a predictive distribution ------------------------------------
##
## Whether outcomes fall where their predictive distributions say they should:
## the PITs' uniformity, coverage and width by level, and for a binary outcome
## the calibration intercept and slope (Van Calster et al. 2019's "weak
## calibration"). The worker behind ilm_calibration.ilm_pred()
## (R/ilm_calibration.R), limn's method for illume's generic.
##
## Coverage is read from the PITs, not from whether the outcome lies inside a
## quantile interval: for a discrete outcome the quantile interval covers more
## than its nominal level, because of ties, while the randomised PIT's coverage
## is exactly nominal when the model is right.
##
## The rows are treated as independent. Rows that share a group or an origin
## are not, and then the tests here are too confident.

#' @keywords internal
#' @noRd
lmn_pred_calibration <- function(pred, y, bins = 10L,
                                 levels = c(0.5, 0.8, 0.9, 0.95), seed = 1L) {
  lmn_check_response_pred(pred)
  bins <- lmn_check_count(bins, "bins", 2L)
  u <- ilm_pit(pred, y, seed = seed)
  ok <- is.finite(u)
  n <- sum(ok)
  if (n < 20L)
    stop("calibration needs at least 20 outcomes; there are ", n, ".", call. = FALSE)
  uo <- u[ok]

  ## uniformity: at least 5 expected in each bin
  b <- min(bins, max(2L, floor(n / 5)))
  br <- seq(0, 1, length.out = b + 1L)
  cnt <- tabulate(findInterval(uo, br, rightmost.closed = TRUE, all.inside = TRUE), b)
  chi <- sum((cnt - n / b)^2 / (n / b))
  hist <- data.frame(lower = br[-(b + 1L)], upper = br[-1L], count = cnt,
                     expected = n / b)

  ## coverage and width, level by level
  sub <- lmn_subset_pred(pred, which(ok))
  cov <- do.call(rbind, lapply(levels, function(lv) {
    inside <- abs(uo - 0.5) <= lv / 2
    w <- lmn_pred_quantile(sub, (1 + lv) / 2) - lmn_pred_quantile(sub, (1 - lv) / 2)
    data.frame(level = lv, coverage = mean(inside),
               se = sqrt(lv * (1 - lv) / n), width = mean(w))
  }))

  ## bias is judged on the outcomes' own scale, y minus its predictive mean. The
  ## mean PIT is not a bias measure: a predictive distribution with the right
  ## mean and the wrong shape moves it too (Poisson on overdispersed counts
  ## puts too many outcomes below its centre).
  yo <- as.numeric(y)[ok]
  r <- yo - lmn_pred_mean(sub)
  res <- list(pit = u, n = n, histogram = hist,
              uniformity = c(statistic = chi, df = b - 1L,
                             p_value = stats::pchisq(chi, b - 1L, lower.tail = FALSE)),
              bias = c(mean_residual = mean(r), se = stats::sd(r) / sqrt(n)),
              coverage = cov, family = pred$dist$family)
  if (pred$dist$discrete && all(yo %in% c(0, 1)) &&
      all(abs(lmn_pred_cdf(sub, rep(1, n)) - 1) < 1e-10)) {
    p1 <- 1 - lmn_pred_cdf(sub, rep(0, n))
    res$binary <- lmn_binary_calibration(p1, yo)
    res$reliability <- lmn_reliability(p1, yo, bins)
  }
  res$checks <- lmn_calibration_checks(res)
  res$verdict <- lmn_verdict_lines(res$checks)
  structure(res, class = "ilm_pred_calibration")
}

## Calibration intercept (calibration-in-the-large: a logistic model with the
## predicted log-odds as an offset) and slope (the predicted log-odds as the
## only covariate), with Wald 95% intervals. 0 and 1 are ideal.
#' @keywords internal
#' @noRd
lmn_binary_calibration <- function(p, y) {
  lp <- stats::qlogis(pmin(pmax(p, 1e-12), 1 - 1e-12))
  f0 <- suppressWarnings(stats::glm(y ~ offset(lp), family = stats::binomial()))
  f1 <- suppressWarnings(stats::glm(y ~ lp, family = stats::binomial()))
  ci <- function(est, se) est + c(-1, 1) * stats::qnorm(0.975) * se
  s0 <- summary(f0)$coefficients; s1 <- summary(f1)$coefficients
  data.frame(measure = c("intercept", "slope"), ideal = c(0, 1),
             estimate = c(s0[1, 1], s1[2, 1]),
             lower = c(ci(s0[1, 1], s0[1, 2])[1], ci(s1[2, 1], s1[2, 2])[1]),
             upper = c(ci(s0[1, 1], s0[1, 2])[2], ci(s1[2, 1], s1[2, 2])[2]),
             observed_rate = mean(y), predicted_rate = mean(p))
}

## The reliability table of a binary outcome: equal-count bins of the predicted
## probability, with each bin's mean prediction, observed rate, and Wilson
## 95% interval. Bins hold at least 20 outcomes, so there are fewer than `bins`
## when the outcomes are few.
#' @keywords internal
#' @noRd
lmn_reliability <- function(p, y, bins = 10L) {
  b <- max(1L, min(bins, floor(length(p) / 20)))
  cuts <- stats::quantile(p, seq_len(b - 1L) / b, names = FALSE)
  bin <- findInterval(p, cuts, left.open = TRUE) + 1L
  do.call(rbind, lapply(sort(unique(bin)), function(k) {
    i <- bin == k; n <- sum(i); r <- mean(y[i]); z <- stats::qnorm(0.975)
    ctr <- (r + z^2 / (2 * n)) / (1 + z^2 / n)
    half <- z * sqrt(r * (1 - r) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
    data.frame(bin = k, n = n, predicted = mean(p[i]), observed = r,
               lower = ctr - half, upper = ctr + half)
  }))
}

## A plain verdict, and a remedy that exists: in illume, where most remedies
## for poor prediction are, or illume's ilm_remedies() as the general route.
## The words are the registry's (R/diagnostics.R).
#' @keywords internal
#' @noRd
lmn_calibration_checks <- function(x) {
  if (!is.null(x$binary)) return(lmn_binary_checks(x))
  out <- list()
  cv <- x$coverage
  z <- (cv$coverage - cv$level) / cv$se
  bz <- x$bias[["mean_residual"]] / x$bias[["se"]]
  fam <- if (is.null(x$family)) "" else x$family
  pct <- lmn_fmt_pct                                    # a share, to one decimal
  ## a biased mean comes first: it alone can make the intervals cover too
  ## little, so a shortfall beside it waits on the mean
  biased <- abs(bz) > 2.58
  if (biased)
    out <- c(out, list(lmn_diag("calib_bias", size = lmn_fmt_sig(abs(x$bias[["mean_residual"]])),
                                direction = if (bz > 0) "above" else "below",
                                se = lmn_fmt_sig(x$bias[["se"]]))))
  if (any(z < -2.58)) {
    k <- which.min(z)
    id <- if (biased) "calib_narrow_bias" else if (fam == "poisson") "calib_narrow_counts"
      else if (fam %in% c("nbinom", "gaussian", "beta")) "calib_narrow_dispersion" else "calib_narrow"
    out <- c(out, list(lmn_diag(id, coverage = pct(cv$coverage[k]), level = round(100 * cv$level[k]))))
  }
  if (any(z > 2.58)) {
    k <- which.max(z)
    out <- c(out, list(lmn_diag("calib_wide", coverage = pct(cv$coverage[k]),
                                level = round(100 * cv$level[k]))))
  }
  if (!length(out)) out <- list(lmn_diag("calib_ok", p = lmn_fmt_p(x$uniformity[["p_value"]])))
  do.call(rbind, out)
}

## A binary outcome's predictive distribution is one probability, so it is
## judged by its calibration intercept (are the probabilities right on
## average?) and slope (are they too extreme or too moderate?), which together
## assess what Van Calster et al. (2019, p. 4) call weak calibration. Interval
## coverage adds nothing for a 0/1 outcome.
#' @keywords internal
#' @noRd
lmn_binary_checks <- function(x) {
  out <- list()
  b <- x$binary
  it <- b[b$measure == "intercept", ]; s <- b[b$measure == "slope", ]
  ci <- function(r) list(est = lmn_fmt_sig(r$estimate), lo = lmn_fmt_sig(r$lower),
                         hi = lmn_fmt_sig(r$upper))
  if (it$lower > 0 || it$upper < 0)
    out <- c(out, list(do.call(lmn_diag, c(list("binary_intercept",
      observed = lmn_fmt_pct(it$observed_rate),
      predicted = lmn_fmt_pct(it$predicted_rate)), ci(it)))))
  ## a slope sure to be below 0 is reversed, not just too extreme
  if (s$upper < 0) out <- c(out, list(do.call(lmn_diag, c(list("binary_reversed"), ci(s)))))
  else if (s$upper < 1) out <- c(out, list(do.call(lmn_diag, c(list("binary_overconfident"), ci(s)))))
  if (s$lower > 1) out <- c(out, list(do.call(lmn_diag, c(list("binary_underconfident"), ci(s)))))
  if (!length(out)) out <- list(lmn_diag("binary_ok"))
  do.call(rbind, out)
}

#' @export
print.ilm_pred_calibration <- function(x, ...) {
  cat("<calibration> ", lmn_fmt_count(x$n), " outcomes", if (!is.null(x$family))
    paste0(", ", x$family, " predictive distributions"), "\n", sep = "")
  cat(sprintf("  PIT uniformity: chi-square %s on %d df, %s\n",
              lmn_fmt_sig(x$uniformity[["statistic"]]), as.integer(x$uniformity[["df"]]),
              lmn_fmt_p(x$uniformity[["p_value"]])))
  if (is.null(x$binary)) {
    cat("  coverage by level (from the PITs):\n")
    for (k in seq_len(nrow(x$coverage)))
      ## the coverage's standard error, an estimate, in percentage points
      cat(sprintf("    %3.0f%%: %6s (SE %s percentage points), mean width %s\n",
                  100 * x$coverage$level[k], lmn_fmt_pct(x$coverage$coverage[k]),
                  lmn_fmt_sig(100 * x$coverage$se[k]), lmn_fmt_sig(x$coverage$width[k])))
  } else {
    for (k in seq_len(nrow(x$binary)))
      cat(sprintf("  calibration %s: %s (95%% interval %s to %s; ideal %g)\n",
                  x$binary$measure[k], lmn_fmt_sig(x$binary$estimate[k]),
                  lmn_fmt_sig(x$binary$lower[k]), lmn_fmt_sig(x$binary$upper[k]),
                  x$binary$ideal[k]))
  }
  cat(paste0("  ", x$verdict), sep = "\n")
  invisible(x)
}
