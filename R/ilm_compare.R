#' Compare models' predictive scores, with standard errors that respect the data
#'
#' Compares the scores of two or more models on the same rows: the mean
#' difference in a score, its standard error, and a test of whether it is zero.
#' The first model is the reference. Every other model is compared with it, as
#' that model's mean score minus the reference's. Scores are losses, so a
#' negative difference means the model predicts better than the reference.
#'
#' @section Why the standard error depends on `by`:
#' Scores on rows from the same group, fold or forecast origin are correlated,
#' because those rows share whatever the model got right or wrong about them.
#' Treating them as independent understates the standard error, and makes a
#' difference look more certain than it is. So the standard error is computed
#' at the level the data are dependent at:
#'
#' * `by = "row"`: rows are independent. Right for rows with no shared
#'   structure.
#' * `by` = a column such as `"group"` or `"fold"`: a cluster-robust standard
#'   error, with the rows of a cluster summed first. This is the same as
#'   `sandwich::vcovCL()` with its default HC1 adjustment, applied to the
#'   differences. The test uses a t distribution with one fewer degrees of
#'   freedom than there are clusters.
#' * `by = "origin"`: rolling-origin forecasts. The differences are averaged
#'   within each origin, and those origin means are tested with a
#'   Diebold-Mariano test: pooled over the horizons, and at each horizon on its
#'   own when there are several. Forecasts from neighbouring origins overlap
#'   when the horizon is longer than one, so the origin means are
#'   autocorrelated, and the test's variance allows for it.
#'   - `method = "wpe"` (the default): the fixed-smoothing test of Coroneo and
#'     Iacone (2020). The variance is the periodogram averaged over the first
#'     m Fourier frequencies of the T origin means (their weighted periodogram
#'     estimate with the Daniell kernel, section 4.2, equation (8)), and the
#'     statistic is referred to t on 2m degrees of freedom (equation (9)).
#'     m = floor(T^(1/3)) is the bandwidth they use in their simulations
#'     (section 5). In limn's
#'     simulation (`studies/findings/compare_dm_fixed.md`) it held its size
#'     from 10 origins up, where the usual test rejected up to 12% of the time
#'     at 5%.
#'   - `method = "dm"`: the usual test, a long-run variance truncated at `lag`
#'     with equal weights, the uniform lag window Diebold and Mariano (1995,
#'     section 1.1) recommend, with Harvey,
#'     Leybourne and Newbold's (1997) small-sample correction (equation (9),
#'     p. 283), referred to t on T - 1 degrees of freedom (pp. 283-284). When
#'     no fallback is needed it agrees with `forecast::dm.test()`, as limn's
#'     study `studies/scripts/compare_dm.R` checks; that package is an outside
#'     comparator, not a source. By default `lag` is the longest horizon minus
#'     one, the overlap when the origin moves one time point at a time.
#'     Diebold and Mariano take k - 1 as a benchmark for k-step-ahead
#'     forecasts, because optimal forecasts' errors are at most
#'     (k - 1)-dependent (section 1.1); they note that this can fail in
#'     practice. With
#'     fewer than about 100 origins over several horizons it rejects too often,
#'     and warns.
#'
#' @eval lmn_diag_rd(c("dm_few_origins", "dm_negative_variance"))
#' @section Catching the checks:
#' Each check above is raised as a warning of class `lmn_<id>`, beside
#' `lmn_diagnostic`, and is kept on the result as a row of `attr(x, "checks")`,
#' a data frame whose `id` column names it:
#' * `dm_few_origins` (warning class `lmn_dm_few_origins`): few origins over
#'   several horizons for the usual test;
#' * `dm_negative_variance` (warning class `lmn_dm_negative_variance`): the
#'   truncated long-run variance was negative, so it was taken as 0 and the
#'   null hypothesis rejected.
#'
#' So `withCallingHandlers(ilm_compare(...), lmn_dm_negative_variance =
#' function(w) ...)` handles that warning alone, and `attr(cmp, "checks")$id` says which checks
#' a comparison raised.
#' @param ... Two or more tables of scores for the same rows, in the same
#'   order, such as the output of the scoring of an `"ilm_pred"`. Name them
#'   (`ilm_compare(ar = s1, iid = s2)`) to label the comparison.
#' @param rule The score to compare: the name of a score column. By default
#'   the first score in the table, which for a plain data frame is its first
#'   column named after a rule (`log`, `crps`, `brier`, `pinball_*`,
#'   `interval_*`, `squared`, `absolute`, `mase` or `rmsse`).
#' @param by `"row"`, `"origin"`, or the name of a column the rows are
#'   clustered by, such as `"group"` or `"fold"`. See *Why the standard error
#'   depends on `by`*.
#' @param lag Only for `by = "origin"` with `method = "dm"`: the truncation
#'   lag of the long-run variance, that is, how many autocovariances of the
#'   per-origin mean differences it keeps.
#'   - **The default** suits origins one time step apart. Forecasts h steps
#'     ahead from neighbouring origins then overlap in h - 1 steps, so the
#'     pooled row uses the longest horizon minus one, and horizon h's own row
#'     uses h - 1.
#'   - **Set it yourself when the origins are further apart.** With origins s
#'     time steps apart, forecasts h ahead overlap in ceiling(h / s) - 1 of
#'     them, so use `lag = ceiling(H / s) - 1` for a longest horizon H. It is
#'     0 when the origins are at least H apart, and the forecasts do not
#'     overlap.
#'   - A value you give applies to the pooled row and to every horizon's row.
#'     It must be between 0 and the number of origins minus 1.
#'   - If the variance it gives is negative, it is taken as 0 and the null
#'     hypothesis rejected (statistic infinite, p-value 0), as Diebold and
#'     Mariano (1995, section 1.1) prescribe, and the test warns.
#'
#'   With `method = "wpe"`, or any `by` other than `"origin"`, giving `lag` is
#'   an error. The fixed-smoothing test chooses its own bandwidth from the
#'   number of origins, whatever their spacing.
#' @param method Only for `by = "origin"`: which Diebold-Mariano test.
#'   - `"wpe"`, the default: the fixed-smoothing test of Coroneo and Iacone
#'     (2020, section 4.2). The variance of the T per-origin mean differences is
#'     their periodogram averaged over the first m = floor(T^(1/3)) Fourier
#'     frequencies (equation (8); the bandwidth of their simulations, section
#'     5). The statistic is referred to t on 2m degrees of freedom (equation
#'     (9)).
#'     It needs at least 3 origins and takes no `lag`. In limn's simulation it
#'     held its size from 10 origins up.
#'   - `"dm"`: the usual test. It uses a long-run variance truncated at `lag`
#'     with equal weights (Diebold and Mariano 1995, section 1.1), Harvey,
#'     Leybourne and Newbold's (1997) correction (equation (9), p. 283), and t
#'     on T - 1 degrees of freedom (pp. 283-284). Use it to compare with
#'     published results. With fewer than about 100 origins over several
#'     horizons it rejects too often, and warns.
#'
#'   Giving `method` with any `by` other than `"origin"` is an error. See *Why
#'   the standard error depends on `by`*.
#' @return A data frame of class `"ilm_compare"`, one row per comparison, with
#'   the two mean scores, their difference, its standard error, the test
#'   statistic, degrees of freedom and p-value, the number of rows and of
#'   clusters, and the method. For `by = "origin"` with several horizons, each
#'   comparison has a row pooled over the horizons (`horizon` is `NA`) and a row
#'   for each horizon.
#'
#'   The test is also given as data, so nothing needs reading from `method`,
#'   which describes it in words:
#'   - `test`: `"independent"` (`by = "row"`), `"clustered"`, `"wpe"` or `"dm"`
#'     (`by = "origin"`).
#'   - `bandwidth`: for `"wpe"`, the number m of Fourier frequencies averaged;
#'     otherwise `NA`.
#'   - `lag`: for `"dm"`, the truncation lag used; otherwise `NA`.
#' @references Diebold, F. X. and Mariano, R. S. (1995). Comparing predictive
#'   accuracy. *Journal of Business and Economic Statistics* 13, 253-263.
#'   Section 1.1: the statistic, the truncated long-run variance with the
#'   uniform (rectangular) lag window they recommend, the lag k - 1 as a
#'   benchmark for k-step-ahead forecasts, and the rule
#'   for a negative estimate.
#'
#'   Coroneo, L. and Iacone, F. (2020). Comparing predictive accuracy in small
#'   samples using fixed-smoothing asymptotics. *Journal of Applied
#'   Econometrics* 35(4), 391-409. Section 4.2, equations (8) and (9): the
#'   weighted periodogram estimate and its t limit on 2m degrees of freedom;
#'   section 5: the bandwidth m = floor(T^(1/3)).
#'
#'   Harvey, D., Leybourne, S. and Newbold, P. (1997). Testing the equality of
#'   prediction mean squared errors. *International Journal of Forecasting* 13,
#'   281-291. Equation (9), p. 283: the corrected statistic; pp. 283-284:
#'   Student's t on n - 1 degrees of freedom.
#' @seealso [ilm_folds()] for folds that make scores comparable.
#' @examples
#' set.seed(1)
#' g <- rep(1:10, each = 5)
#' a <- data.frame(group = g, log = rnorm(50, 1.2) + rnorm(10)[g])
#' b <- data.frame(group = g, log = a$log - 0.1 + rnorm(50, 0, 0.2))
#' ilm_compare(reference = a, candidate = b, by = "group")
#' @export
ilm_compare <- function(..., rule = NULL, by = "row", lag = NULL, method = NULL) {
  tabs <- list(...)
  if (length(tabs) < 2L)
    stop("give at least two tables of scores to compare.", call. = FALSE)
  nm <- names(tabs)
  labels <- vapply(as.list(substitute(list(...)))[-1L], function(e)
    paste(deparse(e), collapse = ""), "")
  if (!is.null(nm)) labels[nzchar(nm)] <- nm[nzchar(nm)]
  if (!all(vapply(tabs, is.data.frame, TRUE)))
    stop("each argument must be a table of scores (a data frame).", call. = FALSE)
  n <- nrow(tabs[[1]])
  if (!all(vapply(tabs, nrow, 0L) == n))
    stop("the tables have different numbers of rows; compare scores for the ",
         "same rows.", call. = FALSE)
  ## scores from limn say which columns are scores; a plain table is read by
  ## the rules' own column names
  rule <- lmn_default_rule(tabs[[1]], rule)
  for (k in seq_along(tabs))
    if (!rule %in% names(tabs[[k]]))
      stop("'", labels[k], "' has no '", rule, "' column.", call. = FALSE)
  if (!identical(by, "row")) {
    for (k in seq_along(tabs))
      if (!by %in% names(tabs[[k]]))
        stop("'", labels[k], "' has no '", by, "' column to group the rows by.",
             call. = FALSE)
    for (k in seq_along(tabs)[-1L])
      if (!identical(as.vector(tabs[[k]][[by]]), as.vector(tabs[[1]][[by]])))
        stop("'", labels[k], "' and '", labels[1], "' disagree on '", by,
             "'; the rows must be the same rows, in the same order.", call. = FALSE)
  }

  if (!identical(by, "origin") && (!is.null(method) || !is.null(lag)))
    stop("`method` and `lag` apply to by = \"origin\" only; with by = \"", by,
         "\" the standard error is ", if (identical(by, "row")) "that of independent rows"
         else "clustered", ".", call. = FALSE)
  method <- if (is.null(method)) lmn_dm_default else match.arg(method, c("wpe", "dm"))
  if (!is.null(lag) && method != "dm")
    stop("`lag` applies to method = \"dm\" only; the fixed-smoothing test sets its own ",
         "bandwidth from the number of origins.", call. = FALSE)
  ref <- tabs[[1]][[rule]]
  ## the checks the comparison raises are warned as before, and kept on it
  col <- lmn_collect(lapply(seq_along(tabs)[-1L], function(k) {
    x <- tabs[[k]][[rule]]
    ok <- is.finite(x) & is.finite(ref)
    d <- x[ok] - ref[ok]
    row <- function(res, keep, horizon = NA_integer_)
      data.frame(model = labels[k], reference = labels[1], rule = rule, horizon = horizon,
                 model_mean = mean(x[ok][keep]), reference_mean = mean(ref[ok][keep]),
                 difference = res$estimate, se = res$se, statistic = res$statistic,
                 df = res$df, p_value = res$p, n = sum(keep), clusters = res$clusters,
                 by = by, method = res$method, test = res$test, bandwidth = res$bandwidth,
                 lag = res$lag, dropped = sum(!ok))
    all_rows <- rep(TRUE, length(d))
    if (identical(by, "row")) return(row(lmn_se_rows(d), all_rows))
    if (!identical(by, "origin")) return(row(lmn_se_clustered(d, tabs[[1]][[by]][ok]), all_rows))
    origin <- tabs[[1]]$origin[ok]; hz <- tabs[[1]]$horizon[ok]
    test <- function(keep)
      lmn_se_origin(d[keep], origin[keep], method,
                    if (is.null(lag)) lmn_default_lag(hz[keep]) else lag)
    pooled <- row(test(all_rows), all_rows)
    H <- if (is.null(hz)) 1L else length(unique(hz))
    lmn_warn_dm(method, pooled$clusters, H)
    if (H < 2L) return(pooled)
    rbind(pooled, do.call(rbind, lapply(sort(unique(hz)), function(h)
      row(test(hz == h), hz == h, as.integer(h)))))
  }))
  structure(do.call(rbind, col$value), class = c("ilm_compare", "data.frame"), checks = col$checks)
}

## The default Diebold-Mariano test for by = "origin".
lmn_dm_default <- "wpe"

## The usual test over-rejects with few origins over several
## horizons (studies/findings/compare_dm_fixed.md), so it warns there.
#' @keywords internal
#' @noRd
lmn_warn_dm <- function(method, Tn, H) {
  if (method == "dm" && H > 1L && Tn < 100L)
    lmn_signal("dm_few_origins", T = Tn, H = H)
  invisible(NULL)
}

## The Diebold-Mariano test over origins, by the chosen method.
#' @keywords internal
#' @noRd
lmn_se_origin <- function(d, origin, method, lag) {
  if (method == "dm") return(lmn_se_dm(d, origin, lag))
  if (is.null(origin))
    stop("`by = \"origin\"` needs an `origin` column.", call. = FALSE)
  dt <- as.numeric(tapply(d, factor(origin, levels = sort(unique(origin))), mean))
  Tn <- length(dt)
  if (Tn < 3L) stop("at least three origins are needed for the Diebold-Mariano test.",
                    call. = FALSE)
  w <- lmn_dm_wpe(dt)
  list(estimate = mean(dt), se = abs(mean(dt) / w$statistic), statistic = w$statistic,
       df = w$df, p = w$p, clusters = Tn,
       method = paste0("fixed-smoothing Diebold-Mariano over origins (m = ", w$m, ")"),
       test = "wpe", bandwidth = as.integer(w$m), lag = NA_integer_)
}

#' @keywords internal
#' @noRd
lmn_se_rows <- function(d) {
  n <- length(d)
  if (n < 2L) stop("at least two rows are needed for a standard error.", call. = FALSE)
  se <- stats::sd(d) / sqrt(n)
  lmn_t_result(mean(d), se, n - 1L, n, "independent rows", test = "independent")
}

## Cluster-robust (HC1) standard error of a mean: rows summed within clusters,
## with the G / (G - 1) small-sample factor, as sandwich::vcovCL() applies.
#' @keywords internal
#' @noRd
lmn_se_clustered <- function(d, cl) {
  n <- length(d); G <- length(unique(cl))
  if (G < 2L) stop("at least two clusters are needed for a standard error.",
                   call. = FALSE)
  u <- rowsum(d - mean(d), as.character(cl), reorder = FALSE)[, 1]
  se <- sqrt(G / (G - 1) * sum(u^2)) / n
  lmn_t_result(mean(d), se, G - 1L, G, "clustered", test = "clustered")
}

## The Diebold-Mariano test on the per-origin mean differences. The long-run
## variance is the sample autocovariances to `lag` with equal weights (Diebold
## and Mariano 1995). When it is negative, it is taken as 0 and the null
## rejected, as Diebold and Mariano (1995, Section 1.1) prescribe: the
## statistic is infinite, the p-value 0, and a warning is raised. Otherwise the
## statistic takes Harvey, Leybourne and Newbold's (1997) small-sample
## correction, with h = lag + 1, and is referred to t on T - 1 degrees of
## freedom.
#' @keywords internal
#' @noRd
lmn_se_dm <- function(d, origin, lag) {
  if (is.null(origin))
    stop("`by = \"origin\"` needs an `origin` column.", call. = FALSE)
  dt <- as.numeric(tapply(d, factor(origin, levels = sort(unique(origin))), mean))
  Tn <- length(dt)
  lag <- as.integer(lag)
  if (Tn < 3L) stop("at least three origins are needed for the Diebold-Mariano test.",
                    call. = FALSE)
  if (lag < 0L || lag >= Tn) stop("`lag` must be between 0 and the number of origins minus 1.",
                                  call. = FALSE)
  gam <- function(k) sum((dt[(k + 1L):Tn] - mean(dt)) * (dt[1L:(Tn - k)] - mean(dt))) / Tn
  v <- gam(0L) + 2 * sum(vapply(seq_len(lag), gam, 0))
  if (v < 0) {
    lmn_signal("dm_negative_variance", lag = lag)
    return(lmn_t_result(mean(dt), 0, Tn - 1L, Tn,
                        paste0("Diebold-Mariano over origins (lag ", lag,
                               ", variance taken as 0)"), test = "dm", lag = lag))
  }
  h <- lag + 1L
  se <- sqrt(v / Tn)
  hln <- sqrt((Tn + 1 - 2 * h + h * (h - 1) / Tn) / Tn)
  lmn_t_result(mean(dt), se / hln, Tn - 1L, Tn,
               paste0("Diebold-Mariano over origins (lag ", lag, ")"), test = "dm", lag = lag)
}

## Fixed-smoothing versions of the Diebold-Mariano test on the per-origin mean
## differences `dt` (Coroneo and Iacone 2020), for few origins. Each returns
## the statistic sqrt(T) mean / sigma-hat and its reference distribution.
##   WCE: the Bartlett-weighted autocovariances to bandwidth M = floor(T^(1/2)),
##        with Kiefer and Vogelsang's (2005) fixed-b critical values, through
##        the cubic in b = M / T that Coroneo and Iacone (2020, p. 5) give for
##        the 0.950 and 0.975 quantiles, so set levels only.
##   WPE: the Daniell-weighted periodogram at the first m = floor(T^(1/3))
##        Fourier frequencies, with the t distribution on 2m degrees of
##        freedom (Hualde and Iacone 2017), so any p-value.
#' @keywords internal
#' @noRd
lmn_dm_wce <- function(dt, M = floor(sqrt(length(dt)))) {
  Tn <- length(dt); e <- dt - mean(dt)
  gam <- function(k) sum(e[(k + 1L):Tn] * e[1L:(Tn - k)]) / Tn
  k <- seq_len(M - 1L)
  v <- gam(0L) + 2 * sum((1 - k / M) * vapply(k, gam, 0))
  b <- M / Tn
  ## two-sided 5% and 10% critical values (the 0.975 and 0.95 quantiles)
  cv <- c(`0.05` = 1.9600 + 2.9694 * b + 0.4160 * b^2 - 0.5324 * b^3,
          `0.10` = 1.6449 + 2.1859 * b + 0.3142 * b^2 - 0.3427 * b^3)
  list(statistic = sqrt(Tn) * mean(dt) / sqrt(v), M = M, b = b, cv = cv)
}

#' @keywords internal
#' @noRd
lmn_dm_wpe <- function(dt, m = floor(length(dt)^(1 / 3))) {
  Tn <- length(dt)
  lam <- 2 * pi * seq_len(m) / Tn
  tt <- seq_len(Tn)
  per <- vapply(lam, function(l) (sum(dt * cos(l * tt))^2 + sum(dt * sin(l * tt))^2) / Tn, 0)
  v <- mean(per)                                  # 2 pi times the averaged periodogram
  stat <- sqrt(Tn) * mean(dt) / sqrt(v)
  list(statistic = stat, m = m, df = 2L * m, p = 2 * stats::pt(-abs(stat), 2L * m))
}

#' @keywords internal
#' @noRd
lmn_default_lag <- function(horizon) {
  if (is.null(horizon)) 0L else max(0L, as.integer(max(horizon)) - 1L)
}

#' @keywords internal
#' @noRd
lmn_t_result <- function(est, se, df, clusters, method, test, bandwidth = NA_integer_,
                         lag = NA_integer_) {
  stat <- est / se
  list(estimate = est, se = se, statistic = stat, df = df,
       p = 2 * stats::pt(-abs(stat), df), clusters = clusters, method = method,
       test = test, bandwidth = as.integer(bandwidth), lag = as.integer(lag))
}

#' @export
print.ilm_compare <- function(x, ...) {
  for (k in seq_len(nrow(x))) {
    r <- x[k, ]
    dir <- if (r$difference < 0) "lower" else "higher"
    how <- switch(r$by, row = "rows treated as independent",
                  origin = paste0(sub(" over origins", "", r$method), " over ", lmn_fmt_count(r$clusters),
                                  " origins", if (!is.na(r$horizon)) paste0(", horizon ", lmn_fmt_data(r$horizon))
                                  else if (sum(x$model == r$model & x$rule == r$rule) > 1L)
                                    ", pooled over horizons" else ""),
                  paste0("clustered by ", r$by, ", ", lmn_fmt_count(r$clusters), " clusters"))
    cat(sprintf("%s's %s score is %s than %s's by %s a row (SE %s, %s; %s).\n",
                r$model, r$rule, dir, r$reference, lmn_fmt_sig(abs(r$difference)),
                lmn_fmt_sig(r$se), how, lmn_fmt_p(r$p_value)))
    if (r$dropped > 0)
      cat("  ", lmn_fmt_count(r$dropped), " row(s) with a missing or infinite score were left out.\n",
          sep = "")
  }
  cat("Scores are losses: lower is better.\n")
  invisible(x)
}
