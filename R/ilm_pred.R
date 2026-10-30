## ---- the predictive distribution --------------------------------------------
##
## An "ilm_pred" holds, for each of n rows, a predictive distribution that is a
## mixture over S draws of the fitted parameters:
##
##     F_i(y) = sum_s w_s F(y | lp[[.]][i, s], par[[.]][s])
##
## `dist` supplies F's d, p, q, r and mean in the family's own
## parameterisation, in the shape illume's ilm_dist() (export E4) returns.
## Every function takes `lp`, a named list of link-scale linear predictors, and
## `par`, a named list of the family's other parameters on illume's internal
## scale. limn passes them through and never transforms them itself, so it
## cannot drift from the fit's own parameterisation.
##
## Scores, quantiles and PITs are all read from the mixture: averaging the
## family's CDF or density over the draws, rather than simulating a response
## from each, removes the simulation's noise (prototype P4a,
## studies/prototypes/p4_scoring_integration.R).

#' Build a predictive distribution from per-draw parameters
#'
#' @param dist The family's distribution functions, in the shape of
#'   `illume::ilm_dist()`: a list with functions `d(y, lp, par, log)`,
#'   `p(q, lp, par)`, `q(prob, lp, par)`, `r(lp, par)` and `mean(lp, par)`, a
#'   logical `discrete`, and character vectors `lp_names` and `par_names`.
#' @param lp Named list of n x S matrices: each row's linear predictors at each
#'   draw.
#' @param par Named list of the family's other parameters: a length-S vector
#'   (one value per draw) or an n x S matrix.
#' @param weights Draw weights; equal by default.
#' @param rows Data frame describing the n rows. Columns such as `group`,
#'   `fold`, `origin` and `horizon` are carried into scores, where
#'   [ilm_compare()] can cluster by them.
#' @param meta Named list recording how the prediction was made.
#' @keywords internal
#' @noRd
lmn_new_pred <- function(dist, lp, par = list(), weights = NULL, rows = NULL,
                         meta = list(), pvec = list()) {
  lmn_check_dist(dist)
  if (!is.list(lp) || !length(lp) || is.null(names(lp)))
    stop("`lp` must be a named list of matrices, one per linear predictor.",
         call. = FALSE)
  lp <- lapply(lp, function(m) { m <- as.matrix(m); storage.mode(m) <- "double"; m })
  dm <- dim(lp[[1]])
  if (!all(vapply(lp, function(m) identical(dim(m), dm), TRUE)))
    stop("every linear predictor in `lp` must have the same rows and draws.",
         call. = FALSE)
  miss <- setdiff(dist$lp_names, names(lp))
  if (length(miss))
    stop("the family needs linear predictor(s) ", paste(miss, collapse = ", "),
         ", which `lp` does not have.", call. = FALSE)
  n <- dm[1]; S <- dm[2]
  par <- lapply(par, function(v) {
    if (is.matrix(v)) {
      if (!identical(dim(v), dm))
        stop("a matrix in `par` must be rows x draws, like `lp`.", call. = FALSE)
      return(v)
    }
    if (length(v) == 1L) v <- rep(v, S)
    if (length(v) != S)
      stop("a vector in `par` must have one value per draw (", S, ").", call. = FALSE)
    v
  })
  ## a parameter that is a vector in every draw, such as an ordinal response's
  ## thresholds, is a k x S matrix passed whole, draw by draw
  pvec <- lapply(pvec, function(v) {
    v <- as.matrix(v)
    if (ncol(v) != S)
      stop("a matrix in `pvec` must have one column per draw (", S, ").", call. = FALSE)
    v
  })
  missp <- setdiff(dist$par_names, c(names(par), names(pvec)))
  if (length(missp))
    stop("the family needs parameter(s) ", paste(missp, collapse = ", "),
         ", which `par` does not have.", call. = FALSE)
  w <- if (is.null(weights)) rep(1 / S, S) else {
    if (length(weights) != S || any(!is.finite(weights)) || any(weights < 0) ||
        sum(weights) <= 0)
      stop("`weights` must be ", S, " non-negative numbers, not all zero.",
           call. = FALSE)
    weights / sum(weights)
  }
  if (is.null(rows)) rows <- data.frame(row = seq_len(n))
  if (!is.data.frame(rows) || nrow(rows) != n)
    stop("`rows` must be a data frame with one row per predicted row (", n, ").",
         call. = FALSE)
  structure(list(kind = "response", dist = dist, lp = lp, par = par, pvec = pvec,
                 weights = w, rows = rows, meta = meta),
            class = "ilm_pred")
}

## A prediction of a mean or a linear predictor rather than of a response: its
## distribution is the parameter draws' values themselves, an n x S matrix, and
## its point estimate is the value at the estimates, as predict() gives it.
#' @keywords internal
#' @noRd
lmn_new_draw_pred <- function(draws, estimate, kind = c("mean", "link"),
                              rows = NULL, meta = list()) {
  kind <- match.arg(kind)
  draws <- as.matrix(draws); storage.mode(draws) <- "double"
  n <- nrow(draws)
  if (length(estimate) != n)
    stop("`estimate` must have one value per row (", n, ").", call. = FALSE)
  if (is.null(rows)) rows <- data.frame(row = seq_len(n))
  structure(list(kind = kind, draws = draws, estimate = as.numeric(estimate),
                 weights = rep(1 / ncol(draws), ncol(draws)), rows = rows, meta = meta),
            class = "ilm_pred")
}

## A predictive distribution of an observable quantity held as samples: a
## total of new observations, say, summed from outcomes simulated from each
## draw. Its CDF, quantiles and mean are the samples' own, and it is scored
## from them. `discrete` says whether its values are whole numbers.
#' @keywords internal
#' @noRd
lmn_new_sample_pred <- function(samples, discrete, rows = NULL, meta = list(),
                                family = NULL, what = "simulated outcomes") {
  samples <- as.matrix(samples); storage.mode(samples) <- "double"
  n <- nrow(samples); M <- ncol(samples)
  if (is.null(rows)) rows <- data.frame(row = seq_len(n))
  est <- apply(samples, 1L, stats::quantile, probs = 0.5, names = FALSE,
               type = if (discrete) 1L else 7L)
  structure(list(kind = "sample", draws = samples, estimate = est, weights = rep(1 / M, M),
                 dist = list(family = family, discrete = discrete, categories = NULL, zero = NULL),
                 what = what, rows = rows, meta = meta),
            class = "ilm_pred")
}

## Is this prediction draws of a mean or a linear predictor (no outcome to
## score it against)?
#' @keywords internal
#' @noRd
lmn_is_draws <- function(pred) !is.null(pred$kind) && pred$kind %in% c("mean", "link")

## Is it held as a matrix of values whose empirical distribution it is: draws
## of a mean, or samples of an outcome?
#' @keywords internal
#' @noRd
lmn_is_empirical <- function(pred) !is.null(pred$kind) && pred$kind %in% c("mean", "link", "sample")

#' @keywords internal
#' @noRd
lmn_pred_n <- function(pred) if (lmn_is_empirical(pred)) nrow(pred$draws) else nrow(pred$lp[[1]])

#' @keywords internal
#' @noRd
lmn_check_dist <- function(dist) {
  fns <- c("d", "p", "q", "r", "mean")
  ok <- is.list(dist) && all(vapply(fns, function(f) is.function(dist[[f]]), TRUE)) &&
    is.logical(dist$discrete) && length(dist$discrete) == 1L &&
    is.character(dist$lp_names)
  if (!ok)
    stop("`dist` must supply functions d, p, q, r and mean, a logical `discrete` ",
         "and `lp_names`, as illume's ilm_dist() does.", call. = FALSE)
  invisible(dist)
}

#' @export
print.ilm_pred <- function(x, ...) {
  n <- lmn_pred_n(x); S <- length(x$weights)
  what <- if (!is.null(x$what)) x$what else switch(if (is.null(x$kind)) "response" else x$kind,
                 mean = "the mean", link = "the linear predictor", "a new observation")
  unit <- if (identical(x$kind, "sample")) " sample" else " draw"
  cat("<ilm_pred> predictive distributions of ", what, " for ", lmn_fmt_count(n), " row",
      if (n != 1L) "s", ", from ", lmn_fmt_count(S), unit, if (S != 1L) "s", "\n", sep = "")
  if (!is.null(x$dist$family)) cat("  family: ", x$dist$family, "\n", sep = "")
  ## the settings: counts (draws asked for, covariate draws) whole and grouped,
  ## the rest as given
  shown <- function(v) if (is.numeric(v) && all(is.finite(v) & v == round(v))) lmn_fmt_count(v) else format(v)
  for (k in names(x$meta))
    cat("  ", k, ": ", paste(shown(x$meta[[k]]), collapse = ", "), "\n", sep = "")
  nl <- sum(x$rows$new_level %in% TRUE)
  if (nl) cat("  note: ", lmn_diag_text(lmn_diag("unseen_levels", n = nl)), "\n", sep = "")
  invisible(x)
}

## ---- evaluating the mixture ---------------------------------------------------

## Rows are processed in chunks, so that one call to a family function sees at
## most this many values.
lmn_chunk <- 4e6

## The linear predictors and parameters of rows `i` at draw `s`, each repeated
## to length length(i) * L, to go with an evaluation matrix of L points a row.
#' @keywords internal
#' @noRd
lmn_draw_args <- function(pred, i, s, L = 1L) {
  lp <- lapply(pred$lp, function(M) rep(M[i, s], L))
  par <- lapply(pred$par, function(v)
    if (is.matrix(v)) rep(v[i, s], L) else rep(v[s], length(i) * L))
  for (nm in names(pred$pvec)) par[[nm]] <- pred$pvec[[nm]][, s]
  list(lp = lp, par = par)
}

## Draws are evaluated in blocks: one family call covers a block of draws at
## once, with at most this many values, so that small evaluations (a few rows,
## a few points) are not dominated by one call per draw.
lmn_block <- 2e6

## The draws of block `ss`, laid out to go with rep(as.vector(X), length(ss))
## for an m x L matrix X: rows fastest, then points, then draws.
#' @keywords internal
#' @noRd
lmn_block_args <- function(pred, i, ss, L = 1L) {
  cols <- rep(ss, each = L); n <- length(i) * L
  lp <- lapply(pred$lp, function(M) as.vector(M[i, cols, drop = FALSE]))
  par <- lapply(pred$par, function(v)
    if (is.matrix(v)) as.vector(v[i, cols, drop = FALSE]) else rep(v[ss], each = n))
  list(lp = lp, par = par)
}

## Blocks of draws for an evaluation of `n` values a draw. A family with
## vector parameters per draw (pvec, such as ordinal cutpoints) takes one draw
## at a time.
#' @keywords internal
#' @noRd
lmn_draw_blocks <- function(pred, n) {
  S <- length(pred$weights)
  size <- if (length(pred$pvec)) 1L else max(1L, floor(lmn_block / max(1L, n)))
  split(seq_len(S), ceiling(seq_len(S) / size))
}

## The mixture's CDF ("p") or density ("d") at an m x L matrix X of points for
## rows i: sum_s w_s f(X | draw s), a block of draws per family call.
#' @keywords internal
#' @noRd
lmn_mix <- function(pred, i, X, what = c("p", "d", "s")) {
  what <- match.arg(what)
  X <- as.matrix(X)
  ## "s", the survival function, from the family's own upper tail, never 1 - F
  f <- if (what == "s") function(q, lp, par) pred$dist$p(q, lp, par, lower.tail = FALSE) else pred$dist[[what]]
  out <- matrix(0, nrow(X), ncol(X))
  x <- as.vector(X)
  for (ss in lmn_draw_blocks(pred, length(x))) {
    if (length(pred$pvec)) {
      a <- lmn_draw_args(pred, i, ss, ncol(X))
      out <- out + pred$weights[ss] * matrix(f(x, a$lp, a$par), nrow(X))
      next
    }
    a <- lmn_block_args(pred, i, ss, ncol(X))
    v <- matrix(f(rep(x, length(ss)), a$lp, a$par), length(x))
    out <- out + matrix(v %*% pred$weights[ss], nrow(X))
  }
  out
}

## Each draw's own quantile at probability `prob`, for rows i: an m x S matrix.
#' @keywords internal
#' @noRd
lmn_draw_quantiles <- function(pred, i, prob) {
  S <- length(pred$weights); m <- length(i)
  out <- matrix(NA_real_, m, S)
  for (ss in lmn_draw_blocks(pred, m)) {
    if (length(pred$pvec)) {
      a <- lmn_draw_args(pred, i, ss)
      out[, ss] <- pred$dist$q(rep(prob, m), a$lp, a$par)
      next
    }
    a <- lmn_block_args(pred, i, ss)
    out[, ss] <- pred$dist$q(rep(prob, m * length(ss)), a$lp, a$par)
  }
  out
}

## Row chunks for an evaluation of L points a row.
#' @keywords internal
#' @noRd
lmn_row_chunks <- function(n, L = 1L) {
  size <- max(1L, floor(lmn_chunk / max(1L, L)))
  split(seq_len(n), ceiling(seq_len(n) / size))
}

## An average over rows is held as the total's distribution with a `scale`
## (one over the number of rows) that turns the total into the average. The
## distribution functions below report the average; scores and PITs work on
## the total's own values, so a count's average keeps its discreteness.
#' @keywords internal
#' @noRd
lmn_rescale <- function(pred, v) if (is.null(pred$scale)) v else v * pred$scale

## The mixture's mean for every row.
#' @keywords internal
#' @noRd
lmn_pred_mean <- function(pred) {
  if (lmn_is_empirical(pred)) return(lmn_rescale(pred, as.numeric(pred$draws %*% pred$weights)))
  if (!is.null(pred$dist$categories))
    stop("an ordered outcome has no numeric mean: its prediction is a probability ",
         "for each category.", call. = FALSE)
  n <- nrow(pred$lp[[1]])
  out <- numeric(n)
  for (i in lmn_row_chunks(n)) {
    acc <- numeric(length(i))
    for (ss in lmn_draw_blocks(pred, length(i))) {
      if (length(pred$pvec)) {
        a <- lmn_draw_args(pred, i, ss)
        acc <- acc + pred$weights[ss] * pred$dist$mean(a$lp, a$par)
        next
      }
      a <- lmn_block_args(pred, i, ss)
      acc <- acc + as.vector(matrix(pred$dist$mean(a$lp, a$par), length(i)) %*% pred$weights[ss])
    }
    out[i] <- acc
  }
  lmn_rescale(pred, out)
}

## The mixture's CDF at one point per row.
#' @keywords internal
#' @noRd
lmn_pred_cdf <- function(pred, q) {
  if (!is.null(pred$scale)) {                          # the average's q is the total's q / scale,
    q <- q / pred$scale                                # on the total's whole numbers when discrete
    if (isTRUE(pred$dist$discrete)) q <- floor(q + 1e-8)
  }
  if (lmn_is_empirical(pred)) return(as.numeric((pred$draws <= q) %*% pred$weights))
  n <- nrow(pred$lp[[1]])
  out <- numeric(n)
  for (i in lmn_row_chunks(n)) out[i] <- lmn_mix(pred, i, matrix(q[i]), "p")[, 1]
  out
}

## The mixture's quantile at probability `prob` for every row. The mixture's
## quantile lies between the smallest and largest of the draws' own quantiles,
## so it is found by bisection inside that bracket: on the integers for a
## discrete family, and to a relative width of 1e-12 otherwise.
#' @keywords internal
#' @noRd
lmn_pred_quantile <- function(pred, prob) {
  ## draws of a mean, or samples: their empirical quantile. For whole-number
  ## samples it is the smallest value k with F(k) >= prob, as for a discrete
  ## mixture.
  if (lmn_is_empirical(pred))
    return(lmn_rescale(pred, apply(pred$draws, 1L, stats::quantile, probs = prob, names = FALSE,
                                   type = if (isTRUE(pred$dist$discrete)) 1L else 7L)))
  n <- nrow(pred$lp[[1]])
  out <- numeric(n)
  for (i in lmn_row_chunks(n)) {
    Q <- lmn_draw_quantiles(pred, i, prob)
    lo <- apply(Q, 1L, min); hi <- apply(Q, 1L, max)
    if (pred$dist$discrete) {
      out[i] <- lmn_mix_qint(pred, i, prob, lo, hi)
    } else {
      b <- lmn_widen_bracket(pred, i, prob, lo, hi, FALSE)
      out[i] <- lmn_quantile_root(pred, i, prob, b$lo, b$hi, b$flo, b$fhi)
    }
  }
  lmn_rescale(pred, out)
}

## A quantile's bracket from the draws' own quantiles holds the mixture's
## quantile only if those quantiles are right, and a family's quantile function
## can fail (R's qbeta at tiny shape parameters is off by up to 0.5 in
## probability). So the bracket is checked against the mixture's CDF before the
## bisection: F(lo) <= p (for whole numbers F(lo - 1) < p) and F(hi) >= p. An
## end that fails is moved out, doubling the step, until it holds. A row where
## it cannot be made to hold (a CDF that is not a number) is NA, with a message.
#' @keywords internal
#' @noRd
lmn_widen_bracket <- function(pred, i, prob, lo, hi, discrete) {
  Fx <- function(rows, x) lmn_mix(pred, rows, matrix(x), "p")[, 1]
  ## F at the ends as checked: F(lo) (F(lo - 1) for whole numbers) and F(hi)
  flo <- Fx(i, if (discrete) lo - 1 else lo); fhi <- Fx(i, hi)
  lo_bad <- function(v) { b <- !(if (discrete) v < prob else v <= prob); b[is.na(b)] <- TRUE; b }
  hi_bad <- function(v) { b <- !(v >= prob); b[is.na(b)] <- TRUE; b }
  bad_lo <- lo_bad(flo); bad_hi <- hi_bad(fhi)
  step <- pmax(hi - lo, 1e-8 * pmax(1, abs(lo), abs(hi)), if (discrete) 1 else 0)
  for (it in seq_len(1100L)) {
    if (!any(bad_lo | bad_hi)) break
    if (any(bad_lo)) {
      k <- which(bad_lo)
      lo[k] <- lo[k] - step[k]
      if (discrete) lo[k] <- floor(lo[k])
      flo[k] <- Fx(i[k], if (discrete) lo[k] - 1 else lo[k]); bad_lo[k] <- lo_bad(flo[k])
    }
    if (any(bad_hi)) {
      k <- which(bad_hi)
      hi[k] <- hi[k] + step[k]
      if (discrete) hi[k] <- ceiling(hi[k])
      fhi[k] <- Fx(i[k], hi[k]); bad_hi[k] <- hi_bad(fhi[k])
    }
    step <- 2 * step
  }
  stuck <- bad_lo | bad_hi | !is.finite(lo) | !is.finite(hi)
  if (any(stuck)) {
    message("a quantile could not be bracketed for ", sum(stuck), " row(s): ",
            "their predictive CDF is not a number there, so the quantile is NA.")
    lo[stuck] <- NA_real_; hi[stuck] <- NA_real_
  }
  list(lo = lo, hi = hi, flo = flo, fhi = fhi)
}

## The mixture's quantile for rows i, continuous: the bracket [lo, hi], with
## F(lo) <= p <= F(hi), closed to the bisection's own stopping width,
## 1e-12 max(1, |lo|, |hi|), and its midpoint returned, so the answer is
## inf{x: F(x) >= p} to that tolerance, as the bisection gave it
## (studies/designs/quantile_b1.md, with its addendum). The next point is the
## Illinois secant, the stale end's value halved; a point within a tolerance of
## an end moves a tolerance inside; three steps that do not halve the bracket
## are followed by a bisection. The cap: a row takes a secant step only while,
## were that step to gain nothing, bisection could still finish within its
## budget, the bisection's own count plus 2; otherwise it
## bisects that pass, and the guard is asked again on the next.
##
## The budget: bisection from the starting bracket, of width w0, stops at a
## width no more than the stopping width of any bracket that holds the
## quantile, and every bracket here holds it. So from the current bracket,
## ceiling(log2(w0 / (1e-12 max(1, |lo|, |hi|)))) is a lower bound on the
## bisection's own count (less 1e-3 for the midpoints' rounding), and grows as
## the bracket closes. The guard asks how many passes bisection could need
## from where the row stands, counting from the smallest stopping width the
## bracket can still reach, from its point nearest zero. So a row never takes
## more than the bisection's passes plus 2 (the proof is in
## studies/designs/quantile_b1_retry.md).
#' @keywords internal
#' @noRd
lmn_quantile_root <- function(pred, i, prob, lo, hi, flo, fhi) {
  n <- length(i)
  tol <- function(k) 1e-12 * pmax(1, abs(lo[k]), abs(hi[k]))
  ## the smallest stopping width reachable inside [lo, hi], and the passes
  ## bisection needs to reach it; the 1e-3 covers the midpoints' rounding
  tol_floor <- function(k) 1e-12 * pmax(1, ifelse(lo[k] <= 0 & hi[k] >= 0, 0, pmin(abs(lo[k]), abs(hi[k]))))
  need_safe <- function(k) pmax(0, ceiling(log2((hi[k] - lo[k]) / tol_floor(k)) + 1e-3))
  ## the budget: a lower bound on the bisection's count, plus 2
  w_start <- hi - lo
  budget <- function(k) pmax(0, ceiling(log2(w_start[k] / tol(k)) - 1e-3)) + 2
  glo <- flo - prob; ghi <- fhi - prob
  all_k <- seq_len(n)
  open <- !is.na(lo) & (hi - lo) > tol(all_k)
  used <- integer(n); side <- integer(n); stall <- integer(n); w0 <- hi - lo
  for (it in seq_len(400L)) {
    k <- which(open)
    if (!length(k)) break
    w <- hi[k] - lo[k]; tk <- tol(k)
    bisect_now <- used[k] + 1L + need_safe(k) > budget(k)
    x <- lo[k] - glo[k] * w / (ghi[k] - glo[k])
    mid <- bisect_now | stall[k] >= 3L | !is.finite(x) | x <= lo[k] | x >= hi[k]
    x[mid] <- (lo[k][mid] + hi[k][mid]) / 2
    nl <- !mid & x - lo[k] < tk; x[nl] <- pmin(lo[k][nl] + tk[nl], (lo[k][nl] + hi[k][nl]) / 2)
    nh <- !mid & hi[k] - x < tk; x[nh] <- pmax(hi[k][nh] - tk[nh], (lo[k][nh] + hi[k][nh]) / 2)
    g <- lmn_mix(pred, i[k], matrix(x), "p")[, 1] - prob
    used[k] <- used[k] + 1L
    up <- !is.na(g) & g >= 0
    ku <- k[up]; kd <- k[!up]
    hi[ku] <- x[up]; ghi[ku] <- g[up]
    rh <- ku[side[ku] == 1L & !mid[up]]; glo[rh] <- glo[rh] / 2
    side[ku] <- 1L
    lo[kd] <- x[!up]; glo[kd] <- g[!up]
    rl <- kd[side[kd] == -1L & !mid[!up]]; ghi[rl] <- ghi[rl] / 2
    side[kd] <- -1L
    wn <- hi[k] - lo[k]
    reset <- wn <= w0[k] / 2 | mid
    stall[k] <- ifelse(reset, 0L, stall[k] + 1L)
    w0[k][reset] <- wn[reset]
    open[k] <- (hi[k] - lo[k]) > tol(k)
  }
  (lo + hi) / 2
}

## A discrete mixture's quantile for rows i: the smallest integer k in
## [lo, hi] with F(k) >= prob, by bisection. lo and hi bracket it, as the
## smallest and largest of the draws' own quantiles do.
#' @keywords internal
#' @noRd
lmn_mix_qint <- function(pred, i, prob, lo, hi) {
  b <- lmn_widen_bracket(pred, i, prob, lo, hi, TRUE); lo <- b$lo; hi <- b$hi
  repeat {
    open <- !is.na(lo) & hi > lo
    if (!any(open)) break
    mid <- floor((lo + hi) / 2)
    Fm <- lmn_mix(pred, i[open], matrix(mid[open]), "p")[, 1]
    up <- Fm >= prob
    hi[open][up] <- mid[open][up]
    lo[open][!up] <- mid[open][!up] + 1
  }
  lo
}
