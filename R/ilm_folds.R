#' Cross-validation folds that follow the data's structure
#'
#' Splits the rows of `data` into training and test sets for cross-validation.
#' The folds are row indices, so they can be inspected, and reused for every
#' model being compared: two models scored on the same folds can be compared
#' fold by fold.
#'
#' @section Which kind of fold:
#' The folds decide which question the cross-validation answers, so they should
#' match the prediction being made.
#'
#' * `"random"`: rows are held out at random. Every held-out row's group is
#'   still in the training data, so this scores predictions for **units the
#'   model has seen**. On data with a time structure it scores interpolation,
#'   not forecasting.
#' * `"grouped"`: whole groups are held out together, so this scores
#'   predictions for **new units**. Random folds would overstate how well a
#'   model predicts a unit it has never seen.
#' * `"blocked"`: contiguous blocks of time are held out, optionally with a gap
#'   either side so that rows next to the block do not leak into training.
#' * `"rolling"`: rolling-origin evaluation, also called time-series
#'   cross-validation. Each fold trains on the times up to an origin and tests
#'   on the `horizon` times after it, so every test row is a genuine forecast.
#'   For a model with an autoregressive term, this is the scheme that scores
#'   forecasts. Random and blocked folds let the model see both sides of a
#'   held-out time, so they score interpolation and overstate forecast skill.
#'
#' Rolling origins are the reference's "evaluation on a rolling forecasting
#' origin": every training set holds only rows before its test rows, and the
#' scores are averaged over the test sets. Which question each of the other
#' kinds answers is limn's own reasoning.
#'
#' @section Time points and horizons:
#' Blocked and rolling folds count **time points**: the distinct values of
#' `time` in `data`, in order. Every row at a time point goes into the same set,
#' so a panel of several series is split at the same times for all of them. A
#' test row's horizon is how many time points after the origin it lies, so with
#' `gap = 0` the first test time point has horizon 1.
#'
#' @section Ex-ante and ex-post:
#' A forecast is *ex-ante* when it uses only what would have been known at the
#' origin. Covariates the model uses at their realised future values, such as
#' the weather on the day being forecast, make it *ex-post*, and flatter the
#' forecast. Name in `known` the covariates whose future values would genuinely
#' be known, such as a calendar or a planned dose. The folds record it so the
#' evaluation can say which kind it is.
#'
#' @param data A data frame: the rows the models are fitted to.
#' @param type The kind of fold; see *Which kind of fold*.
#' @param K Number of folds, for random, grouped and blocked folds.
#' @param group For grouped folds, the name of the grouping column in `data`.
#' @param time For blocked and rolling folds, the name of the time column in
#'   `data`: numeric, a `Date` or a date-time.
#' @param initial For rolling folds, the number of time points in the first
#'   training window.
#' @param horizon For rolling folds, the number of time points in each test set.
#' @param step For rolling folds, how many time points the origin moves from one
#'   fold to the next.
#' @param window For rolling folds: `"expanding"` trains on every time point up
#'   to the origin, `"sliding"` on the last `size` of them.
#' @param size For a sliding window, the number of time points it trains on.
#'   Defaults to `initial`.
#' @param gap Time points left out between training and test rows: after the
#'   origin for rolling folds, and on both sides of the test block for blocked
#'   folds.
#' @param known For rolling folds, the names of the covariates whose future
#'   values would be known at the origin; see *Ex-ante and ex-post*.
#' @param seed Seed for random and grouped folds. The caller's random-number
#'   state is restored afterwards.
#' @return An object of class `"ilm_folds"`: a list with one element per fold,
#'   each holding `train` and `test` row indices, and for rolling folds the
#'   fold's `origin` and each test row's `horizon`. `as.data.frame()` gives one
#'   row per data row and fold.
#' @references Hyndman, R. J. and Athanasopoulos, G. (2021). *Forecasting:
#'   Principles and Practice*, 3rd edition. OTexts. Section 5.10, "Time series
#'   cross-validation" (online edition,
#'   <https://otexts.com/fpp3/toolbox.html#sec-tscv>, read 27 October 2026).
#' @seealso [ilm_cv()], which scores a model on the folds; [ilm_compare()],
#'   which compares two models' scores on the same folds.
#' @examples
#' d <- data.frame(site = rep(letters[1:6], each = 8), month = rep(1:8, 6),
#'                 y = rnorm(48))
#' ilm_folds(d, "grouped", K = 3, group = "site")
#' f <- ilm_folds(d, "rolling", time = "month", initial = 5, horizon = 2)
#' f
#' head(as.data.frame(f))
#' @export
ilm_folds <- function(data, type = c("random", "grouped", "blocked", "rolling"),
                      K = 5L, group = NULL, time = NULL, initial = NULL,
                      horizon = 1L, step = 1L,
                      window = c("expanding", "sliding"), size = NULL, gap = 0L,
                      known = NULL, seed = 1L) {
  if (!is.data.frame(data))
    stop("`data` must be a data frame; it is ", class(data)[1], ".", call. = FALSE)
  type <- match.arg(type)
  window <- match.arg(window)
  n <- nrow(data)
  if (n < 2L) stop("`data` needs at least two rows to split.", call. = FALSE)
  gap <- lmn_check_count(gap, "gap", min = 0L)

  folds <- switch(type,
    random  = lmn_folds_random(n, lmn_check_count(K, "K", 2L), seed),
    grouped = lmn_folds_grouped(lmn_check_column(data, group, "group"),
                                lmn_check_count(K, "K", 2L), seed),
    blocked = lmn_folds_blocked(lmn_check_column(data, time, "time"),
                                lmn_check_count(K, "K", 2L), gap),
    rolling = lmn_folds_rolling(lmn_check_column(data, time, "time"), initial,
                                lmn_check_count(horizon, "horizon"),
                                lmn_check_count(step, "step"), window, size, gap))
  if (!is.null(known) && !is.character(known))
    stop("`known` must be the names of covariates, as strings.", call. = FALSE)

  structure(folds, class = "ilm_folds", type = type, n = n,
            K = length(folds), group = group, time = time,
            initial = if (type == "rolling") as.integer(initial) else NULL,
            horizon = if (type == "rolling") as.integer(horizon) else NULL,
            step = if (type == "rolling") as.integer(step) else NULL,
            window = if (type == "rolling") window else NULL,
            size = if (type == "rolling" && window == "sliding")
              as.integer(if (is.null(size)) initial else size) else NULL,
            gap = gap, known = known,
            seed = if (type %in% c("random", "grouped")) seed else NULL)
}

## Rows shuffled into K folds whose sizes differ by at most one.
#' @keywords internal
#' @noRd
lmn_folds_random <- function(n, K, seed) {
  if (K > n) stop("`K` = ", K, " is more folds than the ", n, " rows.", call. = FALSE)
  id <- lmn_with_seed(seed, sample(rep_len(seq_len(K), n)))
  lapply(seq_len(K), function(k) list(train = which(id != k), test = which(id == k)))
}

## Groups shuffled into K folds whose group counts differ by at most one.
#' @keywords internal
#' @noRd
lmn_folds_grouped <- function(g, K, seed) {
  lev <- unique(as.character(g))
  G <- length(lev)
  if (K > G)
    stop("`K` = ", K, " is more folds than the ", G, " groups; with grouped ",
         "folds each fold holds out whole groups.", call. = FALSE)
  gid <- lmn_with_seed(seed, sample(rep_len(seq_len(K), G)))
  id <- gid[match(as.character(g), lev)]
  lapply(seq_len(K), function(k) list(train = which(id != k), test = which(id == k)))
}

## K contiguous blocks of time points, each held out in turn; `gap` time points
## either side of the block are left out of training.
#' @keywords internal
#' @noRd
lmn_folds_blocked <- function(tm, K, gap) {
  ut <- sort(unique(tm)); nt <- length(ut)
  if (K > nt)
    stop("`K` = ", K, " is more blocks than the ", nt, " time points.", call. = FALSE)
  pos <- match(tm, ut)
  blk <- ceiling(seq_len(nt) * K / nt)          # block of each time point
  lapply(seq_len(K), function(k) {
    tp <- which(blk == k)
    lo <- min(tp) - gap; hi <- max(tp) + gap
    list(train = which(pos < lo | pos > hi), test = which(blk[pos] == k))
  })
}

## Rolling origins. The origin is a position among the sorted time points; each
## fold trains up to it (all of it, or the last `size`) and tests on the
## `horizon` time points starting `gap` + 1 after it. Only full test sets are
## made, so every origin is scored over the same horizons.
#' @keywords internal
#' @noRd
lmn_folds_rolling <- function(tm, initial, horizon, step, window, size, gap) {
  initial <- lmn_check_count(initial, "initial")
  ut <- sort(unique(tm)); nt <- length(ut)
  if (window == "sliding") {
    size <- if (is.null(size)) initial else lmn_check_count(size, "size")
    if (size > initial)
      stop("`size` = ", size, " is longer than the first training window (`initial` = ",
           initial, ").", call. = FALSE)
  }
  last <- nt - gap - horizon                     # the last origin with a full test set
  if (initial > last)
    stop("with ", nt, " time points there is no room for a training window of ",
         initial, ", a gap of ", gap, " and a horizon of ", horizon, ".", call. = FALSE)
  origins <- seq.int(initial, last, by = step)
  pos <- match(tm, ut)
  lapply(origins, function(o) {
    first <- if (window == "sliding") o - size + 1L else 1L
    test <- which(pos > o + gap & pos <= o + gap + horizon)
    list(train = which(pos >= first & pos <= o), test = test,
         origin = ut[o], horizon = pos[test] - o)
  })
}

#' @export
print.ilm_folds <- function(x, ...) {
  a <- attributes(x)
  ntr <- vapply(x, function(f) length(f$train), 0L)
  nte <- vapply(x, function(f) length(f$test), 0L)
  ## rows, horizons and time points are counts: whole, grouped from 1,000
  rng <- function(v) if (min(v) == max(v)) lmn_fmt_count(min(v)) else
    paste(lmn_fmt_count(min(v)), "to", lmn_fmt_count(max(v)))
  label <- switch(a$type,
    random = "random folds of rows",
    grouped = paste0("grouped folds (whole groups of '", a$group, "' held out)"),
    blocked = paste0("blocked folds of time ('", a$time, "')"),
    rolling = paste0("rolling origins over time ('", a$time, "'), ", a$window,
                     " window"))
  cat("<ilm_folds> ", a$K, " ", label, "\n", sep = "")
  cat("  ", lmn_fmt_count(a$n), " rows; training rows ", rng(ntr), ", test rows ", rng(nte),
      " per fold\n", sep = "")
  if (a$type == "rolling") {
    cat("  first origin after ", a$initial, " time points, then every ", a$step,
        if (!is.null(a$size)) paste0("; training on the last ", a$size) else "",
        "\n", sep = "")
    h <- unlist(lapply(x, `[[`, "horizon"))
    cat("  horizons ", rng(h), if (a$gap > 0) paste0(" (a gap of ", a$gap, " first)")
        else "", "\n", sep = "")
    cat(if (length(a$known))
      paste0("  known in advance: ", paste(a$known, collapse = ", "),
             "; any other covariate makes the evaluation ex-post\n")
      else "  known in advance: none named, so any covariate makes the evaluation ex-post\n")
  } else if (a$type == "blocked" && a$gap > 0) {
    cat("  ", a$gap, " time points left out either side of each block\n", sep = "")
  }
  invisible(x)
}

#' @export
as.data.frame.ilm_folds <- function(x, row.names = NULL, optional = FALSE, ...) {
  rolling <- identical(attr(x, "type"), "rolling")
  out <- lapply(seq_along(x), function(k) {
    f <- x[[k]]
    d <- data.frame(fold = k, row = c(f$train, f$test),
                    set = rep(c("train", "test"), c(length(f$train), length(f$test))))
    if (rolling) {
      d$origin <- rep(f$origin, nrow(d))
      d$horizon <- c(rep(NA_integer_, length(f$train)), f$horizon)
    }
    d
  })
  do.call(rbind, out)
}
