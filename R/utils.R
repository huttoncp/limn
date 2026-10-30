## ---- numbers in prose, as illume's prose writes them ---------------------------
##
## Every number limn prints comes through these, and through ilm_disp() in
## shared-helpers.R, the same file in each of the family's packages. Each
## number declares what it is:
##   an estimate (a score, a difference, its standard error, a width, a test
##     statistic) to three significant figures with its trailing zeros (51.0,
##     2.50, 0.610), never whole for happening to be whole: lmn_fmt_sig();
##   a count (rows, draws, folds, clusters) whole: lmn_fmt_count();
##   a data value (an outcome, a time, a horizon) whole where it is whole,
##     three figures otherwise: lmn_fmt_data();
##   a share (a coverage, an observed or predicted rate) as a percentage to
##     three significant figures, as illume's prose writes them (81.0%, 4.64%,
##     100%): lmn_fmt_pct(). A share's standard error is an estimate in
##     percentage points ("SE 1.77 percentage points"). An interval's level
##     (95%) is a setting the user chose, not an estimate, and is written as
##     given;
##   a p-value by the family's own rule: lmn_fmt_p().
## Quantities from 1,000 up are grouped ("1,230"); years, identifiers and codes
## are data and never come through here.

## A p-value by the family's rule, to three significant figures, and
## "< 0.001" below that.
#' @keywords internal
#' @noRd
lmn_fmt_p <- function(p, prefix = "p") {
  if (is.na(p)) return(paste(prefix, "not available"))
  t <- ilm_disp(p, "p")$text
  if (startsWith(t, "<")) paste(prefix, t) else paste(prefix, "=", t)
}

## An ESTIMATE: three significant figures, trailing zeros kept, a thousands
## separator from 1,000 up.
#' @keywords internal
#' @noRd
lmn_fmt_sig <- function(v) ilm_disp(v, "signif", 3L)$text

## A COUNT: whole, with a thousands separator.
#' @keywords internal
#' @noRd
lmn_fmt_count <- function(v) ilm_disp(v, "fixed", 0L)$text

## A DATA VALUE: whole where it is a whole number, three figures otherwise.
#' @keywords internal
#' @noRd
lmn_fmt_data <- function(v)
  ifelse(is.finite(v) & v == round(v), ilm_disp(v, "fixed", 0L)$text, ilm_disp(v, "signif", 3L)$text)

## A SHARE as a percentage, to three significant figures, trailing zeros kept.
#' @keywords internal
#' @noRd
lmn_fmt_pct <- function(v) {
  t <- ilm_disp(100 * v, "signif", 3L)$text
  ifelse(is.na(t), NA_character_, paste0(t, "%"))
}

## ---- random numbers ---------------------------------------------------------

## Evaluate `expr` with the RNG seeded, then put the caller's RNG state back.
## A `seed` argument then makes a result reproducible without resetting the
## user's own random stream, which set.seed() inside a function would do.
#' @keywords internal
#' @noRd
lmn_with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  genv <- globalenv()
  had <- exists(".Random.seed", envir = genv, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = genv, inherits = FALSE)
  on.exit(if (had) assign(".Random.seed", old, envir = genv)
          else if (exists(".Random.seed", envir = genv, inherits = FALSE))
            rm(".Random.seed", envir = genv), add = TRUE)
  set.seed(seed)
  expr
}

## A prediction's seed as it was set: the value, all three of RNGkind() at
## the time, and the stream it drove, so the prediction can be rerun.
#' @keywords internal
#' @noRd
lmn_seed_record <- function(value, stream) {
  k <- RNGkind()
  list(value = value, kind = k[1], normal_kind = k[2], sample_kind = k[3], stream = stream)
}

## ---- argument checks --------------------------------------------------------

#' @keywords internal
#' @noRd
lmn_check_column <- function(data, name, arg) {
  if (is.null(name))
    stop("`", arg, "` is needed: the name of a column in `data`.", call. = FALSE)
  if (!is.character(name) || length(name) != 1L)
    stop("`", arg, "` must be one column name, as a string.", call. = FALSE)
  if (!name %in% names(data))
    stop("`", arg, "` names column '", name, "', which is not in `data`.",
         call. = FALSE)
  v <- data[[name]]
  if (anyNA(v))
    stop("column '", name, "' (`", arg, "`) has missing values; every row needs one.",
         call. = FALSE)
  v
}

#' @keywords internal
#' @noRd
lmn_check_count <- function(x, arg, min = 1L) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x != round(x) || x < min)
    stop("`", arg, "` must be a whole number of at least ", min, ".", call. = FALSE)
  as.integer(x)
}

## ---- quadrature ---------------------------------------------------------------

## Gauss-Legendre nodes and weights on [-1, 1] (Golub-Welsch). Used to
## integrate a predictive CDF for the CRPS of a continuous outcome.
#' @keywords internal
#' @noRd
lmn_gauss_legendre <- function(M) {
  i <- seq_len(M - 1L); b <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, M, M); J[cbind(i, i + 1L)] <- b; J[cbind(i + 1L, i)] <- b
  e <- eigen(J, symmetric = TRUE)
  o <- order(e$values)
  list(x = e$values[o], w = 2 * e$vectors[1, o]^2)
}
