#' The rows a forecast will need: the next time points of every series
#'
#' Builds a data frame with the next `h` time points after each series' last
#' observed one, ready to be filled in and forecast. Columns that are constant
#' within each series, such as a site's region, are copied forward. Columns
#' that vary over time are not guessed: they are left missing and named, so
#' that the forecast uses values someone chose, not ones invented here.
#'
#' @section Covariates that vary over time:
#' A forecast needs the future values of every covariate the model uses. For
#' some they are known in advance, such as a calendar, a planned dose or a
#' scheduled price. For others they are not, and there are three honest
#' choices:
#'
#' * plug in a scenario, and say that is what the forecast assumes;
#' * forecast the covariate itself first;
#' * supply several plausible paths, so that the forecast's uncertainty
#'   includes the covariate's. [ilm_covariate_paths()] draws them from an
#'   illume model of the covariate, and [ilm_as_draws()] takes them from
#'   another tool.
#'
#' `ilm_future()` makes none of these choices for you. It leaves those columns
#' missing and names them in a message. A forecast made with realised future
#' values of a covariate is *ex-post*: in fpp3's terms (section 7.6), made
#' with later information on the predictors, not a genuine forecast. It
#' flatters the model, which is limn's own reading; see [ilm_folds()].
#'
#' @section Time steps:
#' The step between time points is taken from the data unless you give
#' `step`.
#' * For numeric time, it is the one positive difference the series share.
#' * For dates, it is a calendar unit when every step is one: a day, a week,
#'   a month, a quarter or a year. So monthly dates advance to the same day of
#'   the next month, not by a fixed 30 days.
#'
#' Irregular spacing has no single next time, so then `step` must be given.
#'
#' @param data A data frame: the observed history.
#' @param time The name of the time column: numeric, a `Date` or a date-time.
#' @param group The name of the column identifying each series, or `NULL` for
#'   a single series.
#' @param h How many time points ahead.
#' @param step The step between time points: a number for numeric time or
#'   date-times (in seconds), and for dates a number of days or one of `"day"`,
#'   `"week"`, `"month"`, `"quarter"` or `"year"`. By default it is taken from
#'   the data.
#' @param response The name of the response column, if `data` has one. It is
#'   left missing without being listed, since it is what the forecast predicts.
#' @return A data frame with one row per series and future time point, in the
#'   columns of `data`, plus `horizon` (1 to `h`). The names of the columns
#'   left missing are in `attr(, "unknown")`, and the time and series columns
#'   in `attr(, "time")` and `attr(, "group")`, which [ilm_forecast()]'s
#'   checks read.
#' @references Hyndman, R. J. and Athanasopoulos, G. (2021). *Forecasting:
#'   Principles and Practice*, 3rd edition. OTexts. Section 7.6, "Forecasting
#'   with regression", ex-ante and ex-post forecasts (online edition,
#'   <https://otexts.com/fpp3/regression.html#sec-forecasting-regression>,
#'   read 27 October 2026).
#' @seealso [ilm_forecast()], which forecasts the rows; [ilm_covariate_paths()]
#'   and [ilm_as_draws()] for covariates not known in advance; [ilm_folds()]
#'   for rolling-origin evaluation of forecasts.
#' @examples
#' d <- data.frame(site = rep(c("a", "b"), each = 4), month = rep(1:4, 2),
#'                 region = rep(c("north", "south"), each = 4),
#'                 rain = rnorm(8), y = rpois(8, 3))
#' ilm_future(d, time = "month", group = "site", h = 2, response = "y")
#' @export
ilm_future <- function(data, time, group = NULL, h = 1L, step = NULL,
                       response = NULL) {
  if (!is.data.frame(data))
    stop("`data` must be a data frame; it is ", class(data)[1], ".", call. = FALSE)
  tm <- lmn_check_column(data, time, "time")
  h <- lmn_check_count(h, "h")
  if (!is.numeric(tm) && !inherits(tm, c("Date", "POSIXt")))
    stop("the time column '", time, "' must be numeric, a Date or a date-time.",
         call. = FALSE)
  g <- if (is.null(group)) rep("1", nrow(data)) else
    as.character(lmn_check_column(data, group, "group"))
  if ("horizon" %in% names(data))
    stop("`data` already has a column called 'horizon', which the result adds.",
         call. = FALSE)
  if (!is.null(response) && !response %in% names(data))
    stop("`response` names column '", response, "', which is not in `data`.",
         call. = FALSE)
  adv <- lmn_time_step(tm, g, step)

  keys <- unique(g)
  last <- vapply(keys, function(k) which(g == k)[which.max(tm[g == k])], 1L)
  out <- data[rep(last, each = h), , drop = FALSE]
  rownames(out) <- NULL
  out[[time]] <- do.call(c, lapply(seq_along(keys), function(i)
    adv(tm[last[i]], h)))
  out$horizon <- rep(seq_len(h), length(keys))

  ## columns that vary within a series are not carried forward
  other <- setdiff(names(data), c(time, group))
  varies <- vapply(other, function(v) {
    x <- data[[v]]
    any(tapply(seq_along(x), g, function(ix) length(unique(x[ix])) > 1L))
  }, TRUE)
  unknown <- other[varies]
  for (v in unknown) out[[v]] <- out[[v]][NA_integer_]
  listed <- setdiff(unknown, response)
  if (length(listed))
    message("These columns vary over time, so their future values are not known: ",
            paste(listed, collapse = ", "), ". They are left missing: fill them ",
            "in with the values the forecast should assume.")
  attr(out, "unknown") <- listed
  attr(out, "time") <- time
  attr(out, "group") <- group
  out
}

## A function giving the h time points after t, with the step taken from the
## data (the one difference every series shares, or a calendar unit for dates)
## or from `step`.
#' @keywords internal
#' @noRd
lmn_time_step <- function(tm, g, step) {
  units <- c(day = "day", week = "week", month = "month", quarter = "quarter",
             year = "year")
  if (inherits(tm, "Date")) {
    if (is.null(step)) step <- lmn_date_unit(tm, g)
    if (is.character(step)) {
      if (!step %in% units)
        stop("a date step must be a number of days or one of ",
             paste(units, collapse = ", "), ".", call. = FALSE)
      if (step %in% c("day", "week")) {
        k <- if (step == "day") 1 else 7
        return(function(t, h) t + k * seq_len(h))
      }
      k <- c(month = 1L, quarter = 3L, year = 12L)[[step]]
      return(function(t, h) lmn_add_months(t, k * seq_len(h)))
    }
    return(function(t, h) t + step * seq_len(h))
  }
  if (is.null(step)) {
    d <- unlist(tapply(as.numeric(tm), g, function(x) diff(sort(unique(x))),
                       simplify = FALSE))
    d <- unique(signif(d[d > 0], 10))
    if (length(d) != 1L)
      stop("the time steps are not all equal, so there is no single next time; ",
           "give `step`.", call. = FALSE)
    step <- d
  }
  if (!is.numeric(step) || length(step) != 1L || step <= 0)
    stop("`step` must be one positive number here.", call. = FALSE)
  function(t, h) t + step * seq_len(h)
}

## Calendar months after a date. seq(by = "month") rolls 31 March on to 1 May,
## so this works from the first of the month instead: a month-end stays a
## month-end, and any other day is capped at the length of its month.
#' @keywords internal
#' @noRd
lmn_add_months <- function(t, m) {
  lt <- as.POSIXlt(t)
  y <- lt$year + 1900L; mo <- lt$mon + m                       # months from January
  first <- as.Date(sprintf("%04d-%02d-01", y + mo %/% 12L, mo %% 12L + 1L))
  next_first <- as.Date(sprintf("%04d-%02d-01", y + (mo + 1L) %/% 12L,
                                (mo + 1L) %% 12L + 1L))
  month_len <- as.integer(next_first - first)
  own_first <- as.Date(sprintf("%04d-%02d-01", y, lt$mon + 1L))
  own_len <- as.integer(lmn_add_first(own_first) - own_first)
  day <- if (lt$mday == own_len) month_len else pmin(lt$mday, month_len)
  first + day - 1L
}

## The first day of the month after `first`, itself a first of the month.
#' @keywords internal
#' @noRd
lmn_add_first <- function(first) {
  lt <- as.POSIXlt(first)
  as.Date(sprintf("%04d-%02d-01", lt$year + 1900L + (lt$mon + 1L) %/% 12L,
                  (lt$mon + 1L) %% 12L + 1L))
}

## The calendar unit every step of a date series is one of, if there is one.
#' @keywords internal
#' @noRd
lmn_date_unit <- function(tm, g) {
  d <- unlist(tapply(as.numeric(tm), g, function(x) diff(sort(unique(x))),
                     simplify = FALSE))
  if (!length(d)) stop("each series has one date, so there is no step to take ",
                       "from the data; give `step`.", call. = FALSE)
  if (all(d == 1)) return("day")
  if (all(d == 7)) return("week")
  if (all(d >= 28 & d <= 31)) return("month")
  if (all(d >= 89 & d <= 92)) return("quarter")
  if (all(d >= 365 & d <= 366)) return("year")
  d <- unique(d)
  if (length(d) == 1L) return(d)
  stop("the dates are not evenly spaced in days, weeks, months, quarters or ",
       "years, so there is no single next date; give `step`.", call. = FALSE)
}
