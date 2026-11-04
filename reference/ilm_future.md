# The rows a forecast will need: the next time points of every series

Builds a data frame with the next `h` time points after each series'
last observed one, ready to be filled in and forecast. Columns that are
constant within each series, such as a site's region, are copied
forward. Columns that vary over time are not guessed: they are left
missing and named, so that the forecast uses values someone chose, not
ones invented here.

## Usage

``` r
ilm_future(data, time, group = NULL, h = 1L, step = NULL, response = NULL)
```

## Arguments

- data:

  A data frame: the observed history.

- time:

  The name of the time column: numeric, a `Date` or a date-time.

- group:

  The name of the column identifying each series, or `NULL` for a single
  series.

- h:

  How many time points ahead.

- step:

  The step between time points: a number for numeric time or date-times
  (in seconds), and for dates a number of days or one of `"day"`,
  `"week"`, `"month"`, `"quarter"` or `"year"`. By default it is taken
  from the data.

- response:

  The name of the response column, if `data` has one. It is left missing
  without being listed, since it is what the forecast predicts.

## Value

A data frame with one row per series and future time point, in the
columns of `data`, plus `horizon` (1 to `h`). The names of the columns
left missing are in `attr(, "unknown")`, and the time and series columns
in `attr(, "time")` and `attr(, "group")`, which
[`ilm_forecast()`](https://huttoncp.github.io/limn/reference/ilm_forecast.md)'s
checks read.

## Covariates that vary over time

A forecast needs the future values of every covariate the model uses.
For some they are known in advance, such as a calendar, a planned dose
or a scheduled price. For others they are not, and there are three
honest choices:

- plug in a scenario, and say that is what the forecast assumes;

- forecast the covariate itself first;

- supply several plausible paths, so that the forecast's uncertainty
  includes the covariate's.
  [`ilm_covariate_paths()`](https://huttoncp.github.io/limn/reference/ilm_covariate_paths.md)
  draws them from an illume model of the covariate, and
  [`ilm_as_draws()`](https://huttoncp.github.io/limn/reference/ilm_as_draws.md)
  takes them from another tool.

`ilm_future()` makes none of these choices for you. It leaves those
columns missing and names them in a message. A forecast made with
realised future values of a covariate is *ex-post*: in fpp3's terms
(section 7.6), made with later information on the predictors, not a
genuine forecast. It flatters the model, which is limn's own reading;
see
[`ilm_folds()`](https://huttoncp.github.io/limn/reference/ilm_folds.md).

## Time steps

The step between time points is taken from the data unless you give
`step`.

- For numeric time, it is the one positive difference the series share.

- For dates, it is a calendar unit when every step is one: a day, a
  week, a month, a quarter or a year. So monthly dates advance to the
  same day of the next month, not by a fixed 30 days.

Irregular spacing has no single next time, so then `step` must be given.

## References

Hyndman, R. J. and Athanasopoulos, G. (2021). *Forecasting: Principles
and Practice*, 3rd edition. OTexts. Section 7.6, "Forecasting with
regression", ex-ante and ex-post forecasts (online edition,
<https://otexts.com/fpp3/regression.html#sec-forecasting-regression>,
read 27 October 2026).

## See also

[`ilm_forecast()`](https://huttoncp.github.io/limn/reference/ilm_forecast.md),
which forecasts the rows;
[`ilm_covariate_paths()`](https://huttoncp.github.io/limn/reference/ilm_covariate_paths.md)
and
[`ilm_as_draws()`](https://huttoncp.github.io/limn/reference/ilm_as_draws.md)
for covariates not known in advance;
[`ilm_folds()`](https://huttoncp.github.io/limn/reference/ilm_folds.md)
for rolling-origin evaluation of forecasts.

## Examples

``` r
d <- data.frame(site = rep(c("a", "b"), each = 4), month = rep(1:4, 2),
                region = rep(c("north", "south"), each = 4),
                rain = rnorm(8), y = rpois(8, 3))
ilm_future(d, time = "month", group = "site", h = 2, response = "y")
#> These columns vary over time, so their future values are not known: rain. They are left missing: fill them in with the values the forecast should assume.
#>   site month region rain  y horizon
#> 1    a     5  north   NA NA       1
#> 2    a     6  north   NA NA       2
#> 3    b     5  south   NA NA       1
#> 4    b     6  south   NA NA       2
```
