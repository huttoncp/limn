# Cross-validation folds that follow the data's structure

Splits the rows of `data` into training and test sets for
cross-validation. The folds are row indices, so they can be inspected,
and reused for every model being compared: two models scored on the same
folds can be compared fold by fold.

## Usage

``` r
ilm_folds(
  data,
  type = c("random", "grouped", "blocked", "rolling"),
  K = 5L,
  group = NULL,
  time = NULL,
  initial = NULL,
  horizon = 1L,
  step = 1L,
  window = c("expanding", "sliding"),
  size = NULL,
  gap = 0L,
  known = NULL,
  seed = 1L
)
```

## Arguments

- data:

  A data frame: the rows the models are fitted to.

- type:

  The kind of fold; see *Which kind of fold*.

- K:

  Number of folds, for random, grouped and blocked folds.

- group:

  For grouped folds, the name of the grouping column in `data`.

- time:

  For blocked and rolling folds, the name of the time column in `data`:
  numeric, a `Date` or a date-time.

- initial:

  For rolling folds, the number of time points in the first training
  window.

- horizon:

  For rolling folds, the number of time points in each test set.

- step:

  For rolling folds, how many time points the origin moves from one fold
  to the next.

- window:

  For rolling folds: `"expanding"` trains on every time point up to the
  origin, `"sliding"` on the last `size` of them.

- size:

  For a sliding window, the number of time points it trains on. Defaults
  to `initial`.

- gap:

  Time points left out between training and test rows: after the origin
  for rolling folds, and on both sides of the test block for blocked
  folds.

- known:

  For rolling folds, the names of the covariates whose future values
  would be known at the origin; see *Ex-ante and ex-post*.

- seed:

  Seed for random and grouped folds. The caller's random-number state is
  restored afterwards.

## Value

An object of class `"ilm_folds"`: a list with one element per fold, each
holding `train` and `test` row indices, and for rolling folds the fold's
`origin` and each test row's `horizon`.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) gives one
row per data row and fold.

## Which kind of fold

The folds decide which question the cross-validation answers, so they
should match the prediction being made.

- `"random"`: rows are held out at random. Every held-out row's group is
  still in the training data, so this scores predictions for **units the
  model has seen**. On data with a time structure it scores
  interpolation, not forecasting.

- `"grouped"`: whole groups are held out together, so this scores
  predictions for **new units**. Random folds would overstate how well a
  model predicts a unit it has never seen.

- `"blocked"`: contiguous blocks of time are held out, optionally with a
  gap either side so that rows next to the block do not leak into
  training.

- `"rolling"`: rolling-origin evaluation, also called time-series
  cross-validation. Each fold trains on the times up to an origin and
  tests on the `horizon` times after it, so every test row is a genuine
  forecast. For a model with an autoregressive term, this is the scheme
  that scores forecasts. Random and blocked folds let the model see both
  sides of a held-out time, so they score interpolation and overstate
  forecast skill.

Rolling origins are the reference's "evaluation on a rolling forecasting
origin": every training set holds only rows before its test rows, and
the scores are averaged over the test sets. Which question each of the
other kinds answers is limn's own reasoning.

## Time points and horizons

Blocked and rolling folds count **time points**: the distinct values of
`time` in `data`, in order. Every row at a time point goes into the same
set, so a panel of several series is split at the same times for all of
them. A test row's horizon is how many time points after the origin it
lies, so with `gap = 0` the first test time point has horizon 1.

## Ex-ante and ex-post

A forecast is *ex-ante* when it uses only what would have been known at
the origin. Covariates the model uses at their realised future values,
such as the weather on the day being forecast, make it *ex-post*, and
flatter the forecast. Name in `known` the covariates whose future values
would genuinely be known, such as a calendar or a planned dose. The
folds record it so the evaluation can say which kind it is.

## References

Hyndman, R. J. and Athanasopoulos, G. (2021). *Forecasting: Principles
and Practice*, 3rd edition. OTexts. Section 5.10, "Time series
cross-validation" (online edition,
<https://otexts.com/fpp3/toolbox.html#sec-tscv>, read 27 October 2026).

## See also

[`ilm_cv()`](https://huttoncp.github.io/limn/reference/ilm_cv.md), which
scores a model on the folds;
[`ilm_compare()`](https://huttoncp.github.io/limn/reference/ilm_compare.md),
which compares two models' scores on the same folds.

## Examples

``` r
d <- data.frame(site = rep(letters[1:6], each = 8), month = rep(1:8, 6),
                y = rnorm(48))
ilm_folds(d, "grouped", K = 3, group = "site")
#> <ilm_folds> 3 grouped folds (whole groups of 'site' held out)
#>   48 rows; training rows 32, test rows 16 per fold
f <- ilm_folds(d, "rolling", time = "month", initial = 5, horizon = 2)
f
#> <ilm_folds> 2 rolling origins over time ('month'), expanding window
#>   48 rows; training rows 30 to 36, test rows 12 per fold
#>   first origin after 5 time points, then every 1
#>   horizons 1 to 2
#>   known in advance: none named, so any covariate makes the evaluation ex-post
head(as.data.frame(f))
#>   fold row   set origin horizon
#> 1    1   1 train      5      NA
#> 2    1   2 train      5      NA
#> 3    1   3 train      5      NA
#> 4    1   4 train      5      NA
#> 5    1   5 train      5      NA
#> 6    1   9 train      5      NA
```
