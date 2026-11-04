# Predictive distributions from a fitted illume model

The engine every other prediction in `limn` goes through. For each row
of `newdata` it returns the distribution of the predicted quantity: the
mean, a new observation, or the linear predictor. That distribution
carries the fit's own uncertainty, drawn jointly from the Laplace
approximation that gives the fit its standard errors.

## Usage

``` r
ilm_predict(
  fit,
  newdata,
  target = c("mean", "response", "link"),
  groups,
  uncertainty = c("joint", "fixed", "latent", "none"),
  set = NULL,
  nsim = 1000L,
  seed = 1L,
  time = NULL,
  group = NULL
)
```

## Arguments

- fit:

  A fitted `"ilm_model"`.

- newdata:

  A data frame of rows to predict. Or a named list of data frames, one
  per scenario, stacked with a `.scenario` column. A `.draw` column
  gives covariates as draws; see *Covariates given as draws*.

- target:

  `"mean"`, the expected response; `"response"`, a new observation; or
  `"link"`, the linear predictor.

- groups:

  How the random effects enter; see *How the random effects enter*.
  Required when the model has random-effect or correlation-over-time
  terms.

- uncertainty:

  `"joint"`, `"fixed"`, `"latent"` or `"none"`; see *Uncertainty*.

- set:

  A named list of values to give every row, such as
  `list(dose = c(0, 10))`. Each row is repeated once per combination.

- nsim:

  Number of parameter draws.

- seed:

  Seed for the draws. The caller's random-number state is restored
  afterwards.

- time, group:

  The names of the columns of `newdata` holding time and the series,
  such as `time = "week", group = "region"`; see *Time and series
  columns*. Needed only where neither the model nor the rows say.

## Value

A predictive distribution, class `"ilm_pred"`: see *The prediction
object*.
[summary()](https://huttoncp.github.io/limn/reference/as.data.frame.ilm_pred.md)
and
[as.data.frame()](https://huttoncp.github.io/limn/reference/as.data.frame.ilm_pred.md)
give point estimates and intervals.

- For `"mean"` and `"link"`, the point estimate is the value at the
  estimates, the number
  [`predict()`](https://rdrr.io/r/stats/predict.html) gives.

- For `"response"`, it is the predictive median, with the mean beside
  it.

## How the random effects enter (`groups`)

A mixed model predicts different things for the same row, depending on
the question asked. So `groups` has no default for a model with random
effects or a correlation over time; it must be named. The choices:

- `"fitted"`: the unit the row names, with that unit's own estimated
  effects. A level the fit has not seen is predicted as a new unit, and
  marked in the `new_level` column.

- `"new"`: a new, unseen unit. Its effects are drawn from the fitted
  distribution of effects, afresh in every parameter draw. Rows naming
  the same unit share its effects within a draw.

- `"typical"`: a typical unit, with its effects at zero. For a nonlinear
  link this is not the average unit.

- `"population"`: the average over units. For `target = "mean"`, the
  mean is averaged over the distribution of effects, as illume's
  `predict(groups = "population")` does. For a new observation, a
  randomly chosen unit's outcome is the same thing as a new unit's, so
  it equals `"new"`.

Give one choice for every term, or a list per term, as in
`list(site = "fitted", ar = "new")`. The correlation over time is the
term `"ar"`. Smooth terms are part of the mean structure and always
enter as fitted.

## Over time

With a correlation over time (`ilm_ar1()`, `ilm_car1()` or `ilm_rw1()`),
a seen series under `"fitted"` is continued from its fitted cells:

- **after its last cell**, a forecast carried forward from that cell;

- **between two cells**, drawn from the bridge between them;

- **before its first cell**, carried backward from it.

A new series, or any series under `"new"`, is a fresh path. AR(1) and
CAR(1) paths start from their stationary distribution, and a random walk
starts from zero at the series' first time. Each series' new times are
drawn in sequence, each given the one before, so a draw is one path
through them, and rows naming the same series and time share one value,
across scenarios too. The model's own parameters are used throughout;
nothing is refitted.

For a mean or the linear predictor, the estimate at a seen series' new
time is the point forecast: the fitted cells carried forward with no new
shocks, at the estimates, as a Kalman filter's point forecast is. For a
new series, as for a new unit, it is the predictive median.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) gives the
expected value beside it, as `mean`; for counts it can be well above the
estimate at long horizons (see
[`as.data.frame.ilm_pred()`](https://huttoncp.github.io/limn/reference/as.data.frame.ilm_pred.md)).

## Time and series columns

A prediction records which of its columns hold time and the series, for
the trend check (see *Checks*) and for scaled errors read from the
history (`ilm_scores()`'s `history =`). It takes them from, in turn:

- a correlation over time built by name, such as
  `ilm_ar1(~ week | region)`: its own columns. Naming a different column
  in `time` or `group` is an error, since the model can use only its
  own, and so are rows from
  [`ilm_future()`](https://huttoncp.github.io/limn/reference/ilm_future.md)
  built around a different column that leave the model's own missing;

- `time` and `group`, which name columns of `newdata`;

- the rows
  [`ilm_future()`](https://huttoncp.github.io/limn/reference/ilm_future.md)
  built, which record their columns.

A model with time only in its formula, such as a trend or a random slope
in time, has no correlation over time to say. Rows from
[`ilm_future()`](https://huttoncp.github.io/limn/reference/ilm_future.md)
say for it; rows made otherwise, or changed by
[`merge()`](https://rdrr.io/r/base/merge.html) or
[`subset()`](https://rdrr.io/r/base/subset.html), which drop what
[`ilm_future()`](https://huttoncp.github.io/limn/reference/ilm_future.md)
recorded, need `time` and `group`. A correlation over time built from
vectors has no columns of its own: it takes the ones `time` and `group`
name, or else the ones
[`ilm_future()`](https://huttoncp.github.io/limn/reference/ilm_future.md)
recorded, and with neither, the prediction stops.

## Covariates given as draws

A covariate that is not known in advance, such as next month's
temperature, can be given as draws from its own forecast. `newdata` then
carries a `.draw` column: one block of rows per covariate draw, every
block listing the same rows in the same order, with only the drawn
covariates differing. Parameter draw s uses covariate draw ((s - 1) mod
C) + 1 of the C, so the prediction mixes over both. `nsim` is rounded up
to a multiple of C, with a message, so that each covariate draw is used
equally often.

The prediction has one row per row of a block. A drawn covariate shows
as `NA` in the rows where it varies, `meta` records C and the drawn
columns, and a drawn row's estimate is its predictive median, as a new
unit's is. Intervals conditional on a covariate's point forecast leave
out the uncertainty in forecasting it (fpp3, section 10.3). In limn's
own study (`studies/findings/covariate_draws.md`) they covered too
little: 0.75 to 0.83 for 90% gaussian intervals, against about 0.90 with
the covariate given as draws.
[`ilm_covariate_paths()`](https://huttoncp.github.io/limn/reference/ilm_covariate_paths.md)
builds the blocks from an illume model of the covariate, and
[`ilm_as_draws()`](https://huttoncp.github.io/limn/reference/ilm_as_draws.md)
from paths simulated by another tool.

## Uncertainty

- `"joint"` (the default): every parameter is drawn jointly, the
  variances included. In limn's prototype study of few groups
  (studies/prototypes/), it came closest to nominal coverage for a new
  unit's observation, though still below it: 0.88 for 90% intervals with
  8 groups, against 0.85-0.87 for the other methods it tried, each with
  a simulation standard error of about 0.01.

- `"fixed"`: the variance parameters are held at their estimates, and
  everything else is drawn.

- `"latent"`: every parameter is held at its estimate, the dispersion
  included, and only the random effects and the cells of a correlation
  over time are drawn, from their distribution given the parameters. For
  a gaussian model with a correlation over time this is a Kalman
  filter's forecast, which treats the parameters as known; limn's tests
  check it against
  [`stats::KalmanForecast()`](https://rdrr.io/r/stats/KalmanLike.html).
  ARIMA prediction intervals likewise leave out the variation in the
  parameter estimates (fpp3, section 9.8).

- `"none"`: the plug-in prediction at the estimates, with no parameter
  uncertainty. A new unit's effects, and a series' cells at times it was
  not fitted at (a forecast, a bridge or a backcast), are still drawn,
  and covariates given as draws are still mixed over; otherwise a new
  observation's distribution is the family's alone.

## Checks

What each check means and what to do about it, in the words it is
printed with. A reference supports what was found and why; a remedy
draws on it only where it says so, and anything marked as limn's own
reasoning is not the source's:

- **few groups, means** (warn): a random term has k levels, fewer than
  about 20. With so few groups the fit's normal approximation
  understates the uncertainty in the fixed effects, since it takes no
  account of the uncertainty in the variances (the reference shows this
  for linear mixed models), so intervals for means cover as z-intervals
  do: in limn's coverage study with 8 groups, a 90% interval for a
  typical unit's mean covered 0.84 to 0.86 for gaussian, Poisson,
  negative binomial and beta outcomes (binary outcomes covered 0.92).
  Reference: Kenward and Roger (1997), Biometrics 53, 983-997, pp.
  983-984. No remedy for these intervals has been measured yet, so read
  them as somewhat too narrow; more groups is the sure cure.

- **few groups, new units** (warn): a random term has k levels, fewer
  than about 20. With so few groups the variance between them is poorly
  determined: the reference found its standard error too small with 50
  groups or fewer, for linear models fitted by restricted maximum
  likelihood. In limn's coverage study with 8 groups, intervals for a
  new unit's observation covered 0.87 to 0.90 at 90%. Reference: Maas
  and Hox (2005), Methodology 1, 86-92, pp. 86 and 89-90. Refit with
  reml = TRUE in illume::ilm_model(): restricted maximum likelihood is
  almost always at least as good as maximum likelihood, and sometimes
  better, especially at estimating variances (the reference, p. 89,
  citing Browne 1998). In limn's prototype study of few groups
  (studies/prototypes/) it raised 8-group coverage from 0.872 to 0.892
  (Poisson) and from 0.891 to 0.903 (gaussian).

- **few groups, new units, REML fit** (warn): a random term has k
  levels, fewer than about 20. With so few groups the variance between
  them is poorly determined: the reference found its standard error too
  small with 50 groups or fewer, for linear models fitted by restricted
  maximum likelihood. In limn's coverage study with 8 groups, intervals
  for a new unit's observation covered 0.87 to 0.90 at 90%. Reference:
  Maas and Hox (2005), Methodology 1, 86-92, pp. 86 and 89-90. The fit
  already used restricted maximum likelihood (REML), the refit this
  check otherwise suggests. In limn's prototype study of few groups
  (studies/prototypes/), with REML and 8 groups of 6, 90% intervals for
  a new unit's observation covered 0.892 (Poisson) and 0.903 (gaussian),
  with simulation standard errors of 0.004 to 0.005. More groups is the
  sure cure.

- **unseen levels** (note): Some of the rows name a level the fit has
  not seen. Under groups = "fitted" a unit with no estimated effects can
  only be predicted as a new one. They are predicted as new units, with
  effects drawn from the fitted distribution, and marked in the
  new_level column. If they were meant to be units in the fit, check how
  their levels are spelt.

- **deterministic trend** (warn): The model trends in `time` (its trend
  terms). Forecasts past the data assume that the trend's slope does not
  change (the reference), so their intervals leave out any change in it
  (limn's own reasoning). A random slope in time gives each unit a trend
  of its own, and forecasts assume that each unit's slope does not
  change either (limn's own reasoning). Reference: Hyndman and
  Athanasopoulos (2021), Forecasting: Principles and Practice, 3rd
  edition, section 10.4. For longer horizons a stochastic trend is safer
  (the reference): a random walk, ar = ilm_rw1(~ time \| series), lets
  the level move, and the intervals widen with the horizon.

- **boundary draws** (fail): Some of the draws of a variance are not
  finite. Its estimate (SD near zero) is at the edge of its range, where
  the data cannot resolve it, and the fit did not hold it there. Predict
  with uncertainty = "fixed", which holds the variances at their
  estimates, or refit with boundary = "avoid" in illume::ilm_model(),
  which keeps the variance off the edge with a small penalty.

- **draws not computed** (fail): Some of the draws of the mean could not
  be computed. Averaging over units failed where a draw's variance
  reaches a large value, as happens when the variance is poorly
  determined, for instance with few groups. Predict with uncertainty =
  "fixed", which holds the variances at their estimates.

- **draws not computed** (fail): Some of the draws of the mean could not
  be computed. A draw's value overflowed or was undefined. Predict with
  uncertainty = "fixed", which holds the variances at their estimates.

- **draws not computed** (fail): Some of the draws of a new
  observation's distribution could not be computed. In those draws the
  linear predictor or the dispersion ran so far out that the family
  overflowed. That happens when a parameter sits at the edge of its
  range (a variance or residual SD near zero, a correlation near 1) and
  the fit did not hold it there, so its draws spread without limit.
  Predict with uncertainty = "fixed", which holds the variances at their
  estimates, or uncertainty = "latent", which holds the dispersion too.

## What is not here yet

- **Multinomial outcomes** come later, with forecasts of their latent
  state.

- **Binomial outcomes with more than one trial** are predicted per
  trial.

## The prediction object

Read a prediction through
[as.data.frame()](https://huttoncp.github.io/limn/reference/as.data.frame.ilm_pred.md)
or
[summary()](https://huttoncp.github.io/limn/reference/as.data.frame.ilm_pred.md)
where you can. Its elements are listed here so that nothing in it is
hidden. Every prediction, class `"ilm_pred"`, has:

- `kind`: `"response"` for a new observation's distribution, held as the
  family's parameters in each draw; `"mean"` or `"link"` for draws of
  the mean or the linear predictor; or `"sample"` for simulated
  outcomes, which is how
  [`ilm_aggregate()`](https://huttoncp.github.io/limn/reference/ilm_aggregate.md)
  holds some totals.

- `rows`: a data frame with one row per predicted row: `newdata`'s
  columns (one block's, with drawn covariates `NA`; see *Covariates
  given as draws*), and `new_level`, `TRUE` where a row names a level
  the fit has not seen. A total or an average has its `by` columns and
  `n_rows` instead, and an average also `weight`.

- `weights`: the draws' weights, which sum to 1.

- `meta`: a list recording how the prediction was made. It has `target`,
  `groups` (one `"term = choice"` per term), `uncertainty`, `draws` (the
  number of parameter draws), `nsim`, `family`, `seed`, and `response`,
  the name of the response column. With covariates given as draws it
  also has `covariate_draws` and `drawn`. Where its time and series
  columns are known (see *Time and series columns*), it has `time` and
  `series`, their names.

- `checks`: a data frame of the checks the prediction raised (see
  *Checks*), one row each, with columns `id`, `check`, `status`,
  `detail`, `cause`, `suggestion`, `reference` and `values`. It has no
  rows when nothing was found.

- `origin`: each row's forecast origin, its series' last observed time,
  or `NA` for a row that is not ahead of its series' data.

A prediction straight from `ilm_predict()` or
[`ilm_forecast()`](https://huttoncp.github.io/limn/reference/ilm_forecast.md)
also has `new_unit`, `TRUE` where a row's effects are drawn as a new
unit's or a new series'. A total or an average has `what` instead,
saying how it was taken, and `meta` records the same in `aggregate`.

What else it holds depends on `kind`:

- `"response"`: `dist`, the family's distribution functions as
  [`illume::ilm_dist()`](https://rdrr.io/pkg/illume/man/ilm_dist.html)
  gives them; `lp`, a named list of rows-by-draws matrices of linear
  predictors; `par`, the family's other parameters in each draw; and
  `pvec`, any parameter that is a vector in each draw, such as an
  ordered outcome's thresholds.

- `"mean"` and `"link"`: `draws`, a rows-by-draws matrix of the
  predicted values, and `estimate`, the value at the estimates (a new
  unit's or a drawn row's predictive median; see *Value*).

- `"sample"`: `draws`, a rows-by-samples matrix of simulated outcomes;
  `estimate`, their median; and `dist`, the family and whether it is
  discrete.

## References

Hyndman, R. J. and Athanasopoulos, G. (2021). *Forecasting: Principles
and Practice*, 3rd edition. OTexts. Section 9.8, "Forecasting" (ARIMA),
and section 10.3, "Forecasting" (dynamic regression) (online edition,
<https://otexts.com/fpp3/>, read 17 October 2026).

## See also

[`ilm_forecast()`](https://huttoncp.github.io/limn/reference/ilm_forecast.md)
for forecasts;
[`ilm_aggregate()`](https://huttoncp.github.io/limn/reference/ilm_aggregate.md)
and
[`ilm_average()`](https://huttoncp.github.io/limn/reference/ilm_average.md)
for totals and averages over rows;
[`ilm_pit()`](https://huttoncp.github.io/limn/reference/ilm_pit.md) and
[`ilm_scores.ilm_pred()`](https://huttoncp.github.io/limn/reference/ilm_scores.ilm_pred.md)
to check it against outcomes;
[`ilm_plot_fan()`](https://huttoncp.github.io/limn/reference/ilm_plot_fan.md).

## Examples

``` r
set.seed(1)
d <- data.frame(site = factor(rep(1:24, each = 5)), x = rnorm(120))
d$y <- rpois(120, exp(0.5 + 0.4 * d$x + rnorm(24, 0, 0.5)[d$site]))
fit <- illume::ilm_model(y ~ x + (1 | site), data = d, family = "poisson",
                         verbose = FALSE)
nd <- data.frame(site = c("1", "2", "new"), x = 0)
as.data.frame(ilm_predict(fit, nd, groups = "fitted", nsim = 200))
#>   site x new_level estimate     mean lower_0.9 upper_0.9
#> 1    1 0     FALSE 1.248385 1.309497 0.8112306  2.001837
#> 2    2 0     FALSE 2.091386 2.160922 1.5311525  3.048817
#> 3  new 0      TRUE 1.545643 1.676481 0.8835637  2.972637
as.data.frame(ilm_predict(fit, nd, target = "response", groups = "new",
                          nsim = 200))
#>   site x new_level estimate     mean lower_0.9 upper_0.9
#> 1    1 0     FALSE        1 1.672268         0         4
#> 2    2 0     FALSE        1 1.598920         0         4
#> 3  new 0     FALSE        1 1.659086         0         4
```
