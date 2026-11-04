# Scores by forecast horizon

Plots each model's mean score at each forecast horizon, with confidence
intervals, as in the reference's plot of forecast accuracy against the
horizon (section 5.10, Figure 5.25). Forecasts usually get worse with
the horizon. The question is how fast, and whether one model's advantage
holds further out.

## Usage

``` r
ilm_plot_horizon(..., rule = NULL, level = 0.95)
```

## Arguments

- ...:

  One or more tables of scores with a `horizon` column, such as
  rolling-origin cross-validation scores. Name them to label the models.

- rule:

  The score to plot; by default the first score in the first table, as
  in
  [`ilm_compare()`](https://huttoncp.github.io/limn/reference/ilm_compare.md).

- level:

  Confidence level of the intervals.

## Value

The data behind the plot, invisibly: one row per model and horizon, with
the mean score, its standard error and interval, and the numbers of rows
and origins.

## Reading the plot

The horizontal axis is the forecast horizon; the vertical axis is the
mean score there, and lower is better. Each model's points are joined by
a line.

- A line that rises is a model whose forecasts get worse further ahead,
  as most do. A steeper rise than another model's means its advantage,
  if any, fades with the horizon.

- Where two models' intervals at a horizon do not overlap, their scores
  there differ. Where they overlap, the plot cannot say:
  [`ilm_compare()`](https://huttoncp.github.io/limn/reference/ilm_compare.md)
  tests the difference directly, with a standard error of the difference
  itself, which allows for the two models' scores on the same rows being
  correlated.

## Standard errors

Rows forecast from the same origin share that origin's surprises, so
when the scores have an `origin` column, each horizon's standard error
is clustered by origin. This is the same cluster-robust estimator
[`ilm_compare()`](https://huttoncp.github.io/limn/reference/ilm_compare.md)
uses. Otherwise the rows are treated as independent.

## References

Hyndman, R. J. and Athanasopoulos, G. (2021). *Forecasting: Principles
and Practice*, 3rd edition. OTexts. Section 5.10, "Time series
cross-validation", Figure 5.25 (online edition,
<https://otexts.com/fpp3/toolbox.html#sec-tscv>, read 27 October 2026).

## See also

[`ilm_compare()`](https://huttoncp.github.io/limn/reference/ilm_compare.md)
to test whether two models' scores differ;
[`ilm_folds()`](https://huttoncp.github.io/limn/reference/ilm_folds.md)
for rolling-origin folds;
[`ilm_cv()`](https://huttoncp.github.io/limn/reference/ilm_cv.md) for
the scores.

## Examples

``` r
set.seed(1)
s <- expand.grid(horizon = 1:4, origin = 1:12, site = 1:5)
ar <- data.frame(s, log = 1 + 0.08 * s$horizon + rnorm(nrow(s), 0, 0.3))
iid <- data.frame(s, log = ar$log + 0.03 * s$horizon + rnorm(nrow(s), 0, 0.1))
ilm_plot_horizon(ar = ar, iid = iid)
```
