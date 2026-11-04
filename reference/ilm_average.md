# Weighted averages of predicted means: standardisation and post-stratification

Averages a prediction of the mean over groups of its rows, with weights,
draw by draw. With the rows a reference population and equal weights, it
is standardisation: the population's average mean under the model, for
example with a covariate set to each of several values. With the rows
population cells and their sizes as weights, it is post-stratification
(as in multilevel regression and post-stratification, MRP).

## Usage

``` r
ilm_average(pred, by = NULL, weights = NULL)
```

## Arguments

- pred:

  A prediction of the mean, from
  [`ilm_predict()`](https://huttoncp.github.io/limn/reference/ilm_predict.md)
  or
  [`ilm_forecast()`](https://huttoncp.github.io/limn/reference/ilm_forecast.md)
  with `target = "mean"`.

- by:

  Names of columns in the prediction's rows, as
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) shows
  them. There is one average for each combination of their values.
  `NULL`, the default, averages every row. Rows of different scenarios
  are never averaged together.

- weights:

  The rows' weights: a non-negative number per row, or the name of a
  column in the prediction's rows. They are normalised within each
  group, so only their ratios matter. `NULL`, the default, weights the
  rows equally.

## Value

A prediction of the mean (see *The prediction object* in
[`ilm_predict()`](https://huttoncp.github.io/limn/reference/ilm_predict.md)),
with one row per combination of `by` (and of scenario), `n_rows`, the
number of rows each average combines, and `weight`, their total weight.

## Details

The average is taken within each draw, so its uncertainty is the joint
uncertainty of the rows' means: whatever they share within a draw (the
parameters, a unit's effect) carries into the average.

Post-stratification weights each cell's predicted mean by the cell's
population count, as Gelman and Hill (2007) do for state opinion from
national polls (section 14.1, equation (14.1)); their code in the same
section computes the weighted average in each simulation draw, as here.

## Means, not new observations

Standardisation and post-stratification average what the model expects,
so `pred` must be a prediction of the mean (`target = "mean"`). A
weighted average of new observations is a different question, whose
answer depends on the sizes of the groups being averaged. For an
unweighted average of new observations, use
[`ilm_aggregate()`](https://huttoncp.github.io/limn/reference/ilm_aggregate.md)
with `fun = "mean"`.

## References

Gelman, A. and Hill, J. (2007). *Data Analysis Using Regression and
Multilevel/Hierarchical Models*. Cambridge University Press. Section
14.1, "State-level opinions from national polls": equation (14.1), and
the weighted average computed in each simulation draw.

## See also

[`ilm_predict()`](https://huttoncp.github.io/limn/reference/ilm_predict.md)
with `target = "mean"` for the rows' means;
[`ilm_aggregate()`](https://huttoncp.github.io/limn/reference/ilm_aggregate.md)
for totals, and averages of new observations.

## Examples

``` r
set.seed(1)
d <- data.frame(clinic = factor(rep(sprintf("c%02d", 1:20), each = 10)),
                age = rep(c("young", "old"), 100), x = rnorm(200))
d$visits <- rpois(200, exp(0.5 + 0.4 * d$x + 0.3 * (d$age == "old") +
                           rnorm(20, 0, 0.4)[d$clinic]))
fit <- illume::ilm_model(visits ~ x + age + (1 | clinic), data = d,
                         family = "poisson", verbose = FALSE)
## post-stratification: the mean for a population that is 70% old (the
## population average does not use the clinic, but the rows need its column)
cells <- data.frame(age = c("young", "old"), x = 0, size = c(300, 700), clinic = NA)
p <- ilm_predict(fit, cells, target = "mean", groups = "population", nsim = 500)
as.data.frame(ilm_average(p, weights = "size"))
#>   n_rows weight estimate     mean lower_0.9 upper_0.9
#> 1      2   1000 2.604135 2.645949  2.211712  3.192428
```
