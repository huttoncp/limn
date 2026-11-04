# Totals and averages of predictions over rows

Sums or averages a prediction over groups of its rows: a series'
forecasts over a horizon, every series at each time, or all rows. The
sum is taken draw by draw. So whatever the rows share within a draw (the
parameters, a unit's effects, a series' path over time) carries into the
total, and totals at different levels add up.

## Usage

``` r
ilm_aggregate(pred, by = NULL, fun = c("sum", "mean"))
```

## Arguments

- pred:

  A prediction from
  [`ilm_predict()`](https://huttoncp.github.io/limn/reference/ilm_predict.md)
  or
  [`ilm_forecast()`](https://huttoncp.github.io/limn/reference/ilm_forecast.md).

- by:

  Names of columns in the prediction's rows, as
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) shows
  them. There is one total for each combination of their values. `NULL`,
  the default, totals every row.

- fun:

  `"sum"`, the default, or `"mean"`.

## Value

A prediction for the totals (see *The prediction object* in
[`ilm_predict()`](https://huttoncp.github.io/limn/reference/ilm_predict.md)),
with one row per combination of `by` (and of scenario), and `n_rows`,
the number of rows each one combines. Its `what` and `meta$aggregate`
say how the total was taken.

## How a total is computed

- **A prediction of the mean:** each draw's means are summed, and the
  estimate is the sum of the rows' estimates.

- **New observations from a Poisson or gaussian model with no zero
  part:** within a draw, a sum of such outcomes is again Poisson (with
  the summed mean) or gaussian (with the summed mean and variance). So
  the total is an exact mixture, and it is scored like any other
  prediction.

- **New observations from any other family:** outcomes are simulated
  from each draw, four per draw, and summed. The total is then held as
  those samples:

  - its CDF, quantiles, mean and PIT are the samples';

  - its CRPS is the fair estimator from them (Ferro 2014);

  - it has no log score, since samples have no density.

  A fixed internal seed means totals over different groupings of the
  same prediction add up.

- **The linear predictor** has no total.

The method used is printed, and kept in `meta$aggregate`. An average
(`fun = "mean"`) is the total over the number of rows. Where the family
does not stay closed under that division, it is scored through the
total, so an average of counts keeps its whole-number steps. Rows of
different scenarios are never added together.

## Checks

What each check means and what to do about it, in the words it is
printed with. A reference supports what was found and why; a remedy
draws on it only where it says so, and anything marked as limn's own
reasoning is not the source's:

- **totals not simulated** (fail): Some of the simulated totals are
  missing. The family could not draw outcomes from some parameter draws,
  as happens when a dispersion parameter sits at its limit and the fit
  did not hold it there (a negative binomial whose size runs to
  infinity, say). If the data show no overdispersion beyond the model's
  random terms, refit with family = "poisson" in illume::ilm_model(),
  whose totals are exact.

## References

Ferro, C. A. T. (2014). Fair scores for ensemble forecasts. *Quarterly
Journal of the Royal Meteorological Society* 140, 1917-1923.

## See also

[`ilm_forecast()`](https://huttoncp.github.io/limn/reference/ilm_forecast.md),
[`ilm_future()`](https://huttoncp.github.io/limn/reference/ilm_future.md);
[`ilm_average()`](https://huttoncp.github.io/limn/reference/ilm_average.md)
for weighted averages of means.

## Examples

``` r
set.seed(1)
d <- data.frame(clinic = factor(rep(sprintf("c%02d", 1:20), each = 24)),
                month = rep(1:24, 20))
d$visits <- rpois(nrow(d), exp(1.5 + ave(rnorm(nrow(d), 0, 0.15), d$clinic,
                                         FUN = cumsum)))
fit <- illume::ilm_model(visits ~ 1, data = d, family = "poisson",
                         ar = illume::ilm_rw1(~ month | clinic), verbose = FALSE)
fut <- ilm_future(d, time = "month", group = "clinic", h = 3, response = "visits")
fc <- ilm_forecast(fit, fut, nsim = 500)
## each clinic's total over the quarter, and every clinic's total by month
head(as.data.frame(ilm_aggregate(fc, by = "clinic")))
#>   clinic n_rows estimate      mean lower_0.9 upper_0.9
#> 1    c01      3       20 21.095366        10        36
#> 2    c02      3       17 17.597244         8        31
#> 3    c03      3       38 39.768869        22        64
#> 4    c04      3       17 18.215321         8        31
#> 5    c05      3        9  9.940684         4        18
#> 6    c06      3        7  7.420347         2        15
as.data.frame(ilm_aggregate(fc, by = "month"))
#>   month n_rows estimate     mean lower_0.9 upper_0.9
#> 1    25     20      123 123.0562       101       146
#> 2    26     20      124 124.4171       101       150
#> 3    27     20      126 126.4700       102       154
```
