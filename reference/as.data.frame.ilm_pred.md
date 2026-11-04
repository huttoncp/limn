# Point estimates and intervals from a predictive distribution

Point estimates and intervals from a predictive distribution

## Usage

``` r
# S3 method for class 'ilm_pred'
as.data.frame(x, row.names = NULL, optional = FALSE, ..., level = 0.9)

# S3 method for class 'ilm_pred'
summary(object, level = c(0.5, 0.9), ...)
```

## Arguments

- x, object:

  A predictive distribution from
  [`ilm_predict()`](https://huttoncp.github.io/limn/reference/ilm_predict.md).

- row.names, optional, ...:

  Unused; for the generics.

- level:

  Coverage of the central intervals.

## Value

A data frame: the rows predicted, the point estimate, the mean, and for
each level the interval's lower and upper limits.

- For a new observation, the estimate is the predictive median.

- For a mean or a linear predictor, the estimate is the value at the
  estimates: for a seen series' new time, its fitted cells carried
  forward with no new shocks. For a gaussian model that is a Kalman
  filter's point forecast; limn's tests check it against
  [`stats::KalmanForecast()`](https://rdrr.io/r/stats/KalmanLike.html).

- `mean` is the expected value, averaged over everything the prediction
  is uncertain about: the parameters, a new unit's effects, and a
  series' future shocks. On a log or logit scale it can sit well above
  the estimate, most for forecasts at long horizons and for new units.
  Use it for an expected count, for totals (expected values add up,
  medians do not) and for squared error; use the estimate for a central
  value. For a forecast, `mean` is the point forecast as fpp3 and its
  fable package give it: the mean of the forecast distribution, in
  fable's `.mean` column (sections 5.1 and 5.5).

- An ordered outcome is summarised in its categories, with no mean.

[`summary()`](https://rdrr.io/r/base/summary.html) is
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) with the
50% and 90% intervals.

## References

Hyndman, R. J. and Athanasopoulos, G. (2021). *Forecasting: Principles
and Practice*, 3rd edition. OTexts. Section 5.1, the point forecast in
fable's `.mean` column, and section 5.5, "Distributional forecasts and
prediction intervals" (online edition, <https://otexts.com/fpp3/>, read
17 October 2026).

## See also

[`ilm_predict()`](https://huttoncp.github.io/limn/reference/ilm_predict.md),
which makes the prediction;
[`ilm_aggregate()`](https://huttoncp.github.io/limn/reference/ilm_aggregate.md)
for totals.

## Examples

``` r
set.seed(1)
d <- data.frame(site = factor(rep(1:12, each = 6)), x = rnorm(72))
d$y <- rpois(72, exp(0.5 + 0.3 * d$x + rnorm(12, 0, 0.4)[d$site]))
fit <- illume::ilm_model(y ~ x + (1 | site), data = d, family = "poisson",
                         verbose = FALSE)
p <- ilm_predict(fit, data.frame(site = c("1", "new"), x = 0),
                 target = "response", groups = "fitted", nsim = 200)
#> Warning: 'site' has 12 levels, fewer than about 20. With so few groups the variance between them is poorly determined: the reference found its standard error too small with 50 groups or fewer, for linear models fitted by restricted maximum likelihood. In limn's coverage study with 8 groups, intervals for a new unit's observation covered 0.87 to 0.90 at 90%. Reference: Maas and Hox (2005), Methodology 1, 86-92, pp. 86 and 89-90. Refit with reml = TRUE in illume::ilm_model(): restricted maximum likelihood is almost always at least as good as maximum likelihood, and sometimes better, especially at estimating variances (the reference, p. 89, citing Browne 1998). In limn's prototype study of few groups (studies/prototypes/) it raised 8-group coverage from 0.872 to 0.892 (Poisson) and from 0.891 to 0.903 (gaussian).
as.data.frame(p, level = c(0.5, 0.9))
#>   site x new_level estimate    mean lower_0.5 upper_0.5 lower_0.9 upper_0.9
#> 1    1 0     FALSE        2 1.76430         1         3         0         4
#> 2  new 0      TRUE        1 5.93636         1         3         0         5
summary(p)
#>   site x new_level estimate    mean lower_0.5 upper_0.5 lower_0.9 upper_0.9
#> 1    1 0     FALSE        2 1.76430         1         3         0         4
#> 2  new 0      TRUE        1 5.93636         1         3         0         5
```
