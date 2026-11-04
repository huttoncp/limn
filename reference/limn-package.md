# limn: A Projection Engine for Distributional Forecasting

Predictive distributions, forecasts, proper scoring rules, calibration,
cross-validation and model comparison for models fitted by
[`illume::ilm_model()`](https://rdrr.io/pkg/illume/man/ilm_model.html).
`limn` is the third package in a family, after `illumex` for exploration
and `illume` for inference, and shares their `ilm_` prefix, so one
prefix finds a function in any of the three.

## A prediction is a distribution

Every prediction `limn` makes is a predictive distribution: a mixture,
over draws of the fitted parameters, of the model family's distribution
for the response. Intervals, quantiles and scores are all read from that
mixture. Scoring it directly, rather than scoring responses simulated
from it, removes the simulation's noise from the score.

## Scores are losses

Every score `limn` reports is a loss, so lower is better: the log score
is the negative log predictive density, and the CRPS, Brier, pinball,
interval, squared and absolute scores are penalties.

## The workflow

- **Predict.**
  [`ilm_predict()`](https://huttoncp.github.io/limn/reference/ilm_predict.md)
  gives the distribution of a mean, a new observation or the linear
  predictor, for a unit the fit has seen, a new one, a typical one or
  the average over units.
  [`ilm_forecast()`](https://huttoncp.github.io/limn/reference/ilm_forecast.md)
  carries series forward in time.
  [`ilm_future()`](https://huttoncp.github.io/limn/reference/ilm_future.md)
  builds the rows a forecast needs, and
  [`ilm_covariate_paths()`](https://huttoncp.github.io/limn/reference/ilm_covariate_paths.md)
  and
  [`ilm_as_draws()`](https://huttoncp.github.io/limn/reference/ilm_as_draws.md)
  give covariates that are not known in advance as draws.

- **Combine.**
  [`ilm_aggregate()`](https://huttoncp.github.io/limn/reference/ilm_aggregate.md)
  totals a prediction over rows, and
  [`ilm_average()`](https://huttoncp.github.io/limn/reference/ilm_average.md)
  takes weighted averages of means (standardisation and
  post-stratification).

- **Check against outcomes.**
  [ilm_scores()](https://huttoncp.github.io/limn/reference/ilm_scores.ilm_pred.md)
  scores predictions by proper scoring rules;
  [`ilm_pit()`](https://huttoncp.github.io/limn/reference/ilm_pit.md)
  and
  [`ilm_plot_pit()`](https://huttoncp.github.io/limn/reference/ilm_plot_pit.md)
  show where outcomes fell in their distributions; and
  [ilm_calibration()](https://huttoncp.github.io/limn/reference/ilm_calibration.ilm_pred.md)
  gives a verdict on calibration, with
  [`ilm_plot_calibration()`](https://huttoncp.github.io/limn/reference/ilm_plot_calibration.md).

- **Validate out of sample.**
  [`ilm_folds()`](https://huttoncp.github.io/limn/reference/ilm_folds.md)
  builds folds that follow the data's structure, and
  [`ilm_cv()`](https://huttoncp.github.io/limn/reference/ilm_cv.md)
  scores a model on them.

- **Compare models.**
  [`ilm_compare()`](https://huttoncp.github.io/limn/reference/ilm_compare.md)
  compares scores with standard errors that respect how the data hang
  together;
  [`ilm_plot_compare()`](https://huttoncp.github.io/limn/reference/ilm_plot_compare.md)
  and
  [`ilm_plot_horizon()`](https://huttoncp.github.io/limn/reference/ilm_plot_horizon.md)
  draw the comparison.

- **Draw.**
  [`ilm_plot_fan()`](https://huttoncp.github.io/limn/reference/ilm_plot_fan.md)
  draws forecasts as fan charts.

- **Extend.**
  [`ilm_predict_draws()`](https://huttoncp.github.io/limn/reference/ilm_predict_draws.md)
  and
  [`ilm_cv_refit()`](https://huttoncp.github.io/limn/reference/ilm_cv_refit.md)
  let a package that fits illume's models another way supply the draws
  and the refits, so that every function above works on its fits.

## Checks, and what to do about them

Every check `limn` makes takes its words from one registry. Each says
what was found, what it means, and what to do about it, and gives the
source it rests on, or says that the reasoning is `limn`'s own. A check
is warned when it is raised, with class `lmn_<id>`, and kept on the
result in `checks`, so it can be read as data as well as words. Where
the fix is a change to the model, the check names the call:
[`illume::ilm_remedies()`](https://rdrr.io/pkg/illume/man/ilm_remedies.html)
lists the changes the fit's own checks support, and
[`illume::ilm_apply_remedy()`](https://rdrr.io/pkg/illume/man/ilm_apply_remedy.html)
refits the model with the one you choose.

## See also

Useful links:

- <https://github.com/huttoncp/limn>

- <https://huttoncp.github.io/limn/>

- Report bugs at <https://github.com/huttoncp/limn/issues>

## Author

**Maintainer**: Craig P. Hutton <craig.hutton@gmail.com> \[copyright
holder\]

## Examples

``` r
set.seed(1)
d <- data.frame(site = factor(rep(1:20, each = 8)), x = rnorm(160))
d$y <- rpois(160, exp(0.5 + 0.4 * d$x + rnorm(20, 0, 0.4)[d$site]))
train <- d[as.integer(d$site) <= 15, ]; test <- d[as.integer(d$site) > 15, ]
fit <- illume::ilm_model(y ~ x + (1 | site), data = train, family = "poisson",
                         verbose = FALSE)
## the held-out sites are new units: predict their outcomes, then score
## and check the predictions against what happened
p <- ilm_predict(fit, test, target = "response", groups = "new", nsim = 300)
#> Warning: 'site' has 15 levels, fewer than about 20. With so few groups the variance between them is poorly determined: the reference found its standard error too small with 50 groups or fewer, for linear models fitted by restricted maximum likelihood. In limn's coverage study with 8 groups, intervals for a new unit's observation covered 0.87 to 0.90 at 90%. Reference: Maas and Hox (2005), Methodology 1, 86-92, pp. 86 and 89-90. Refit with reml = TRUE in illume::ilm_model(): restricted maximum likelihood is almost always at least as good as maximum likelihood, and sometimes better, especially at estimating variances (the reference, p. 89, citing Browne 1998). In limn's prototype study of few groups (studies/prototypes/) it raised 8-group coverage from 0.872 to 0.892 (Poisson) and from 0.891 to 0.903 (gaussian).
head(as.data.frame(p))
#>     site          x y new_level estimate     mean lower_0.9 upper_0.9
#> 121   16 -0.5059575 2     FALSE        2 1.840419         0         5
#> 122   16  1.3430388 1     FALSE        3 3.439174         0         8
#> 123   16 -0.2145794 3     FALSE        2 2.028546         0         5
#> 124   16 -0.1795565 4     FALSE        2 2.052479         0         5
#> 125   16 -0.1001907 2     FALSE        2 2.107815         0         5
#> 126   16  0.7126663 2     FALSE        2 2.773334         0         7
s <- illume::ilm_scores(p, test$y)
colMeans(s[attr(s, "rules")])
#>       log      crps 
#> 1.9038094 0.9510902 
illume::ilm_calibration(p, test$y)
#> <calibration> 40 outcomes, poisson predictive distributions
#>   PIT uniformity: chi-square 22.4 on 7 df, p = 0.00217
#>   coverage by level (from the PITs):
#>      50%:  35.0% (SE 7.91 percentage points), mean width 2.10
#>      80%:  82.5% (SE 6.32 percentage points), mean width 4.13
#>      90%:  90.0% (SE 4.74 percentage points), mean width 5.28
#>      95%:  97.5% (SE 3.45 percentage points), mean width 6.30
#>   Consistent with calibration: coverage within sampling error at every level, no bias in the mean, and the PITs consistent with uniform (p = 0.00217).
```
