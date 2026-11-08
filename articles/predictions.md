# Predictions: which unit, and how uncertain

A mixed model can answer several different questions about the same row
of data. How many visits will a 50-year-old make at clinic C01? At a
clinic not in the data? At a typical clinic? On average across clinics?
[`ilm_predict()`](https://huttoncp.github.io/limn/reference/ilm_predict.md)
makes you name the question, because the answers differ. This vignette
works through the choices: which unit, whether the mean or a new
observation, and how much of the fit’s uncertainty to carry, and then
the check limn makes when there are few groups, with its remedy.

``` r

library(limn)
#> Loading required package: illume
#> Loading required package: illumex
#> Registered S3 methods overwritten by 'illume':
#>   method             from   
#>   c.ilm_remedies     illumex
#>   print.ilm_remedies illumex
#> 
#> Attaching package: 'illume'
#> The following objects are masked from 'package:illumex':
#> 
#>     ilm_apply_remedy, ilm_remedies, ilm_remedy_table, iml_apply_remedy,
#>     iml_remedies, iml_remedy_table
```

## The data and the model

Patients’ yearly visits to a clinic depend on their age, and clinics
differ. The simulation below makes data for any number of clinics, 25
patients each, with the clinics’ effects drawn from a normal
distribution on the log scale.

``` r

simulate_clinics <- function(n_clinics, seed) {
  set.seed(seed)
  d <- data.frame(clinic = factor(rep(sprintf("C%02d", 1:n_clinics), each = 25)))
  d$age <- round(runif(nrow(d), 20, 80))
  clinic_effect <- rnorm(n_clinics, 0, 0.5)
  d$visits <- rpois(nrow(d), exp(0.3 + 0.02 * (d$age - 50) + clinic_effect[d$clinic]))
  d
}
clinics30 <- simulate_clinics(30, seed = 2026)
fit <- illume::ilm_model(visits ~ I(age - 50) + (1 | clinic), data = clinics30,
                         family = "poisson", verbose = FALSE)
```

## Which unit

The same row, a 50-year-old at clinic C01, predicted four ways:

``` r

row <- data.frame(clinic = "C01", age = 50)
which_unit <- lapply(c("fitted", "new", "typical", "population"), function(g)
  cbind(groups = g, as.data.frame(ilm_predict(fit, row, groups = g, nsim = 1000))))
do.call(rbind, which_unit)
#>       groups clinic age new_level estimate     mean lower_0.9 upper_0.9
#> 1     fitted    C01  50     FALSE 1.332570 1.359288 1.0089211  1.796362
#> 2        new    C01  50     FALSE 1.252211 1.449515 0.5337022  2.976385
#> 3    typical    C01  50     FALSE 1.256268 1.258428 1.0555045  1.475182
#> 4 population    C01  50     FALSE 1.436300 1.447097 1.2004809  1.713918
```

These are predictions of the mean number of visits, and each answers a
different question:

- `"fitted"`: clinic C01 itself, with its own estimated effect.
- `"new"`: a clinic not in the data. Its effect is unknown, so it is
  drawn from the fitted distribution of clinics, and the interval is
  much the widest. Its `estimate` is the predictive median, and `mean`
  the expected value.
- `"typical"`: a clinic whose effect is zero.
- `"population"`: the mean averaged over clinics.

The typical clinic and the population average differ. With a log link,
the mean at the middle clinic is not the average of the clinics’ means:
here the average is higher. It agrees, up to simulation error, with the
new clinic’s expected value, since a new clinic is a clinic drawn at
random. Which of the four to report depends on the question, so a model
with random effects has no default `groups`.

## The mean, or a new observation

`target = "mean"` predicts the expected number of visits.
`target = "response"` predicts one patient’s count, which adds the
Poisson variation around that mean:

``` r

new_clinic <- data.frame(clinic = "new", age = 50)
rbind(cbind(target = "mean",
            as.data.frame(ilm_predict(fit, new_clinic, target = "mean", groups = "new", nsim = 1000))),
      cbind(target = "response",
            as.data.frame(ilm_predict(fit, new_clinic, target = "response", groups = "new", nsim = 1000))))
#>     target clinic age new_level estimate     mean lower_0.9 upper_0.9
#> 1     mean    new  50     FALSE 1.252211 1.449515 0.5337022  2.976385
#> 2 response    new  50     FALSE 1.000000 1.449515 0.0000000  4.000000
```

The mean is the same either way; the interval for one patient is wider,
and in whole visits.

## Averages over a population

A model’s predictions can be averaged over a population of rows, draw by
draw, so the average keeps the uncertainty its rows share.
[`ilm_average()`](https://huttoncp.github.io/limn/reference/ilm_average.md)
does this. Averaged over this sample’s own patients, each at their own
clinic, the model gives:

``` r

p_all <- ilm_predict(fit, clinics30, groups = "fitted", nsim = 500)
as.data.frame(ilm_average(p_all))
#>   n_rows weight estimate     mean lower_0.9 upper_0.9
#> 1    750    750 1.479423 1.496225   1.42749  1.564945
```

With population counts as weights, the same average is
post-stratification. Gelman and Hill (2007, section 14.1, equation
(14.1)) weight each cell’s prediction by its count in the population.
Here the cells are three ages, averaged over clinics, in a population
with more young patients than the sample:

``` r

cells <- data.frame(age = c(30, 50, 70), count = c(500, 300, 200))
p_cells <- ilm_predict(fit, cells, groups = "population", nsim = 500)
as.data.frame(p_cells)
#>   age count new_level estimate     mean lower_0.9 upper_0.9
#> 1  30   500     FALSE 1.007052 1.019340 0.8441404  1.230015
#> 2  50   300     FALSE 1.436300 1.451943 1.2150372  1.732464
#> 3  70   200     FALSE 2.048512 2.071122 1.7188660  2.481511
as.data.frame(ilm_average(p_cells, weights = "count"))
#>   n_rows weight estimate     mean lower_0.9 upper_0.9
#> 1      3   1000 1.344118 1.359477  1.136654  1.625557
```

Half the weight is on the youngest cell, so the average sits below the
middle cell’s prediction.

## How uncertain

A prediction carries the fit’s uncertainty by drawing the parameters.
The `uncertainty` argument says how much is drawn:

- `"joint"` (the default): every parameter, the variances included;
- `"fixed"`: the variances held at their estimates, everything else
  drawn;
- `"latent"`: every parameter held, and only the random effects drawn;
- `"none"`: the plug-in prediction at the estimates.

With 30 clinics the four choices give nearly the same interval for a new
clinic (not shown). With few clinics they don’t, so the rest of this
vignette uses 8. The model is fitted by maximum likelihood
(`reml = FALSE`), to be compared with restricted maximum likelihood
below.

``` r

clinics8 <- simulate_clinics(8, seed = 2026)
fit8 <- illume::ilm_model(visits ~ I(age - 50) + (1 | clinic), data = clinics8,
                          family = "poisson", reml = FALSE, verbose = FALSE)
```

The 90% interval for a new clinic’s mean, under each choice:

``` r

widths <- sapply(c("joint", "fixed", "latent", "none"), function(u) {
  p <- suppressWarnings(ilm_predict(fit8, new_clinic, groups = "new", uncertainty = u, nsim = 2000))
  a <- as.data.frame(p)
  c(lower = a$lower_0.9, upper = a$upper_0.9, width = a$upper_0.9 - a$lower_0.9)
})
round(widths, 2)
#>       joint fixed latent none
#> lower  0.78  0.81   0.85 0.85
#> upper  2.78  2.68   2.62 2.62
#> width  2.00  1.87   1.77 1.77
```

Here the more that is drawn, the wider the interval. For a new clinic,
`"latent"` and `"none"` agree: both hold the parameters, and the new
clinic’s effect is drawn either way. In limn’s prototype study of few
groups (`studies/prototypes/`), `"joint"` came closest to nominal
coverage for a new unit’s observation, though still below it: 0.88 for
90% intervals with 8 groups, against 0.85 to 0.87 for the other methods
it tried, each with a simulation standard error of about 0.01.
([`suppressWarnings()`](https://rdrr.io/r/base/warning.html) hides a
warning here, the few-groups check for means; the next section turns to
the few-groups checks.)

## The few-groups check, and its remedy

With so few clinics, the variance between them is poorly determined, and
intervals for a new clinic cover less often than they say.
[`ilm_predict()`](https://huttoncp.github.io/limn/reference/ilm_predict.md)
warns, says why, and names the remedy:

``` r

p8 <- ilm_predict(fit8, new_clinic, target = "response", groups = "new", nsim = 1000)
#> Warning: 'clinic' has 8 levels, fewer than about 20. With so few groups the
#> variance between them is poorly determined: the reference found its standard
#> error too small with 50 groups or fewer, for linear models fitted by restricted
#> maximum likelihood. In limn's coverage study with 8 groups, intervals for a new
#> unit's observation covered 0.87 to 0.90 at 90%. Reference: Maas and Hox (2005),
#> Methodology 1, 86-92, pp. 86 and 89-90. Refit with reml = TRUE in
#> illume::ilm_model(): restricted maximum likelihood is almost always at least as
#> good as maximum likelihood, and sometimes better, especially at estimating
#> variances (the reference, p. 89, citing Browne 1998). In limn's prototype study
#> of few groups (studies/prototypes/) it raised 8-group coverage from 0.872 to
#> 0.892 (Poisson) and from 0.891 to 0.903 (gaussian).
```

Maas and Hox (2005, pp. 86 and 89-90) found the standard errors of the
variances between groups too small with 50 groups or fewer, in linear
multilevel models. The remedy the warning names is to refit by
restricted maximum likelihood:

``` r

fit8_reml <- illume::ilm_model(visits ~ I(age - 50) + (1 | clinic), data = clinics8,
                               family = "poisson", reml = TRUE, verbose = FALSE)
p8_reml <- ilm_predict(fit8_reml, new_clinic, target = "response", groups = "new", nsim = 1000)
#> Warning: 'clinic' has 8 levels, fewer than about 20. With so few groups the
#> variance between them is poorly determined: the reference found its standard
#> error too small with 50 groups or fewer, for linear models fitted by restricted
#> maximum likelihood. In limn's coverage study with 8 groups, intervals for a new
#> unit's observation covered 0.87 to 0.90 at 90%. Reference: Maas and Hox (2005),
#> Methodology 1, 86-92, pp. 86 and 89-90. The fit already used restricted maximum
#> likelihood (REML), the refit this check otherwise suggests. In limn's prototype
#> study of few groups (studies/prototypes/), with REML and 8 groups of 6, 90%
#> intervals for a new unit's observation covered 0.892 (Poisson) and 0.903
#> (gaussian), with simulation standard errors of 0.004 to 0.005. More groups is
#> the sure cure.
```

The check still warns, since 8 clinics are still few, but it no longer
offers the refit: the fit already used it. It says what limn’s study
found with it, and that more groups is the sure cure.

For a mean rather than a new observation, a related check warns that the
intervals are somewhat too narrow, and that no remedy for them has been
measured yet; see the Checks section of
[`?ilm_predict`](https://huttoncp.github.io/limn/reference/ilm_predict.md).
With the 30 clinics above, neither check warns.

## Reading a prediction

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) gives
each row’s estimate, mean and 90% interval;
[`summary()`](https://rdrr.io/r/base/summary.html) gives the 50% and 90%
intervals. Their `new_level` column is `TRUE` only for a row naming a
clinic the fit has not seen, predicted under `groups = "fitted"` as a
new one.

``` r

summary(ilm_predict(fit, new_clinic, groups = "new", nsim = 1000))
#>   clinic age new_level estimate     mean lower_0.5 upper_0.5 lower_0.9
#> 1    new  50     FALSE 1.252211 1.449515 0.8891517  1.752097 0.5337022
#>   upper_0.9
#> 1  2.976385
```

## Where next

- [`?ilm_predict`](https://huttoncp.github.io/limn/reference/ilm_predict.md)
  for every choice, the prediction object, and the checks.
- [`?ilm_average`](https://huttoncp.github.io/limn/reference/ilm_average.md)
  for standardisation and post-stratification.
- [`vignette("limn")`](https://huttoncp.github.io/limn/articles/limn.md)
  for forecasting, scoring and comparing models.

## References

Gelman, A. and Hill, J. (2007). *Data Analysis Using Regression and
Multilevel/Hierarchical Models*. Cambridge University Press. Section
14.1.

Maas, C. J. M. and Hox, J. J. (2005). Sufficient sample sizes for
multilevel modeling. *Methodology* 1(3), 86-92. Pp. 86 and 89-90.
