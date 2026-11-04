# Compare models' predictive scores, with standard errors that respect the data

Compares the scores of two or more models on the same rows: the mean
difference in a score, its standard error, and a test of whether it is
zero. The first model is the reference. Every other model is compared
with it, as that model's mean score minus the reference's. Scores are
losses, so a negative difference means the model predicts better than
the reference.

## Usage

``` r
ilm_compare(..., rule = NULL, by = "row", lag = NULL, method = NULL)
```

## Arguments

- ...:

  Two or more tables of scores for the same rows, in the same order,
  such as the output of the scoring of an `"ilm_pred"`. Name them
  (`ilm_compare(ar = s1, iid = s2)`) to label the comparison.

- rule:

  The score to compare: the name of a score column. By default the first
  score in the table, which for a plain data frame is its first column
  named after a rule (`log`, `crps`, `brier`, `pinball_*`, `interval_*`,
  `squared`, `absolute`, `mase` or `rmsse`).

- by:

  `"row"`, `"origin"`, or the name of a column the rows are clustered
  by, such as `"group"` or `"fold"`. See *Why the standard error depends
  on `by`*.

- lag:

  Only for `by = "origin"` with `method = "dm"`: the truncation lag of
  the long-run variance, that is, how many autocovariances of the
  per-origin mean differences it keeps.

  - **The default** suits origins one time step apart. Forecasts h steps
    ahead from neighbouring origins then overlap in h - 1 steps, so the
    pooled row uses the longest horizon minus one, and horizon h's own
    row uses h - 1.

  - **Set it yourself when the origins are further apart.** With origins
    s time steps apart, forecasts h ahead overlap in ceiling(h / s) - 1
    of them, so use `lag = ceiling(H / s) - 1` for a longest horizon H.
    It is 0 when the origins are at least H apart, and the forecasts do
    not overlap.

  - A value you give applies to the pooled row and to every horizon's
    row. It must be between 0 and the number of origins minus 1.

  - If the variance it gives is negative, it is taken as 0 and the null
    hypothesis rejected (statistic infinite, p-value 0), as Diebold and
    Mariano (1995, section 1.1) prescribe, and the test warns.

  With `method = "wpe"`, or any `by` other than `"origin"`, giving `lag`
  is an error. The fixed-smoothing test chooses its own bandwidth from
  the number of origins, whatever their spacing.

- method:

  Only for `by = "origin"`: which Diebold-Mariano test.

  - `"wpe"`, the default: the fixed-smoothing test of Coroneo and Iacone
    (2020, section 4.2). The variance of the T per-origin mean
    differences is their periodogram averaged over the first m =
    floor(T^(1/3)) Fourier frequencies (equation (8); the bandwidth of
    their simulations, section 5). The statistic is referred to t on 2m
    degrees of freedom (equation (9)). It needs at least 3 origins and
    takes no `lag`. In limn's simulation it held its size from 10
    origins up.

  - `"dm"`: the usual test. It uses a long-run variance truncated at
    `lag` with equal weights (Diebold and Mariano 1995, section 1.1),
    Harvey, Leybourne and Newbold's (1997) correction (equation (9), p.
    283), and t on T - 1 degrees of freedom (pp. 283-284). Use it to
    compare with published results. With fewer than about 100 origins
    over several horizons it rejects too often, and warns.

  Giving `method` with any `by` other than `"origin"` is an error. See
  *Why the standard error depends on `by`*.

## Value

A data frame of class `"ilm_compare"`, one row per comparison, with the
two mean scores, their difference, its standard error, the test
statistic, degrees of freedom and p-value, the number of rows and of
clusters, and the method. For `by = "origin"` with several horizons,
each comparison has a row pooled over the horizons (`horizon` is `NA`)
and a row for each horizon.

The test is also given as data, so nothing needs reading from `method`,
which describes it in words:

- `test`: `"independent"` (`by = "row"`), `"clustered"`, `"wpe"` or
  `"dm"` (`by = "origin"`).

- `bandwidth`: for `"wpe"`, the number m of Fourier frequencies
  averaged; otherwise `NA`.

- `lag`: for `"dm"`, the truncation lag used; otherwise `NA`.

## Why the standard error depends on `by`

Scores on rows from the same group, fold or forecast origin are
correlated, because those rows share whatever the model got right or
wrong about them. Treating them as independent understates the standard
error, and makes a difference look more certain than it is. So the
standard error is computed at the level the data are dependent at:

- `by = "row"`: rows are independent. Right for rows with no shared
  structure.

- `by` = a column such as `"group"` or `"fold"`: a cluster-robust
  standard error, with the rows of a cluster summed first. This is the
  same as
  [`sandwich::vcovCL()`](https://zeileis.codeberg.page/sandwich/reference/vcovCL.html)
  with its default HC1 adjustment, applied to the differences. The test
  uses a t distribution with one fewer degrees of freedom than there are
  clusters.

- `by = "origin"`: rolling-origin forecasts. The differences are
  averaged within each origin, and those origin means are tested with a
  Diebold-Mariano test: pooled over the horizons, and at each horizon on
  its own when there are several. Forecasts from neighbouring origins
  overlap when the horizon is longer than one, so the origin means are
  autocorrelated, and the test's variance allows for it.

  - `method = "wpe"` (the default): the fixed-smoothing test of Coroneo
    and Iacone (2020). The variance is the periodogram averaged over the
    first m Fourier frequencies of the T origin means (their weighted
    periodogram estimate with the Daniell kernel, section 4.2, equation
    (8)), and the statistic is referred to t on 2m degrees of freedom
    (equation (9)). m = floor(T^(1/3)) is the bandwidth they use in
    their simulations (section 5). In limn's simulation
    (`studies/findings/compare_dm_fixed.md`) it held its size from 10
    origins up, where the usual test rejected up to 12% of the time at
    5%.

  - `method = "dm"`: the usual test, a long-run variance truncated at
    `lag` with equal weights, the uniform lag window Diebold and Mariano
    (1995, section 1.1) recommend, with Harvey, Leybourne and
    Newbold's (1997) small-sample correction (equation (9), p. 283),
    referred to t on T - 1 degrees of freedom (pp. 283-284). When no
    fallback is needed it agrees with `forecast::dm.test()`, as limn's
    study `studies/scripts/compare_dm.R` checks; that package is an
    outside comparator, not a source. By default `lag` is the longest
    horizon minus one, the overlap when the origin moves one time point
    at a time. Diebold and Mariano take k - 1 as a benchmark for
    k-step-ahead forecasts, because optimal forecasts' errors are at
    most (k - 1)-dependent (section 1.1); they note that this can fail
    in practice. With fewer than about 100 origins over several horizons
    it rejects too often, and warns.

## Checks

What each check means and what to do about it, in the words it is
printed with. A reference supports what was found and why; a remedy
draws on it only where it says so, and anything marked as limn's own
reasoning is not the source's:

- **few origins** (warn): The comparison has T origins over several
  horizons, fewer than about 100. With so few origins the usual
  Diebold-Mariano test rejects too often when forecasts overlap: in
  limn's simulation, up to 12% of the time at 5% with 50 origins or
  fewer (studies/findings/compare_dm_fixed.md). Reference: Coroneo and
  Iacone (2020), Journal of Applied Econometrics 35(4), 391-409, section
  4.2 (the test) and section 5 (its bandwidth). Use method = "wpe", the
  reference's fixed-smoothing test, which held its size from 10 origins
  up in the same simulation; or compare over more origins.

- **Diebold-Mariano variance** (warn): The long-run variance truncated
  at lag L was negative, so it was taken as 0 and the null hypothesis of
  equal accuracy was rejected. A long-run variance truncated with equal
  weights can come out negative. Diebold and Mariano's rule for that
  case is to take it as 0 and reject the null, so the statistic is
  infinite and the p-value 0. They call a negative estimate a rare
  event. In limn's own check against the simulation design of Harvey,
  Leybourne and Newbold (1997, p. 284), it arose in 35% of samples with
  8 origins at horizon 5 (lag 4), and in 17% with 16, and the rule took
  the modified test's size at nominal 10% to 41.9% there, and to 29.7%
  with 8 origins at horizon 3 (lag 2), against their Table 1's 12.9% and
  18.1% (p. 285) (studies/findings/compare_dm_negative.md). Reference:
  Diebold and Mariano (1995), Journal of Business and Economic
  Statistics 13, 253-263, section 1.1. Read this rejection with care, or
  use method = "wpe", the default, which needs no lag and held its size
  from 10 origins up in limn's simulation.

## Catching the checks

Each check above is raised as a warning of class `lmn_<id>`, beside
`lmn_diagnostic`, and is kept on the result as a row of
`attr(x, "checks")`, a data frame whose `id` column names it:

- `dm_few_origins` (warning class `lmn_dm_few_origins`): few origins
  over several horizons for the usual test;

- `dm_negative_variance` (warning class `lmn_dm_negative_variance`): the
  truncated long-run variance was negative, so it was taken as 0 and the
  null hypothesis rejected.

So
`withCallingHandlers(ilm_compare(...), lmn_dm_negative_variance = function(w) ...)`
handles that warning alone, and `attr(cmp, "checks")$id` says which
checks a comparison raised.

## References

Diebold, F. X. and Mariano, R. S. (1995). Comparing predictive accuracy.
*Journal of Business and Economic Statistics* 13, 253-263. Section 1.1:
the statistic, the truncated long-run variance with the uniform
(rectangular) lag window they recommend, the lag k - 1 as a benchmark
for k-step-ahead forecasts, and the rule for a negative estimate.

Coroneo, L. and Iacone, F. (2020). Comparing predictive accuracy in
small samples using fixed-smoothing asymptotics. *Journal of Applied
Econometrics* 35(4), 391-409. Section 4.2, equations (8) and (9): the
weighted periodogram estimate and its t limit on 2m degrees of freedom;
section 5: the bandwidth m = floor(T^(1/3)).

Harvey, D., Leybourne, S. and Newbold, P. (1997). Testing the equality
of prediction mean squared errors. *International Journal of
Forecasting* 13, 281-291. Equation (9), p. 283: the corrected statistic;
pp. 283-284: Student's t on n - 1 degrees of freedom.

## See also

[`ilm_folds()`](https://huttoncp.github.io/limn/reference/ilm_folds.md)
for folds that make scores comparable.

## Examples

``` r
set.seed(1)
g <- rep(1:10, each = 5)
a <- data.frame(group = g, log = rnorm(50, 1.2) + rnorm(10)[g])
b <- data.frame(group = g, log = a$log - 0.1 + rnorm(50, 0, 0.2))
ilm_compare(reference = a, candidate = b, by = "group")
#> candidate's log score is lower than reference's by 0.0684 a row (SE 0.0354, clustered by group, 10 clusters; p = 0.0854).
#> Scores are losses: lower is better.
```
