# Proper scores of a prediction against observed outcomes

Scores each row of a prediction against its observed outcome, by proper
scoring rules computed from the predictive distribution itself, and by
the errors of the point forecasts it gives. Every score is a loss: lower
is better.

## Usage

``` r
# S3 method for class 'ilm_pred'
ilm_scores(
  object,
  y,
  rules = c("log", "crps"),
  probs = c(0.1, 0.5, 0.9),
  level = 0.9,
  naive = NULL,
  history = NULL,
  period = 1L,
  alpha = c(0.02, 0.05, 1:9/10),
  ...
)
```

## Arguments

- object:

  A prediction from
  [`ilm_predict()`](https://huttoncp.github.io/limn/reference/ilm_predict.md)
  or
  [`ilm_forecast()`](https://huttoncp.github.io/limn/reference/ilm_forecast.md)
  with `target = "response"`.

- y:

  The observed outcomes, one per row of the prediction.

- rules:

  The scoring rules: any of `"log"`, `"crps"`, `"brier"`, `"pinball"`,
  `"interval"`, `"wis"`, `"squared"`, `"absolute"`, `"mase"` and
  `"rmsse"`.

- probs:

  The probabilities of the pinball loss.

- level:

  The coverage of the interval score's central interval.

- naive:

  For `mase` and `rmsse`: each series' scale, from
  [`ilm_naive_scale()`](https://huttoncp.github.io/limn/reference/ilm_naive_scale.md).
  See *Giving MASE and RMSSE their scale*.

- history:

  For `mase` and `rmsse`: the data the model was trained on, from which
  each series' scale is computed. See *Giving MASE and RMSSE their
  scale*.

- period:

  With `history`, the naive forecast's lag: 1 for the naive forecast,
  the season's length for the seasonal naive.

- alpha:

  For `"wis"`: the central intervals it weighs, each at level 1 -
  `alpha`, with `alpha` strictly between 0 and 1. The default is Bracher
  et al. (2021)'s set, 0.02, 0.05 and 0.1 to 0.9 by 0.1.

- ...:

  Unused.

## Value

A data frame: the prediction's rows, with one column per score: `log`,
`crps`, `brier`, `pinball_<p>` for each p, `interval_<level>`, `wis`,
`squared`, `absolute`, `mase` and `rmsse`, as asked for. Attribute
`rules` names the score columns, and `orientation` is "lower is better".
For a discrete outcome whose CRPS was summed approximately, attributes
`crps_bound` and `crps_bracket` give each row's error estimate and its
rigorous bound.

## The scoring rules

For a row with outcome `y`, predictive CDF `F`, density or probability
mass `f`, mean `m` and median `med`, written as R expressions (a
comparison such as `(y < q)` counts as 1 when true and 0 when false):

- `"log"`: `-log(f(y))`, the negative of the logarithmic score (Gneiting
  and Raftery 2007, section 4.1, equation (19), p. 365; for a discrete
  outcome, the probability mass, section 3.1, Example 3, p. 363).

- `"crps"`: the integral over `z` of `(F(z) - (z >= y))^2`, the
  continuous ranked probability score in its negative orientation
  (section 4.2, equation (20), p. 367). For a count or an ordered
  outcome, `F` steps only at whole numbers, so the integral is the sum
  over them of `(F(k) - (k >= y))^2`.

- `"brier"`: `(P(Y = 1) - y)^2`, for a binary outcome coded 0/1: the
  Brier score for a dichotomous event (Table 1, p. 365), as a loss.

- `"pinball"`: `(y - q) * (p - (y < q))`, with `q` the predictive
  `p`-quantile, for each `p` in `probs`: the negative of the quantile
  score (section 6.1, equation (41), p. 370), the "tick" or "check"
  loss.

- `"interval"`:
  `(u - l) + (2 / alpha) * (l - y) * (y < l) + (2 / alpha) * (y - u) * (y > u)`,
  for the central interval from `l` to `u` at `level` = 1 - `alpha`: the
  interval score (section 6.2, equation (43), p. 370).

- `"wis"`: `(abs(y - med) / 2 + sum(alpha / 2 * IS)) / (K + 1/2)`, with
  `IS` the interval score above for the central interval at level 1 -
  `alpha`, for each of the `K` values of `alpha`: the weighted interval
  score of Bracher et al. (2021, equations (1) and (2), p. 4). The
  default `alpha` is theirs: 0.02, 0.05, 0.1, 0.2, ..., 0.9, the median
  and 11 central intervals from 98% to 10%, 23 quantiles in all (p. 4).
  It is a weighted sum of quantile scores, so it is proper, and it needs
  only the prediction's quantiles. With many evenly spaced levels it
  approaches the CRPS (equation (3), p. 4), and with
  `alpha = numeric(0)` it is the absolute error (p. 4). limn's tests
  reproduce their worked example (section 3.2, p. 7) with the values of
  its correction (Bracher et al. 2022, p. 1).

- `"squared"`: `(y - m)^2`. The energy score with exponent 2 reduces to
  this, in the opposite orientation (section 4.3, equation (23), p.
  367), a proper but not strictly proper score.

- `"absolute"`: `abs(y - med)`, twice the pinball loss at `p` = 1/2.

- `"mase"`: `abs(y - med) / s`, with `s` the series' in-sample mean
  absolute naive error at lag `period`. The mean over rows is the mean
  absolute scaled error.

- `"rmsse"`: `(y - m)^2 / s2`, with `s2` the series' in-sample mean
  squared naive error. The square root of the mean over rows is the root
  mean squared scaled error.

The scaled errors and their in-sample naive scales are as fpp3 defines
them (section 5.8, "Scaled errors"), where they are credited to Hyndman
and Koehler (2006). fpp3 scales the error of a point forecast; `limn`
takes the median as the point forecast for the absolute error and the
mean for the squared error, the point forecasts each error rewards (its
own choice).

## Giving MASE and RMSSE their scale

The scaled errors (`mase`, `rmsse`) divide each row's error by its
series' in-sample naive-forecast error, so they need each series'
history. There are two ways to give it, and they give identical scores
wherever both apply:

- **`history = train`**, the common case: the data the model was trained
  on. The prediction already knows its time, series and response columns
  (it records them), so nothing else is needed. `period` sets the naive
  forecast's lag: 1, or the season's length for the seasonal naive.

- **`naive = ilm_naive_scale(train, time =, group =, response =)`**: the
  scale computed once and passed in. Use it when the prediction doesn't
  know its time column (see *Time and series columns* in
  [`ilm_predict()`](https://huttoncp.github.io/limn/reference/ilm_predict.md)),
  or when several models' predictions must share one scale, so their
  scaled errors compare on the same footing.

When `history =` can't work, because the prediction has no time column,
the error says so and gives the `naive =` call to use instead. Give one
of the two, not both.

## References

Bracher, J., Ray, E. L., Gneiting, T. and Reich, N. G. (2021).
Evaluating epidemic forecasts in an interval format. *PLOS Computational
Biology* 17(2), e1008618. Section 2.2, equations (1) to (4) (pp. 4-5);
section 3.2 (p. 7).

Bracher, J., Ray, E. L., Gneiting, T. and Reich, N. G. (2022).
Correction: Evaluating epidemic forecasts in an interval format. *PLOS
Computational Biology* 18(10), e1010592. P. 1.

Gneiting, T. and Raftery, A. E. (2007). Strictly proper scoring rules,
prediction, and estimation. *Journal of the American Statistical
Association* 102(477), 359-378. Section 3.1, Example 3 (p. 363); Table 1
(p. 365); section 4.1, equation (19) (p. 365); section 4.2, equation
(20) and section 4.3, equation (23) (p. 367); section 6.1, equation
(41), and section 6.2, equation (43) (p. 370).

Hyndman, R. J. and Athanasopoulos, G. (2021). *Forecasting: Principles
and Practice*, 3rd edition. OTexts. Section 5.8, "Evaluating point
forecast accuracy" (online edition,
<https://otexts.com/fpp3/toolbox.html#sec-accuracy>, read 27 October
2026).

Hyndman, R. J. and Koehler, A. B. (2006). Another look at measures of
forecast accuracy. *International Journal of Forecasting* 22, 679-688.

## See also

[`ilm_naive_scale()`](https://huttoncp.github.io/limn/reference/ilm_naive_scale.md)
for the scaled errors' scale;
[`ilm_cv()`](https://huttoncp.github.io/limn/reference/ilm_cv.md) for
scores out of sample;
[`ilm_compare()`](https://huttoncp.github.io/limn/reference/ilm_compare.md)
to compare two models' scores;
[`ilm_pit()`](https://huttoncp.github.io/limn/reference/ilm_pit.md) and
[ilm_calibration()](https://huttoncp.github.io/limn/reference/ilm_calibration.ilm_pred.md)
for calibration.

## Examples

``` r
set.seed(1)
d <- data.frame(site = factor(rep(c("a", "b", "c"), each = 16)), t = rep(1:16, 3))
d$y <- rpois(nrow(d), exp(1.5 + ave(rnorm(nrow(d), 0, 0.15), d$site, FUN = cumsum)))
train <- d[d$t <= 12, ]; test <- d[d$t > 12, ]
fit <- illume::ilm_model(y ~ 1, data = train, family = "poisson",
                         ar = illume::ilm_rw1(~ t | site), verbose = FALSE)
fc <- ilm_forecast(fit, test[c("site", "t")], nsim = 200)
## the scale read from the training data, through the prediction's own columns
s <- ilm_scores(fc, test$y, rules = c("crps", "mase"), history = train)
colMeans(s[attr(s, "rules")])
#>     crps     mase 
#> 1.320875 0.727273 
## the same scale computed once, as several models would share it
sc <- ilm_naive_scale(train, time = "t", group = "site", response = "y")
all.equal(s$mase, ilm_scores(fc, test$y, rules = "mase", naive = sc)$mase)
#> [1] TRUE
```
