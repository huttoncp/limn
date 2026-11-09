# Cross-validated scores of a fitted illume model

Scores a model on data it was not fitted to, fold by fold: the model is
refitted to each fold's training rows
([`ilm_cv_refit()`](https://huttoncp.github.io/limn/reference/ilm_cv_refit.md),
which for an illume fit is
[`illume::ilm_refit()`](https://rdrr.io/pkg/illume/man/ilm_refit.html)),
predicts the fold's test rows, and each test row's predictive
distribution is scored against its outcome. Every row is scored by a
model that did not see it.

## Usage

``` r
ilm_cv(
  fit,
  data,
  folds,
  rules = c("log", "crps"),
  groups = NULL,
  uncertainty = c("joint", "fixed", "latent", "none"),
  refit = "each",
  ncores = 1L,
  nsim = 1000L,
  seed = 1L,
  ...
)
```

## Arguments

- fit:

  A fitted `"ilm_model"`.

- data:

  The data the folds index: the rows the model is refitted to and scored
  on.

- folds:

  Folds from
  [`ilm_folds()`](https://huttoncp.github.io/limn/reference/ilm_folds.md)
  built on `data`.

- rules:

  The scores, as losses (lower is better): any of `"log"`, `"crps"`,
  `"brier"` (binary outcomes), `"pinball"` (at `probs`), `"interval"`
  (for the central `level` interval), `"wis"` (the weighted interval
  score, over the central intervals at levels 1 - `alpha`), `"squared"`
  (about the predictive mean) and `"absolute"` (about the median). Each
  is computed from the predictive mixture.
  [`ilm_scores.ilm_pred()`](https://huttoncp.github.io/limn/reference/ilm_scores.ilm_pred.md)
  defines each rule and gives its source.

- groups:

  How the random effects enter the predictions, as in
  [`ilm_predict()`](https://huttoncp.github.io/limn/reference/ilm_predict.md);
  `NULL`, the default, follows the folds (see *Which unit a test row is
  predicted for*).

- uncertainty:

  As in
  [`ilm_predict()`](https://huttoncp.github.io/limn/reference/ilm_predict.md).

- refit:

  `"each"`: refit the model to every fold's training rows.

- ncores:

  Number of worker processes. The folds are shared out as workers become
  free; the results are the same whatever the number.

- nsim:

  Number of parameter draws for each fold's predictions.

- seed:

  Seed; each fold's draws use their own seed derived from it, so the
  results do not depend on `ncores`.

- ...:

  For the scores: `probs`, the pinball score's probabilities (default
  0.1, 0.5, 0.9); `level`, the interval score's coverage (default 0.9);
  and `alpha`, the weighted interval score's intervals (default 0.02,
  0.05, 0.1, ..., 0.9).

## Value

A data frame of class `"ilm_cv"`: one row per test row, with its fold,
its row in `data`, for rolling folds its origin and horizon, and one
column per score. A fold whose refit or prediction failed has `NA`
scores and its error in `failed`.
[`print()`](https://rdrr.io/r/base/print.html) gives the mean of each
score, by horizon for rolling folds; to compare two models' scores with
standard errors that respect the folds, use
[`ilm_compare()`](https://huttoncp.github.io/limn/reference/ilm_compare.md).

## Which unit a test row is predicted for

`groups` says how the random effects and the correlation over time
enter, as in
[`ilm_predict()`](https://huttoncp.github.io/limn/reference/ilm_predict.md).
By default it follows the folds:

- **grouped folds** hold out whole units, so a test row's unit is new to
  the training fit, and it is predicted as a new unit (`"new"`);

- **random, blocked and rolling folds** keep a unit's other rows in the
  training data, so it is predicted as itself (`"fitted"`).

Random or blocked folds on a model with a correlation over time score
interpolation, not forecasting: a test row sits between training rows of
its own series, so its neighbours' cells inform it. `ilm_cv()` warns,
and names rolling-origin folds, which score forecasts.

## Ex-ante and ex-post

Rolling-origin folds score forecasts. A covariate whose future value
would not have been known at the origin, such as the weather on the day
forecast, makes the evaluation *ex-post*: it flatters the model. The
folds record which covariates are known in advance
(`ilm_folds(known = )`); the result says which evaluation it is.

## Checks

What each check means and what to do about it, in the words it is
printed with. A reference supports what was found and why; a remedy
draws on it only where it says so, and anything marked as limn's own
reasoning is not the source's:

- **folds and time** (warn): The model has a correlation over time, and
  the folds are random or blocked. A test row then sits between training
  rows of its own series, whose cells inform it, so the scores are those
  of interpolation, not forecasting, and flatter the model's forecasts
  (limn's own reasoning). Reference: Hyndman and Athanasopoulos (2021),
  Forecasting: Principles and Practice, 3rd edition, section 5.10. To
  score forecasts, use rolling-origin folds, whose training rows all
  come before the rows they score (the reference): ilm_folds(type =
  "rolling", time = , group = , initial = , horizon = ).

- **cross-validation** (fail): None of the K folds could be refitted and
  scored. The first fold's error: (the first error) ilm_cv_refit()
  refits the model to each fold's training rows; for an illume fit it is
  illume::ilm_refit(), which repeats the model's own call. A model it
  cannot refit, or folds whose training rows cannot support it, leave
  nothing to score. Check that every fold's training rows can support
  the model, and fit it with a direct call to illume::ilm_model(): a fit
  made inside a function that passes its arguments on with ... can leave
  a call ilm_refit() cannot repeat.

## References

Hyndman, R. J. and Athanasopoulos, G. (2021). *Forecasting: Principles
and Practice*, 3rd edition. OTexts. Section 5.10, "Time series
cross-validation" (online edition,
<https://otexts.com/fpp3/toolbox.html#sec-tscv>, read 27 October 2026).

## See also

[`ilm_folds()`](https://huttoncp.github.io/limn/reference/ilm_folds.md),
[`ilm_compare()`](https://huttoncp.github.io/limn/reference/ilm_compare.md),
[`ilm_plot_horizon()`](https://huttoncp.github.io/limn/reference/ilm_plot_horizon.md),
[`ilm_scores.ilm_pred()`](https://huttoncp.github.io/limn/reference/ilm_scores.ilm_pred.md)
for the scoring rules, and
[`ilm_cv_refit()`](https://huttoncp.github.io/limn/reference/ilm_cv_refit.md)
for how a fit is refitted to each fold.

## Examples

``` r
set.seed(1)
d <- data.frame(site = factor(rep(1:12, each = 10)), t = rep(1:10, 12), x = rnorm(120))
d$y <- rpois(120, exp(0.5 + 0.3 * d$x + rnorm(12, 0, 0.4)[d$site]))
fit <- illume::ilm_model(y ~ x + (1 | site), data = d, family = "poisson",
                         verbose = FALSE)
folds <- ilm_folds(d, type = "grouped", group = "site", K = 4)
cv <- ilm_cv(fit, d, folds, nsim = 200)
cv
#> <ilm_cv> 120 scored rows over 4 grouped folds; refit each, groups = "new"
#> Mean scores (lower is better):
#>   log  crps 
#>  1.72 0.830 
#> For standard errors, compare models with ilm_compare().
```
