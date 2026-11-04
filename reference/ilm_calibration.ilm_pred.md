# Calibration of a prediction against observed outcomes

Whether outcomes fell where their predictive distributions said they
would: coverage and width level by level, bias in the mean, and the
uniformity of the probability integral transforms (PITs; Gneiting,
Balabdaoui and Raftery 2007, section 3.1, pp. 251-253). A binary outcome
is judged instead by its calibration intercept and slope, which together
assess what Van Calster et al. (2019, p. 4) call weak calibration. The
result ends with a verdict in plain words and, where something is wrong,
a remedy that exists, most often in `illume`.

## Usage

``` r
# S3 method for class 'ilm_pred'
ilm_calibration(
  object,
  y,
  bins = 10L,
  levels = c(0.5, 0.8, 0.9, 0.95),
  seed = 1L,
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

  The observed outcomes, one per row of the prediction. At least 20 must
  be scoreable.

- bins:

  The number of bins for the PIT's uniformity test, and for a binary
  outcome's reliability table. There are fewer when the outcomes are
  few.

- levels:

  The central interval levels at which coverage is checked.

- seed:

  The seed for the randomised PITs of a discrete outcome. The caller's
  random-number state is left as it was.

- ...:

  Unused.

## Value

An object of class `"ilm_pred_calibration"`, a list:

- `pit`, each row's PIT;

- `n`, the outcomes used;

- `histogram` and `uniformity`;

- `bias`;

- `coverage`, a data frame of level, coverage, its standard error and
  mean width;

- `family`;

- for a binary outcome, `binary` (the intercept and slope) and
  `reliability`;

- `checks`, the diagnostics as data;

- `verdict`, their words as printed.

[`ilm_plot_calibration()`](https://huttoncp.github.io/limn/reference/ilm_plot_calibration.md)
draws it.

## Details

This is limn's method for `illume`'s
[`illume::ilm_calibration()`](https://rdrr.io/pkg/illume/man/ilm_calibration.html)
generic, for a prediction. `illume`'s own method checks a fitted model's
category probabilities.

## What is checked

- **Coverage and width.** At each level in `levels`, the share of
  outcomes inside the central interval of their predictive
  distributions, with its standard error, and the intervals' mean width.
  Coverage is read from the PITs, not from whether the outcome lies
  inside a quantile interval. A discrete outcome's quantile interval
  covers more than its level, because of ties, while the randomised
  PIT's coverage is exactly nominal when the model is right.

- **Bias**, on the outcome's own scale: the mean of the outcome minus
  its predictive mean, with its standard error. The mean PIT is not
  used, since a distribution with the right mean and the wrong shape
  moves it too.

- **PIT uniformity**: a chi-square test over `bins` equal bins, with at
  least 5 outcomes expected in each (limn's choice of minimum).

- **A binary outcome**: the calibration intercept (are the probabilities
  right on average?) and slope (are they too extreme, or too moderate?),
  with Wald 95% intervals, and a reliability table of equal-count bins
  with Wilson intervals. Interval coverage adds nothing for a 0/1
  outcome.

A coverage or bias departure is called when it is more than 2.58
standard errors from its ideal (the two-sided 99% point of the normal
distribution, limn's own choice), and a binary measure when its 95%
interval excludes the ideal. The intercept's ideal is 0 and the slope's
is 1 (Van Calster et al. 2019, p. 4).

A biased mean is reported first. It alone can make the intervals cover
too little, since they are centred in the wrong place. So a shortfall in
coverage beside it is reported as one the bias can cause, and the
remedies that widen intervals wait until the mean is fixed. A
calibration slope whose interval lies wholly below 0 is reported as
reversed, a verdict of limn's own: the predictions order the outcomes
the wrong way, which is more than being too extreme.

## Verdicts

What each check means and what to do about it, in the words it is
printed with. A reference supports what was found and why; a remedy
draws on it only where it says so, and anything marked as limn's own
reasoning is not the source's:

- **bias** (warn): Biased: outcomes are on average some amount above or
  below their predicted means (SE s). The mean structure is missing
  something the outcomes follow. illume::ilm_remedies() and illume's
  residual checks point at which term.

- **coverage** (warn): Too narrow: a share of outcomes fell inside the
  central L% of their predictive distributions, against L% expected. A
  biased mean alone can cause this: the intervals are centred in the
  wrong place. Fix the mean first, as above, then check the coverage
  again; widen the intervals only if it still falls short.

- **coverage** (warn): Too narrow: a share of outcomes fell inside the
  central L% of their predictive distributions, against L% expected. The
  predictions are more certain than the data bear out. Reference:
  Gneiting, Balabdaoui and Raftery (2007), Journal of the Royal
  Statistical Society B 69, 243-268, pp. 252-253. Overdispersion is the
  usual cause for counts (limn's own reasoning): refit with family =
  "nbinom" in illume::ilm_model(), which lets the variance exceed the
  mean.

- **coverage** (warn): Too narrow: a share of outcomes fell inside the
  central L% of their predictive distributions, against L% expected. The
  predictions are more certain than the data bear out. Reference:
  Gneiting, Balabdaoui and Raftery (2007), Journal of the Royal
  Statistical Society B 69, 243-268, pp. 252-253. Let the dispersion
  vary with the predictors through a dispformula in illume::ilm_model();
  with few groups and predictions for new units, refit with reml = TRUE.

- **coverage** (warn): Too narrow: a share of outcomes fell inside the
  central L% of their predictive distributions, against L% expected. The
  predictions are more certain than the data bear out. Reference:
  Gneiting, Balabdaoui and Raftery (2007), Journal of the Royal
  Statistical Society B 69, 243-268, pp. 252-253. illume::ilm_remedies()
  lists the model changes the fit's own checks support, and
  illume::ilm_apply_remedy() refits the model with the one you choose.

- **coverage** (warn): Too wide: a share of outcomes fell inside the
  central L% of their predictive distributions, against L% expected. The
  predictions are less certain than the data show. Reference: Gneiting,
  Balabdaoui and Raftery (2007), Journal of the Royal Statistical
  Society B 69, 243-268, pp. 252-253. An overdispersed family, or a
  dispersion model the data do not need, widens every interval (limn's
  own reasoning): illume::ilm_remedies() lists the model changes the
  fit's own checks support, and illume::ilm_apply_remedy() refits the
  model with the one you choose.

- **calibration** (ok): Consistent with calibration: coverage within
  sampling error at every level, no bias in the mean, and the PITs
  consistent with uniform (p).

- **calibration intercept** (warn): Off on average: a share of outcomes
  were 1, against another predicted; the calibration intercept is b (95%
  interval l to u). The probabilities are too high or too low overall.
  Reference: Van Calster and others (2019), BMC Medicine 17, 230, pp. 2
  and 4. A difference between the data the model was fitted to and the
  data scored, such as a different rate of the outcome, can cause this
  (the reference, p. 2): illume::ilm_remedies() lists the model changes
  the fit's own checks support, and illume::ilm_apply_remedy() refits
  the model with the one you choose.

- **calibration slope** (warn): Reversed: the calibration slope is b
  (95% interval l to u); the predictions order the outcomes the wrong
  way. This verdict is limn's own: the reference defines the calibration
  slope, with a target of 1, and reads a slope below 1 as too extreme,
  but does not discuss one below 0. Outcomes the model rates more likely
  turn out less likely; it may have fitted noise, or the scored data may
  differ from the data it was fitted to (limn's own reasoning).
  Reference: Van Calster and others (2019), BMC Medicine 17, 230, p. 4.
  Use fewer terms, or penalised ones: smooths in illume::ilm_model() are
  penalised; and check that the scored data come from the same
  population and coding.

- **calibration slope** (warn): Overconfident: the calibration slope is
  b (95% interval l to u); predicted probabilities are too extreme. The
  model has fitted noise, so its probabilities stray too far from the
  average. Reference: Van Calster and others (2019), BMC Medicine 17,
  230, pp. 2 and 4. Use fewer terms, or penalised ones (the
  reference, p. 4): smooths in illume::ilm_model() are penalised.

- **calibration slope** (warn): Underconfident: the calibration slope is
  b (95% interval l to u); predicted probabilities are too moderate. The
  probabilities stay too close to the average. Reference: Van Calster
  and others (2019), BMC Medicine 17, 230, p. 4. A penalty stronger than
  the data warrant can cause this (limn's own reasoning):
  illume::ilm_remedies() lists the model changes the fit's own checks
  support, and illume::ilm_apply_remedy() refits the model with the one
  you choose.

- **calibration** (ok): Consistent with calibration: the calibration
  intercept includes 0 and the slope includes 1.

## The rows are treated as independent

Rows that share a group, a series or a forecast origin are not
independent, and then the standard errors and tests here are too
confident. For forecasts, check calibration horizon by horizon, or over
rows from different series and origins.

## References

Gneiting, T., Balabdaoui, F. and Raftery, A. E. (2007). Probabilistic
forecasts, calibration and sharpness. *Journal of the Royal Statistical
Society B* 69, 243-268. Section 3.1, pp. 251-253: the PIT, its
uniformity, coverage read from it, and the histogram's shapes.

Van Calster, B., McLernon, D. J., van Smeden, M., Wynants, L. and
Steyerberg, E. W. (2019). Calibration: the Achilles heel of predictive
analytics. *BMC Medicine* 17, 230. p. 2: causes of poor calibration; p.
4: weak calibration, the intercept and slope and their targets, and
remedies.

## See also

[`ilm_pit()`](https://huttoncp.github.io/limn/reference/ilm_pit.md),
[`ilm_plot_calibration()`](https://huttoncp.github.io/limn/reference/ilm_plot_calibration.md),
[`ilm_scores.ilm_pred()`](https://huttoncp.github.io/limn/reference/ilm_scores.ilm_pred.md)

## Examples

``` r
set.seed(1)
d <- data.frame(x = rnorm(400))
d$y <- rnbinom(400, size = 1.5, mu = exp(1 + 0.5 * d$x))
train <- d[1:200, ]; test <- d[201:400, ]
## Poisson on overdispersed counts: intervals too narrow, and the remedy
fp <- illume::ilm_model(y ~ x, data = train, family = "poisson", verbose = FALSE)
illume::ilm_calibration(ilm_predict(fp, test, target = "response", nsim = 200), test$y)
#> <calibration> 200 outcomes, poisson predictive distributions
#>   PIT uniformity: chi-square 53.7 on 9 df, p < 0.001
#>   coverage by level (from the PITs):
#>      50%:  35.5% (SE 3.54 percentage points), mean width 2.15
#>      80%:  61.5% (SE 2.83 percentage points), mean width 4.10
#>      90%:  75.5% (SE 2.12 percentage points), mean width 5.17
#>      95%:  81.0% (SE 1.54 percentage points), mean width 6.05
#>   Too narrow: 81.0% of outcomes fell inside the central 95% of their predictive distributions, against 95% expected. The predictions are more certain than the data bear out. Reference: Gneiting, Balabdaoui and Raftery (2007), Journal of the Royal Statistical Society B 69, 243-268, pp. 252-253.
#>   Remedy: Overdispersion is the usual cause for counts (limn's own reasoning): refit with family = "nbinom" in illume::ilm_model(), which lets the variance exceed the mean.
## the negative binomial the verdict names
fn <- illume::ilm_model(y ~ x, data = train, family = "nbinom", verbose = FALSE)
illume::ilm_calibration(ilm_predict(fn, test, target = "response", nsim = 200), test$y)
#> <calibration> 200 outcomes, nbinom predictive distributions
#>   PIT uniformity: chi-square 14.6 on 9 df, p = 0.103
#>   coverage by level (from the PITs):
#>      50%:  54.0% (SE 3.54 percentage points), mean width 3.31
#>      80%:  82.5% (SE 2.83 percentage points), mean width 6.42
#>      90%:  91.5% (SE 2.12 percentage points), mean width 8.25
#>      95%:  97.0% (SE 1.54 percentage points), mean width 10.1
#>   Consistent with calibration: coverage within sampling error at every level, no bias in the mean, and the PITs consistent with uniform (p = 0.103).
```
