# limn

**A Projection Engine for Distributional Forecasting.**

Prediction and forecasting for models fitted by
[`illume`](https://github.com/huttoncp/illume): predictive
distributions, forecasts, proper scoring rules, calibration,
cross-validation and model comparison.

`limn` is the third package in a family, after
[`illumex`](https://github.com/huttoncp/illumex) for exploration and
`illume` for inference. All three use the `ilm_` prefix, so one prefix
finds a function in any of them.
[`library(limn)`](https://github.com/huttoncp/limn) attaches `illume`,
which attaches `illumex`.

## Status: early development

`limn` is pre-release: this is a development version, published before
its first CRAN release. Its core is built and validated, and its
interfaces may still change before then. The validation studies and
their findings, with the prototypes behind some of its defaults, are in
[`studies`](https://huttoncp.github.io/limn/studies).

What exists now:

- Prediction and forecasting:
  - [`ilm_predict()`](https://huttoncp.github.io/limn/reference/ilm_predict.md):
    the prediction engine. It gives the predictive distribution of the
    mean, a new observation or the linear predictor for new rows, with
    the fit’s own uncertainty. How the random effects enter is always
    named.
  - [`ilm_forecast()`](https://huttoncp.github.io/limn/reference/ilm_forecast.md):
    forecasts past a series’ last observed time. A correlation over time
    (`ilm_ar1()`, `ilm_car1()`, `ilm_rw1()`) is carried forward from the
    fit’s own parameters, with no refit, as one path per draw.
  - [`ilm_future()`](https://huttoncp.github.io/limn/reference/ilm_future.md):
    the rows a forecast needs, the next time points of every series,
    with time-varying covariates left for you to fill in.
  - [`ilm_covariate_paths()`](https://huttoncp.github.io/limn/reference/ilm_covariate_paths.md)
    and
    [`ilm_as_draws()`](https://huttoncp.github.io/limn/reference/ilm_as_draws.md):
    covariates whose future is uncertain, given as draws, so the
    forecast mixes over them. The draws are simulated from an `illume`
    model of the covariate, or taken from elsewhere.
- Totals and averages:
  - [`ilm_aggregate()`](https://huttoncp.github.io/limn/reference/ilm_aggregate.md):
    totals or averages over a series’ horizon, over every series at each
    time, or over all rows. They are taken draw by draw, so a path’s
    correlation carries into them.
  - [`ilm_average()`](https://huttoncp.github.io/limn/reference/ilm_average.md):
    weighted averages of predicted means, for standardisation or
    post-stratification (as in MRP).
- Evaluation:
  - `ilm_scores()` on a prediction, as a method of `illume`’s generic.
    It gives the log score, the CRPS, and the Brier, pinball, interval,
    squared and absolute scores. MASE and RMSSE take their scale from
    the training data or from
    [`ilm_naive_scale()`](https://huttoncp.github.io/limn/reference/ilm_naive_scale.md).
  - [`ilm_pit()`](https://huttoncp.github.io/limn/reference/ilm_pit.md)
    and
    [`ilm_plot_pit()`](https://huttoncp.github.io/limn/reference/ilm_plot_pit.md):
    probability integral transforms, randomised for discrete outcomes,
    and their histogram.
  - [`ilm_folds()`](https://huttoncp.github.io/limn/reference/ilm_folds.md):
    cross-validation folds that follow the data’s structure, with
    random, grouped, blocked and rolling-origin folds.
  - [`ilm_cv()`](https://huttoncp.github.io/limn/reference/ilm_cv.md): a
    model refitted to each fold, then scored on the rows it wasn’t
    fitted to.
  - [`ilm_compare()`](https://huttoncp.github.io/limn/reference/ilm_compare.md):
    two models’ scores compared, with standard errors clustered by group
    or fold, or Diebold-Mariano over forecast origins.
  - `ilm_calibration()` on a prediction, as a method of `illume`’s
    generic: coverage, bias, and a binary outcome’s calibration
    intercept and slope, with a verdict that names a remedy.
- Diagnostic plots:
  - [`ilm_plot_fan()`](https://huttoncp.github.io/limn/reference/ilm_plot_fan.md):
    a fan chart of predictive distributions over time;
  - [`ilm_plot_horizon()`](https://huttoncp.github.io/limn/reference/ilm_plot_horizon.md):
    scores by forecast horizon;
  - [`ilm_plot_compare()`](https://huttoncp.github.io/limn/reference/ilm_plot_compare.md):
    comparisons with their intervals;
  - [`ilm_plot_calibration()`](https://huttoncp.github.io/limn/reference/ilm_plot_calibration.md):
    coverage, or reliability for a binary outcome.

## Principles

- **A prediction is a distribution.** It is a mixture over draws of the
  fitted parameters, so its intervals carry the fit’s own uncertainty.
  Every score is computed from that mixture, not from responses
  simulated from it.
- **How the random effects enter is always named.** The choices are a
  unit seen in the fit, a new unit, a typical unit, or the population
  average. These answer different questions, and there is no default
  that suits every one.
- **Evaluation follows the data’s structure.** Grouped folds score
  predictions for new units, and rolling origins score forecasts.
  Comparisons carry standard errors clustered the same way.
- **Every validity claim is measured** against an independent
  implementation where one exists, and by simulation otherwise. The
  scripts are in `studies/`.

## Installation

`limn` and the two packages it builds on are installed from GitHub:

``` r

# install.packages("remotes")
remotes::install_github("huttoncp/illumex")
remotes::install_github("huttoncp/illume")
remotes::install_github("huttoncp/limn")
```

## How this package was built

limn’s code was written at Craig P. Hutton’s direction by his Claude
Code sessions. Craig sets the design and approves every change that
reaches the main branch.

**What the model did.** It prototyped and measured the design, and
drafted the implementation, the documentation and the validation
studies. It proposed designs, and argued for them when it thought a
choice was wrong.

**What the author did.** Specified what the package should be and what
it should refuse to do. Set the standing constraints:

- `illume`’s exported interface only;
- no benchmark methods inside the package;
- no claim without a measurement behind it.

Made the decisions the design turned on: the shared prefix, that the
user names how random effects enter a prediction, and which parts are
built in `limn` and which in `illume`. Directs the order of the work.

**Status of independent review.** The author reviews the package and
tests it independently of the development process. That work is in
progress, and this paragraph will say when it is complete. Until then
the package is pre-release, and should be treated as such.

The disclosure matters less than it might seem. The question about a
statistical package is whether its claims are checkable and whether
anyone checked them. `limn`’s answer is the same as `illume`’s: each
validity claim is a comparison against an independent implementation
that somebody else wrote, or a simulation study. The scripts and their
findings are in `studies/`.

## Acknowledgements

`limn`’s results are checked against three packages written by others:
`sandwich` for clustered standard errors (Zeileis, Köll and Graham 2020)
and a Newey-West long-run variance (Zeileis 2004), in the tests;
`scoringRules` for the CRPS and the log score (Jordan, Krüger and Lerch
2019), in the tests and studies; and `forecast` for the Diebold-Mariano
test (Hyndman and Khandakar 2008), in the studies. They are comparators:
`limn` shares their methods, not their code.

- Hyndman, R. J. and Khandakar, Y. (2008). Automatic time series
  forecasting: the forecast package for R. *Journal of Statistical
  Software* 27(3), 1-22. <https://doi.org/10.18637/jss.v027.i03>
- Jordan, A., Krüger, F. and Lerch, S. (2019). Evaluating probabilistic
  forecasts with scoringRules. *Journal of Statistical Software* 90(12),
  1-37. <https://doi.org/10.18637/jss.v090.i12>
- Zeileis, A. (2004). Econometric computing with HC and HAC covariance
  matrix estimators. *Journal of Statistical Software* 11(10), 1-17.
  <https://doi.org/10.18637/jss.v011.i10>
- Zeileis, A., Köll, S. and Graham, N. (2020). Various versatile
  variances: an object-oriented implementation of clustered covariances
  in R. *Journal of Statistical Software* 95(1), 1-36.
  <https://doi.org/10.18637/jss.v095.i01>

## License

GPL (\>= 2): the GNU General Public License, version 2 or (at your
option) any later version. The text of version 2 is in `LICENSE.md`.
