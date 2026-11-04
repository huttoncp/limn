# Package index

## Predictions and forecasts

Predictive distributions from a fitted illume model, the rows a forecast
needs, and covariates whose future is not known.

- [`ilm_predict()`](https://huttoncp.github.io/limn/reference/ilm_predict.md)
  : Predictive distributions from a fitted illume model
- [`ilm_forecast()`](https://huttoncp.github.io/limn/reference/ilm_forecast.md)
  : Forecasts from a fitted illume model
- [`ilm_future()`](https://huttoncp.github.io/limn/reference/ilm_future.md)
  : The rows a forecast will need: the next time points of every series
- [`ilm_covariate_paths()`](https://huttoncp.github.io/limn/reference/ilm_covariate_paths.md)
  : Draws of a covariate's future, from its own model
- [`ilm_as_draws()`](https://huttoncp.github.io/limn/reference/ilm_as_draws.md)
  : Covariate paths from another tool, as draws
- [`as.data.frame(`*`<ilm_pred>`*`)`](https://huttoncp.github.io/limn/reference/as.data.frame.ilm_pred.md)
  [`summary(`*`<ilm_pred>`*`)`](https://huttoncp.github.io/limn/reference/as.data.frame.ilm_pred.md)
  : Point estimates and intervals from a predictive distribution

## Totals and averages

Sums and weighted averages of predictions over rows, taken draw by draw.

- [`ilm_aggregate()`](https://huttoncp.github.io/limn/reference/ilm_aggregate.md)
  : Totals and averages of predictions over rows
- [`ilm_average()`](https://huttoncp.github.io/limn/reference/ilm_average.md)
  : Weighted averages of predicted means: standardisation and
  post-stratification

## Scores and calibration

How well predictive distributions match the outcomes: proper scores,
scaled errors, probability integral transforms and calibration.

- [`ilm_scores(`*`<ilm_pred>`*`)`](https://huttoncp.github.io/limn/reference/ilm_scores.ilm_pred.md)
  : Proper scores of a prediction against observed outcomes
- [`ilm_naive_scale()`](https://huttoncp.github.io/limn/reference/ilm_naive_scale.md)
  : Each series' in-sample naive-forecast error, the scale of MASE and
  RMSSE
- [`ilm_pit()`](https://huttoncp.github.io/limn/reference/ilm_pit.md) :
  Probability integral transforms of observed outcomes
- [`ilm_calibration(`*`<ilm_pred>`*`)`](https://huttoncp.github.io/limn/reference/ilm_calibration.ilm_pred.md)
  : Calibration of a prediction against observed outcomes

## Cross-validation and comparison

Folds that follow the data’s structure, cross-validated scores, and
model comparisons with standard errors.

- [`ilm_folds()`](https://huttoncp.github.io/limn/reference/ilm_folds.md)
  : Cross-validation folds that follow the data's structure
- [`ilm_cv()`](https://huttoncp.github.io/limn/reference/ilm_cv.md) :
  Cross-validated scores of a fitted illume model
- [`ilm_compare()`](https://huttoncp.github.io/limn/reference/ilm_compare.md)
  : Compare models' predictive scores, with standard errors that respect
  the data

## Plots

Fan charts, scores by horizon, model comparisons, PIT histograms and
calibration.

- [`ilm_plot_fan()`](https://huttoncp.github.io/limn/reference/ilm_plot_fan.md)
  : Fan chart of predictive distributions over time
- [`ilm_plot_horizon()`](https://huttoncp.github.io/limn/reference/ilm_plot_horizon.md)
  : Scores by forecast horizon
- [`ilm_plot_compare()`](https://huttoncp.github.io/limn/reference/ilm_plot_compare.md)
  : Plot model comparisons with their intervals
- [`ilm_plot_pit()`](https://huttoncp.github.io/limn/reference/ilm_plot_pit.md)
  : Histogram of probability integral transforms
- [`ilm_plot_calibration()`](https://huttoncp.github.io/limn/reference/ilm_plot_calibration.md)
  : Plot the calibration of predictive distributions

## Extending limn

Methods a package defines so that limn can predict from, and refit, its
own models.

- [`ilm_cv_refit()`](https://huttoncp.github.io/limn/reference/ilm_cv_refit.md)
  : Refit a model to a fold's training rows
- [`ilm_predict_draws()`](https://huttoncp.github.io/limn/reference/ilm_predict_draws.md)
  : The parameter draws behind a prediction
