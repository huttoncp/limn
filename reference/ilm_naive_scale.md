# Each series' in-sample naive-forecast error, the scale of MASE and RMSSE

The mean absolute and mean squared error of the naive forecast (the
value `period` steps back) over each series' history: the scale that
makes MASE and RMSSE comparable across series. These are the
denominators of fpp3's scaled errors, the mean over t = `period` + 1 to
T of \|y_t - y\_(t - `period`)\| and of its square (section 5.8, "Scaled
errors"), which fpp3 credits to Hyndman and Koehler (2006). Compute it
from the data the model was trained on, and pass it to
[`ilm_scores.ilm_pred()`](https://huttoncp.github.io/limn/reference/ilm_scores.ilm_pred.md)
as `naive`.

## Usage

``` r
ilm_naive_scale(history, time, group = NULL, response, period = 1L)
```

## Arguments

- history:

  A data frame: the observed history, as the model was trained on it.

- time:

  The name of its time column.

- group:

  The name of the column identifying each series, or `NULL` for a single
  series. The prediction's rows need the same column.

- response:

  The name of the outcome column.

- period:

  The lag of the naive forecast: 1 for the naive forecast, the season's
  length for the seasonal naive.

## Value

A data frame with one row per series: `key`, and the naive error's mean
absolute value `abs` and mean square `sq`.

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

Hyndman, R. J. and Athanasopoulos, G. (2021). *Forecasting: Principles
and Practice*, 3rd edition. OTexts. Section 5.8, "Evaluating point
forecast accuracy" (online edition,
<https://otexts.com/fpp3/toolbox.html#sec-accuracy>, read 27 October
2026).

Hyndman, R. J. and Koehler, A. B. (2006). Another look at measures of
forecast accuracy. *International Journal of Forecasting* 22, 679-688.

## See also

[`ilm_scores.ilm_pred()`](https://huttoncp.github.io/limn/reference/ilm_scores.ilm_pred.md),
which uses it for `mase` and `rmsse`.

## Examples

``` r
d <- data.frame(id = rep(c("a", "b"), each = 5), t = rep(1:5, 2),
                y = c(1, 3, 2, 5, 4, 10, 12, 11, 15, 14))
ilm_naive_scale(d, time = "t", group = "id", response = "y")
#>   key  abs   sq
#> 1   a 1.75 3.75
#> 2   b 2.00 5.50
```
