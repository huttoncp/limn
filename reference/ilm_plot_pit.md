# Histogram of probability integral transforms

Plots PITs from
[`ilm_pit()`](https://huttoncp.github.io/limn/reference/ilm_pit.md) as a
histogram on the density scale, where a correct predictive distribution
gives a flat histogram at height 1. The shaded band is where each bar
falls with probability `level` if the PITs really are uniform.

## Usage

``` r
ilm_plot_pit(x, bins = 10L, level = 0.95, ...)
```

## Arguments

- x:

  PITs from
  [`ilm_pit()`](https://huttoncp.github.io/limn/reference/ilm_pit.md),
  or any numbers between 0 and 1.

- bins:

  Number of bars.

- level:

  Coverage of the shaded band under uniformity.

- ...:

  Passed to
  [`tinyplot::tinyplot()`](https://grantmcdermott.com/tinyplot/man/tinyplot.html).

## Value

The data behind the plot, invisibly: each bar's limits, count and
height, and the band.

## Reading the plot

- **Flat, within the band:** the outcomes fell where their predictive
  distributions said they would. The band holds for each bar on its own,
  so with many bars one may stray outside it by chance.

- **U shape**, high at both ends: too many outcomes fell in the tails.
  This often means predictive distributions too narrow.

- **Hump**, high in the middle: intervals too wide on average.

- **Triangle**, sloping one way: the predictions are biased, too high if
  the bars rise towards 0, too low if they rise towards 1.

The shapes are those Gneiting, Balabdaoui and Raftery (2007, p. 253)
describe; which way a triangle points is limn's own reading of it. For a
verdict with what to do about it, use
[ilm_calibration()](https://huttoncp.github.io/limn/reference/ilm_calibration.ilm_pred.md).

## References

Gneiting, T., Balabdaoui, F. and Raftery, A. E. (2007). Probabilistic
forecasts, calibration and sharpness. *Journal of the Royal Statistical
Society B* 69, 243-268. Section 3.1, p. 253.

## See also

[`ilm_pit()`](https://huttoncp.github.io/limn/reference/ilm_pit.md),
which computes the PITs;
[ilm_calibration()](https://huttoncp.github.io/limn/reference/ilm_calibration.ilm_pred.md)
for a verdict.

## Examples

``` r
u <- runif(500)
ilm_plot_pit(u)
```
