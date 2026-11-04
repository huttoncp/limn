# Covariate paths from another tool, as draws

Turns simulated paths of covariates made elsewhere into `newdata`'s
blocks of draws, marked by a `.draw` column, which
[`ilm_forecast()`](https://huttoncp.github.io/limn/reference/ilm_forecast.md)
and
[`ilm_predict()`](https://huttoncp.github.io/limn/reference/ilm_predict.md)
mix over (see *Covariates given as draws* in
[`ilm_predict()`](https://huttoncp.github.io/limn/reference/ilm_predict.md)).
The paths might come from fable's `generate()`, a weather ensemble, or a
scenario generator. For a covariate modelled with illume,
[`ilm_covariate_paths()`](https://huttoncp.github.io/limn/reference/ilm_covariate_paths.md)
draws them directly.

## Usage

``` r
ilm_as_draws(paths, newdata, by, path = ".rep", values = NULL)
```

## Arguments

- paths:

  A data frame in long form: one row per path and row of `newdata`, with
  a path identifier, the `by` columns and the simulated values.

- newdata:

  The rows to forecast.

- by:

  The columns, in both, that identify a row of `newdata`, such as its
  series and time.

- path:

  The column of `paths` that identifies a path. The default, `".rep"`,
  is the name fable's `generate()` uses (fpp3, section 5.5).

- values:

  The drawn columns, as a named character vector mapping each column of
  `newdata` to its column in `paths`: for fable's `generate()`, whose
  simulated values are in `.sim`, `c(x = ".sim")`. By default every
  column of `paths` that is not in `by` or `path` and has the same name
  in `newdata`.

## Value

`newdata` repeated once per path, with a `.draw` column and the drawn
columns filled in from the paths.

## Details

Every path must give a value for every row of `newdata`, matched on the
`by` columns, and paths are numbered in the order they first appear.
Columns of `newdata` that are not drawn are copied into every block.

## References

Hyndman, R. J. and Athanasopoulos, G. (2021). *Forecasting: Principles
and Practice*, 3rd edition. OTexts. Section 5.5, "Distributional
forecasts and prediction intervals": fable's `generate()` and its `.rep`
and `.sim` columns (online edition,
<https://otexts.com/fpp3/toolbox.html#sec-prediction-intervals>, read 27
October 2026).

## See also

[`ilm_covariate_paths()`](https://huttoncp.github.io/limn/reference/ilm_covariate_paths.md)
for paths from an illume model;
[`ilm_forecast()`](https://huttoncp.github.io/limn/reference/ilm_forecast.md),
which mixes over the blocks.

## Examples

``` r
fut <- data.frame(id = rep(c("a", "b"), each = 2), t = rep(31:32, 2), x = NA)
sims <- expand.grid(t = 31:32, id = c("a", "b"), .rep = 1:3)
sims$.sim <- rnorm(nrow(sims))
ilm_as_draws(sims, fut, by = c("id", "t"), values = c(x = ".sim"))
#>    .draw id  t          x
#> 1      1  a 31 -1.3457937
#> 2      1  a 32  1.0336654
#> 3      1  b 31 -0.8117765
#> 4      1  b 32  1.8017255
#> 5      2  a 31  1.7715420
#> 6      2  a 32 -1.4546914
#> 7      2  b 31 -0.8456543
#> 8      2  b 32 -1.2504797
#> 9      3  a 31  0.6672881
#> 10     3  a 32 -1.2907697
#> 11     3  b 31 -2.0350035
#> 12     3  b 32  2.0213470
```
