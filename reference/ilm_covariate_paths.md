# Draws of a covariate's future, from its own model

When a covariate is not known in advance, its uncertainty belongs in the
forecast. `ilm_covariate_paths()` simulates future paths of the
covariate from an illume model of it, and returns `newdata` stacked into
one block of rows per path, marked by a `.draw` column.
[`ilm_forecast()`](https://huttoncp.github.io/limn/reference/ilm_forecast.md)
and
[`ilm_predict()`](https://huttoncp.github.io/limn/reference/ilm_predict.md)
then mix over the paths; see *Covariates given as draws* in
[`ilm_predict()`](https://huttoncp.github.io/limn/reference/ilm_predict.md).

## Usage

``` r
ilm_covariate_paths(
  fit,
  newdata,
  nsim = 1000L,
  uncertainty = c("joint", "fixed", "latent", "none"),
  seed = 1L
)
```

## Arguments

- fit:

  A fitted `"ilm_model"` of the covariate. Its response names the column
  of `newdata` to draw. Or a list of such fits, one per covariate.

- newdata:

  The rows to forecast, as for
  [`ilm_forecast()`](https://huttoncp.github.io/limn/reference/ilm_forecast.md),
  with the columns the covariate models need. The drawn columns may be
  missing, as
  [`ilm_future()`](https://huttoncp.github.io/limn/reference/ilm_future.md)
  leaves them.

- nsim:

  The number of paths.

- uncertainty:

  How the covariate model's own parameters enter: `"joint"` (the
  default), `"fixed"`, `"latent"` or `"none"`, as in
  [`ilm_predict()`](https://huttoncp.github.io/limn/reference/ilm_predict.md).

- seed:

  Seed for the draws. The caller's random-number state is restored
  afterwards.

## Value

`newdata` repeated once per path, with a `.draw` column (1 to `nsim`)
and the drawn columns filled in.

## How the paths are drawn

Each path is one draw from the covariate model's predictive distribution
of a new observation. That draw has three parts:

- a draw of the model's parameters, under `uncertainty`;

- each series' latent path continued over the new times, from its fitted
  state (`groups = "fitted"`);

- a new observation at each time.

So a path is joint over its times. Together, the paths carry both the
covariate's future shocks and the uncertainty in its model's parameters.
Give `nsim` the same value as the forecast's own `nsim`: parameter draw
s of the forecast is then paired with path s.

## Several covariates

Give a list of fits, one per covariate. Each covariate is drawn from its
own model, independently of the others. Covariates that move together
would need a joint model of them, which is not in limn yet.

## See also

[`ilm_as_draws()`](https://huttoncp.github.io/limn/reference/ilm_as_draws.md)
for paths simulated by another tool;
[`ilm_future()`](https://huttoncp.github.io/limn/reference/ilm_future.md)
for the rows.

## Examples

``` r
set.seed(1)
d <- data.frame(id = factor(rep(c("a", "b"), each = 30)), t = rep(1:30, 2))
d$x <- ave(rnorm(60), d$id, FUN = function(e) stats::filter(e, 0.8, "recursive"))
d$y <- 2 + 0.5 * d$x + rnorm(60, 0, 0.3)
fx <- illume::ilm_model(x ~ 1, data = d, family = "gaussian",
                        ar = illume::ilm_ar1(~ t | id), verbose = FALSE)
#> ilm_model(): the residual SD has run to zero: the model's other terms leave no residual variation, usually because a correlation over time or a random effect has about one observation per cell or level and has taken up the noise. The fixed effects and their standard errors are usable; sigma's own estimate is where the optimiser stopped, and is held there.
fy <- illume::ilm_model(y ~ x, data = d, family = "gaussian", verbose = FALSE)
fut <- ilm_future(d, time = "t", group = "id", h = 3, response = "y")
#> These columns vary over time, so their future values are not known: x. They are left missing: fill them in with the values the forecast should assume.
paths <- ilm_covariate_paths(fx, fut, nsim = 200)
head(paths)
#>   .draw id  t           x  y horizon
#> 1     1  a 31 -0.09111319 NA       1
#> 2     1  a 32  0.22876583 NA       2
#> 3     1  a 33 -0.34277032 NA       3
#> 4     1  b 31  1.48188772 NA       1
#> 5     1  b 32  0.47211154 NA       2
#> 6     1  b 33  0.42517592 NA       3
fc <- ilm_forecast(fy, paths, nsim = 200)
as.data.frame(fc)
#>   id  t  x  y horizon new_level estimate     mean lower_0.9 upper_0.9
#> 1  a 31 NA NA       1     FALSE 1.901752 1.899944  1.082757  2.698562
#> 2  a 32 NA NA       2     FALSE 2.022371 2.018640  1.112669  2.921620
#> 3  a 33 NA NA       3     FALSE 2.119380 2.101238  1.087254  3.050851
#> 4  b 31 NA NA       1     FALSE 2.338161 2.353596  1.539417  3.232534
#> 5  b 32 NA NA       2     FALSE 2.276033 2.281953  1.433222  3.161757
#> 6  b 33 NA NA       3     FALSE 2.292942 2.298043  1.364966  3.252102
```
