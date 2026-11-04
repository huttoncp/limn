# The parameter draws behind a prediction

`ilm_predict_draws()` is the step of
[`ilm_predict()`](https://huttoncp.github.io/limn/reference/ilm_predict.md)
that draws the model's parameters. Every prediction, forecast, score,
calibration and cross-validation in limn is computed from these draws.
You will not normally call it yourself: it is a generic, so that a
package fitting illume's models another way can supply the draws from
its own fits, and limn then works on those fits unchanged.

## Usage

``` r
ilm_predict_draws(fit, nsim, seed, uncertainty, natural, ...)

# S3 method for class 'ilm_model'
ilm_predict_draws(fit, nsim, seed, uncertainty, natural, ...)
```

## Arguments

- fit:

  A fitted model: an `"ilm_model"`, or a class that inherits from it.

- nsim:

  The number of draws to return.

- seed:

  The seed the draws are made from. The same seed must give the same
  draws, whatever else is running.

- uncertainty:

  Which parameters are drawn, as in
  [`ilm_predict()`](https://huttoncp.github.io/limn/reference/ilm_predict.md):
  `"joint"` (all of them), `"fixed"` (the variances held at their
  estimates), `"latent"` (every parameter held, so only the random
  effects and latent cells are drawn) or `"none"` (one draw, at the
  estimates). A method may refuse a value it has no meaning for, with an
  error that says which values it supports.

- natural:

  Whether each draw's variances are also needed on their natural scale.
  limn asks for them only when a unit's effects are drawn afresh or
  averaged over.

- ...:

  Passed on to methods.

## Value

An object of class `"ilm_draws"`, as
[`illume::ilm_draws()`](https://rdrr.io/pkg/illume/man/ilm_draws.html)
returns:

- `draws`: a matrix with one row per parameter, in illume's order and on
  its internal scale, and `nsim` columns (one for
  `uncertainty = "none"`);

- `mode`: the estimates;

- `map`: what each row of `draws` is, as
  [`illume::ilm_draws()`](https://rdrr.io/pkg/illume/man/ilm_draws.html)
  gives it;

- `natural`: when `natural = TRUE`, each draw's variances on their
  natural scale, which at the estimates equal
  [`illume::ilm_varcorr()`](https://rdrr.io/pkg/illume/man/ilm_varcorr.html)'s.

[`ilm_predict()`](https://huttoncp.github.io/limn/reference/ilm_predict.md)
checks this shape, and stops with a message saying what is wrong when a
method returns something else.

## Details

The method for illume's fits draws from the fit's normal approximation,
with
[`illume::ilm_draws()`](https://rdrr.io/pkg/illume/man/ilm_draws.html).
A method for another class of fit must return the same shape: the class
of fit must inherit from `"ilm_model"`, so that illume's
[`illume::ilm_dist()`](https://rdrr.io/pkg/illume/man/ilm_dist.html),
[`illume::ilm_varcorr()`](https://rdrr.io/pkg/illume/man/ilm_varcorr.html),
[`illume::ilm_matrices()`](https://rdrr.io/pkg/illume/man/ilm_matrices.html)
and
[`illume::ilm_cells()`](https://rdrr.io/pkg/illume/man/ilm_cells.html)
accept it, and the draws must be laid out as
[`illume::ilm_draws()`](https://rdrr.io/pkg/illume/man/ilm_draws.html)
lays them out.

## Checks

What each check means and what to do about it, in the words it is
printed with. A reference supports what was found and why; a remedy
draws on it only where it says so, and anything marked as limn's own
reasoning is not the source's:

- **parameter draws** (fail): The parameter draws for this fit are not
  in the shape limn needs: what is wrong. Every prediction is computed
  from the draws ilm_predict_draws() returns, laid out as
  illume::ilm_draws() lays them out. The method for this class of fit
  returned something else. If the fit's class comes from another
  package, its ilm_predict_draws() method needs fixing there:
  ?ilm_predict_draws gives the shape. For an illume fit this is a fault
  in limn; please report it.

## See also

[`ilm_cv_refit()`](https://huttoncp.github.io/limn/reference/ilm_cv_refit.md),
the other extension point.

## Examples

``` r
set.seed(1)
d <- data.frame(g = factor(rep(letters[1:8], each = 5)), x = rnorm(40))
d$y <- rpois(40, exp(0.5 + 0.3 * d$x))
fit <- illume::ilm_model(y ~ x + (1 | g), data = d, family = "poisson",
                         verbose = FALSE)
#> ilm_model(): the random-effect covariance of `g` sits at the edge of its range -- a variance of zero or a correlation of +/-1 -- where the data cannot resolve it. The fixed effects and their standard errors are still usable; summary() says what else is. If the term belongs in the model, boundary = "avoid" keeps it inside its range with a small penalty: it is then assumed nonzero rather than estimated at zero, so do not test whether it is; its variance comes out larger, and for a binary or categorical outcome the fixed effects a little further from zero -- markedly so when a category is rare.
dr <- ilm_predict_draws(fit, nsim = 10, seed = 1, uncertainty = "joint",
                        natural = FALSE)
dim(dr$draws)
#> [1] 11 10
```
