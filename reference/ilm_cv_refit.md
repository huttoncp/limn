# Refit a model to a fold's training rows

`ilm_cv_refit()` is the step of
[`ilm_cv()`](https://huttoncp.github.io/limn/reference/ilm_cv.md) that
refits the model to each fold's training rows. You will not normally
call it yourself: it is a generic, so that a package fitting illume's
models another way can refit its own fits, and
[`ilm_cv()`](https://huttoncp.github.io/limn/reference/ilm_cv.md) then
cross-validates them unchanged.

## Usage

``` r
ilm_cv_refit(fit, data, ...)

# S3 method for class 'ilm_model'
ilm_cv_refit(fit, data, ...)
```

## Arguments

- fit:

  A fitted model: an `"ilm_model"`, or a class that inherits from it.

- data:

  The fold's training rows.

- ...:

  Passed on to methods.

## Value

The model refitted to `data`.

## Details

The method for illume's fits is
[`illume::ilm_refit()`](https://rdrr.io/pkg/illume/man/ilm_refit.html),
which repeats the model's own call on the new rows. A method for another
class of fit must return a fit that
[`ilm_predict()`](https://huttoncp.github.io/limn/reference/ilm_predict.md)
accepts, and must give the same fit whether it runs in the main R
session or in one of the worker processes that `ilm_cv(ncores = )`
starts. Each worker loads limn and the package that registered the
method, and runs one thread.

## See also

[`ilm_predict_draws()`](https://huttoncp.github.io/limn/reference/ilm_predict_draws.md),
the other extension point.

## Examples

``` r
set.seed(1)
d <- data.frame(g = factor(rep(letters[1:8], each = 5)), x = rnorm(40))
d$y <- rpois(40, exp(0.5 + 0.3 * d$x))
fit <- illume::ilm_model(y ~ x + (1 | g), data = d, family = "poisson",
                         verbose = FALSE)
#> ilm_model(): the random-effect covariance of `g` sits at the edge of its range -- a variance of zero or a correlation of +/-1 -- where the data cannot resolve it. The fixed effects and their standard errors are still usable; summary() says what else is. If the term belongs in the model, boundary = "avoid" keeps it inside its range with a small penalty: it is then assumed nonzero rather than estimated at zero, so do not test whether it is; its variance comes out larger, and for a binary or categorical outcome the fixed effects a little further from zero -- markedly so when a category is rare.
f1 <- ilm_cv_refit(fit, d[d$g != "a", ])
#> ilm_model(): the random-effect covariance of `g` sits at the edge of its range -- a variance of zero or a correlation of +/-1 -- where the data cannot resolve it. The fixed effects and their standard errors are still usable; summary() says what else is. If the term belongs in the model, boundary = "avoid" keeps it inside its range with a small penalty: it is then assumed nonzero rather than estimated at zero, so do not test whether it is; its variance comes out larger, and for a binary or categorical outcome the fixed effects a little further from zero -- markedly so when a category is rare.
```
