# Distribution objects in the shape illume's ilm_dist() (export E4) will
# return, built from base R's d/p/q/r functions. They exist only here, to test
# limn's mixture machinery before E4 does. The package itself never builds a
# family's distribution: that is illume's, so the two cannot drift.
#
# Each takes `lp` (named list of link-scale linear predictors) and `par` (named
# list of the family's other parameters, here on a log scale).

test_dist_gaussian <- function() list(
  family = "gaussian", discrete = FALSE, lp_names = "mu", par_names = "logsd",
  d = function(y, lp, par, log = FALSE) stats::dnorm(y, lp$mu, exp(par$logsd), log = log),
  p = function(q, lp, par) stats::pnorm(q, lp$mu, exp(par$logsd)),
  q = function(prob, lp, par) stats::qnorm(prob, lp$mu, exp(par$logsd)),
  r = function(lp, par) stats::rnorm(length(lp$mu), lp$mu, exp(par$logsd)),
  mean = function(lp, par) lp$mu)

test_dist_poisson <- function() list(
  family = "poisson", discrete = TRUE, lp_names = "mu", par_names = character(0),
  d = function(y, lp, par, log = FALSE) stats::dpois(y, exp(lp$mu), log = log),
  p = function(q, lp, par) stats::ppois(q, exp(lp$mu)),
  q = function(prob, lp, par) stats::qpois(prob, exp(lp$mu)),
  r = function(lp, par) stats::rpois(length(lp$mu), exp(lp$mu)),
  mean = function(lp, par) exp(lp$mu))

test_dist_nbinom <- function() list(
  family = "nbinom", discrete = TRUE, lp_names = "mu", par_names = "logsize",
  d = function(y, lp, par, log = FALSE)
    stats::dnbinom(y, size = exp(par$logsize), mu = exp(lp$mu), log = log),
  p = function(q, lp, par) stats::pnbinom(q, size = exp(par$logsize), mu = exp(lp$mu)),
  q = function(prob, lp, par) stats::qnbinom(prob, size = exp(par$logsize), mu = exp(lp$mu)),
  r = function(lp, par) stats::rnbinom(length(lp$mu), size = exp(par$logsize), mu = exp(lp$mu)),
  mean = function(lp, par) exp(lp$mu))

test_dist_binomial <- function() list(
  family = "binomial", discrete = TRUE, lp_names = "mu", par_names = character(0),
  d = function(y, lp, par, log = FALSE) stats::dbinom(y, 1, stats::plogis(lp$mu), log = log),
  p = function(q, lp, par) stats::pbinom(q, 1, stats::plogis(lp$mu)),
  q = function(prob, lp, par) stats::qbinom(prob, 1, stats::plogis(lp$mu)),
  r = function(lp, par) stats::rbinom(length(lp$mu), 1, stats::plogis(lp$mu)),
  mean = function(lp, par) stats::plogis(lp$mu))

# A predictive distribution for n rows from S draws: row means `m` (length n),
# a per-draw shift with sd `spread` (parameter uncertainty), plus `par`.
test_pred <- function(dist, m, S = 50L, spread = 0.3, par = list(), seed = 1) {
  set.seed(seed)
  n <- length(m)
  shift <- stats::rnorm(S, 0, spread)
  lp <- list(mu = outer(m, shift, `+`))
  lmn_new_pred(dist, lp, par = par)
}
