## The continuous CRPS (studies/designs/crps_continuous.md): exact for gaussian
## mixtures, resolved panels on the natural axis for the other families, Inf
## where it is infinite, NA where integration fails, and an error for a family
## with no method.

## families in illume's parameterisations, location e and scale s
crps_fams <- list(
  gaussian = list(p = function(q, m, s, lt) stats::pnorm(q, m, s, lower.tail = lt), q = function(p, m, s) stats::qnorm(p, m, s)),
  lognormal = list(p = function(q, m, s, lt) stats::plnorm(q, m, s, lower.tail = lt), q = function(p, m, s) stats::qlnorm(p, m, s)),
  weibull = list(p = function(q, m, s, lt) stats::pweibull(q, 1 / s, exp(m), lower.tail = lt),
                 q = function(p, m, s) stats::qweibull(p, 1 / s, exp(m))),
  loglogistic = list(p = function(q, m, s, lt) stats::plogis((log(q) - m) / s, lower.tail = lt),
                     q = function(p, m, s) exp(m + s * stats::qlogis(p))),
  beta = list(p = function(q, m, s, lt) stats::pbeta(q, stats::plogis(m) / s, (1 - stats::plogis(m)) / s, lower.tail = lt),
              q = function(p, m, s) stats::qbeta(p, stats::plogis(m) / s, (1 - stats::plogis(m)) / s)))
crps_pred <- function(f, m, s, family = f) {
  F <- crps_fams[[f]]
  dist <- list(family = family, discrete = FALSE, lp_names = "mu", par_names = "logsd",
               d = function(y, lp, par, log = FALSE) stop("unused"),
               p = function(q, lp, par, lower.tail = TRUE) F$p(q, lp$mu, exp(par$logsd), lower.tail),
               q = function(prob, lp, par) F$q(prob, lp$mu, exp(par$logsd)),
               r = function(lp, par) stop("unused"), mean = function(lp, par) stop("unused"))
  lmn_new_pred(dist, list(mu = matrix(m, 1)), par = list(logsd = log(s)))
}
## the reference: integrate() over the whole support on the raw axis
crps_ref <- function(f, m, s, y) {
  F <- crps_fams[[f]]; S <- length(m)
  Fm <- function(z, lt) vapply(z, function(zz) mean(F$p(zz, m, s, lt)), 0)
  lo <- if (f == "gaussian") -Inf else 0; hi <- if (f == "beta") 1 else Inf
  stats::integrate(function(z) Fm(z, TRUE)^2, lo, y, rel.tol = 1e-12, subdivisions = 5000L)$value +
    stats::integrate(function(z) Fm(z, FALSE)^2, y, hi, rel.tol = 1e-12, subdivisions = 5000L)$value
}

test_that("every family of illume's is classified, so one added upstream fails here first", {
  fams <- eval(formals(illume::ilm_family)$family)
  known <- c(lmn_crps_families, lmn_discrete_families, lmn_no_dist_families)
  expect_true(all(fams %in% known), info = paste("unclassified:", paste(setdiff(fams, known), collapse = ", ")))
  expect_length(intersect(lmn_crps_families, lmn_discrete_families), 0L)
  ## every continuous family but the gaussian has its calibrated panel width
  ## (studies/findings/crps_continuous.md, run 4), and the gaussian its fast-path threshold
  expect_setequal(names(lmn_crps_tau_c), setdiff(lmn_crps_families, "gaussian"))
  expect_equal(lmn_crps_tau_c[c("lognormal", "weibull", "loglogistic", "beta")], c(lognormal = 1, weibull = 0.5, loglogistic = 1, beta = 1))
  expect_equal(lmn_crps_tau, 0.5)
})

test_that("a gaussian mixture's CRPS is exact, where the uniform quadrature was not", {
  skip_if_not_installed("scoringRules")
  set.seed(285)
  S <- 300
  m <- rnorm(S, 0, 0.5); s <- 10^runif(S, -2, 1)            # sds over three orders of magnitude
  p <- crps_pred("gaussian", m, s)
  for (y in c(-20, 0.001, 1, 20)) {
    ref <- scoringRules::crps_mixnorm(y, matrix(m, 1), matrix(s, 1))
    expect_equal(as.numeric(lmn_crps_continuous(p, y)), ref, tolerance = 1e-12)
  }
  ## the uniform quadrature is visibly off on this shape
  expect_gt(abs(as.numeric(lmn_crps_uniform(p, 1)) / scoringRules::crps_mixnorm(1, matrix(m, 1), matrix(s, 1)) - 1), 1e-5)
  ## a bimodal mixture, and unequal weights through ilm_scores()'s path
  m2 <- c(rep(0, 50), rep(20, 50)); s2 <- rep(1, 100)
  expect_equal(as.numeric(lmn_crps_continuous(crps_pred("gaussian", m2, s2), 10)),
               scoringRules::crps_mixnorm(10, matrix(m2, 1), matrix(s2, 1)), tolerance = 1e-12)
})

test_that("the positive families and beta match integrate() over the whole support", {
  set.seed(286)
  S <- 40
  cases <- list(lognormal = list(m = rnorm(S, 0.3, 0.3), s = exp(runif(S, log(0.1), log(2))), y = c(0.2, 1.5, 9)),
                weibull = list(m = rnorm(S, 0.3, 0.3), s = c(rep(1 / 0.3, S / 2), rep(0.5, S / 2)), y = c(0.05, 1.5, 9)),
                loglogistic = list(m = rnorm(S, 0, 0.3), s = exp(runif(S, log(0.1), log(1))), y = c(0.2, 1.5, 9)),
                beta = list(m = rnorm(S, 0, 0.5), s = c(rep(1 / 0.6, S / 2), rep(0.02, S / 2)), y = c(0.01, 0.5, 0.99)))
  for (f in names(cases)) {
    cs <- cases[[f]]; p <- crps_pred(f, cs$m, cs$s)
    for (y in cs$y)
      expect_equal(as.numeric(lmn_crps_continuous(p, y)), crps_ref(f, cs$m, cs$s, y), tolerance = 1e-9,
                   info = paste(f, y))
  }
})

test_that("a heavy loglogistic tail is integrated past the 1 - 1e-10 quantile", {
  ## s = 1.67, where stopping at that quantile loses 6.6e-4 of the CRPS
  p <- crps_pred("loglogistic", 0, 1.67)
  F <- function(u) stats::plogis(u / 1.67, log.p = TRUE)
  ref <- stats::integrate(function(u) exp(u + 2 * F(u)), -Inf, 0, rel.tol = 1e-13)$value +
    stats::integrate(function(u) exp(u + 2 * stats::plogis(-u / 1.67, log.p = TRUE)), 0, Inf, rel.tol = 1e-13)$value
  v <- lmn_crps_continuous(p, 1)
  expect_equal(as.numeric(v), ref, tolerance = 1e-10)
  expect_lt(attr(v, "tail_bound"), 1e-12 * ref)
})

test_that("an infinite CRPS is Inf, a failed integration NA, and an unknown family an error", {
  p <- crps_pred("loglogistic", c(0, 0.5), c(0.5, 2.2))
  expect_message(v <- lmn_crps_continuous(p, 1), "infinite")
  expect_identical(as.numeric(v), Inf)
  ## the cap forces integrate(), and a CDF that is not a number makes it fail
  bad <- crps_pred("lognormal", c(0, 0.2), c(0.5, 0.5))
  bad$dist$p <- function(q, lp, par, lower.tail = TRUE) rep(NaN, length(q))
  expect_message(v <- lmn_crps_continuous(bad, 1, cap = 0L), "could not be computed")
  expect_true(is.na(v))
  expect_identical(attr(v, "failed"), 1L)
  odd <- crps_pred("gaussian", 0, 1, family = "student")
  expect_error(lmn_crps_continuous(odd, 0), "no CRPS for a continuous 'student'")
})

test_that("where the family's quantiles are wrong, the row is integrated over the whole axis", {
  ## beta shape parameters near 1e-4, where R's qbeta() is inaccurate
  set.seed(3); S <- 20
  m <- c(rnorm(S / 2, -9, 0.5), rnorm(S / 2, 0, 0.5)); s <- rep(1.6, S)
  p <- crps_pred("beta", m, s)
  v <- suppressWarnings(lmn_crps_continuous(p, 0.3))
  expect_identical(attr(v, "integrated"), 1L)
  F <- crps_fams$beta
  lF <- function(u, lt) vapply(u, function(uu) log(mean(F$p(stats::plogis(uu), m, s, lt))), 0)
  lJ <- function(u) stats::plogis(u, log.p = TRUE) + stats::plogis(-u, log.p = TRUE)
  ref <- stats::integrate(function(u) exp(lJ(u) + 2 * lF(u, TRUE)), -Inf, stats::qlogis(0.3), rel.tol = 1e-12)$value +
    stats::integrate(function(u) exp(lJ(u) + 2 * lF(u, FALSE)), stats::qlogis(0.3), Inf, rel.tol = 1e-12)$value
  expect_equal(as.numeric(v), ref, tolerance = 1e-9)
})

test_that("a quantile is right even where the family's own quantiles are not", {
  ## beta shape parameters near 6e-5, where R's qbeta() is off by up to 0.5 in
  ## probability; the draws' quantiles then bracket the wrong place
  set.seed(3); S <- 20
  m <- c(rnorm(S / 2, -9, 0.5), rnorm(S / 2, 0, 0.5)); s <- rep(1.6, S)
  p <- crps_pred("beta", m, s)
  Fmix <- function(x) mean(crps_fams$beta$p(x, m, s, TRUE))
  for (prob in c(0.6, 0.75, 0.9)) {
    q <- suppressWarnings(lmn_pred_quantile(p, prob))
    expect_equal(Fmix(q), prob, tolerance = 1e-8, info = paste("prob", prob))
  }
  ## half the mass sits within 1e-300 of 0 (shape1 near 6e-5), so the 0.3
  ## quantile underflows: the answer is the bottom of the support to the
  ## bisection's tolerance (1e-12 absolute near 0), with F >= p
  q3 <- suppressWarnings(lmn_pred_quantile(p, 0.3))
  expect_lt(q3, 1e-11)
  expect_gte(Fmix(max(q3, .Machine$double.xmin)), 0.3)
  ## the draws' own quantiles really are wrong here, so this is the path tested
  Q <- suppressWarnings(as.vector(lmn_draw_quantiles(p, 1, 0.3)))
  Fd <- vapply(seq_len(S), function(k) crps_fams$beta$p(Q[k], m[k], s[k], TRUE), 0)
  expect_gt(max(abs(Fd - 0.3)), 0.1)
})

test_that("a real weibull AFT prediction is scored through the axis method", {
  set.seed(1); n <- 120
  d <- data.frame(x = rnorm(n)); tt <- rweibull(n, 1.5, exp(1 + 0.4 * d$x))
  d$time <- tt; d$event <- 1L
  fit <- suppressMessages(illume::ilm_model(time ~ x, data = d, family = "weibull",
                                            censor = illume::ilm_surv(d$time, d$event), verbose = FALSE))
  p <- ilm_predict(fit, d[1:5, ], target = "response", nsim = 30)
  v <- lmn_crps_continuous(p, tt[1:5])
  ## against integrate() of the mixture's own CDF, from illume's p()
  ref <- vapply(1:5, function(r) {
    Fm <- function(z, lt) vapply(z, function(zz) lmn_mix(p, r, matrix(zz), if (lt) "p" else "s")[1, 1], 0)
    stats::integrate(function(z) Fm(z, TRUE)^2, 0, tt[r], rel.tol = 1e-12)$value +
      stats::integrate(function(z) Fm(z, FALSE)^2, tt[r], Inf, rel.tol = 1e-12)$value
  }, 0)
  expect_equal(as.numeric(v), ref, tolerance = 1e-8)
})
