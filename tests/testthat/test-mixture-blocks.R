## The mixture is evaluated a block of draws per family call. These check that
## it gives what one call per draw gives, to rounding, across families.

ref_mix <- function(pred, i, X, what) {
  X <- as.matrix(X); f <- pred$dist[[what]]; out <- matrix(0, nrow(X), ncol(X))
  for (s in seq_along(pred$weights)) {
    a <- lmn_draw_args(pred, i, s, ncol(X))
    out <- out + pred$weights[s] * matrix(f(as.vector(X), a$lp, a$par), nrow(X))
  }
  out
}
ref_quant <- function(pred, i, prob) {
  vapply(seq_along(pred$weights), function(s) {
    a <- lmn_draw_args(pred, i, s); pred$dist$q(rep(prob, length(i)), a$lp, a$par)
  }, numeric(length(i)))
}
ref_mean <- function(pred) {
  i <- seq_len(nrow(pred$lp[[1]]))
  Reduce(`+`, lapply(seq_along(pred$weights), function(s) {
    a <- lmn_draw_args(pred, i, s); pred$weights[s] * pred$dist$mean(a$lp, a$par)
  }))
}
ref_log <- function(pred, y) {
  i <- seq_along(y)
  ld <- vapply(seq_along(pred$weights), function(s) {
    a <- lmn_draw_args(pred, i, s)
    pred$dist$d(y, a$lp, a$par, log = TRUE) + log(pred$weights[s])
  }, numeric(length(y)))
  mx <- apply(ld, 1L, max)
  -(mx + log(rowSums(exp(ld - mx))))
}

check_blocks <- function(p, y) {
  i <- seq_along(y)
  X <- outer(y, c(-1, 0, 2), `+`)
  expect_equal(lmn_mix(p, i, X, "p"), ref_mix(p, i, X, "p"), tolerance = 1e-12)
  expect_equal(lmn_mix(p, i, X, "d"), ref_mix(p, i, X, "d"), tolerance = 1e-12)
  expect_equal(lmn_draw_quantiles(p, i, 0.3), ref_quant(p, i, 0.3), tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_equal(lmn_score_log(p, y), ref_log(p, y), tolerance = 1e-12)
  if (is.null(p$dist$categories)) expect_equal(lmn_pred_mean(p), unname(ref_mean(p)), tolerance = 1e-12)
}

test_that("blocks of draws give what one draw at a time gives", {
  set.seed(11)
  d <- data.frame(g = factor(rep(1:12, each = 6)), x = rnorm(72))
  eta <- 0.4 + 0.5 * d$x + rnorm(12, 0, 0.4)[d$g]
  d$yc <- rpois(72, exp(eta)); d$yn <- rnbinom(72, mu = exp(eta), size = 3)
  d$yz <- ifelse(runif(72) < 0.2, 0L, d$yc)
  d$yg <- eta + rnorm(72, 0, exp(0.2 * d$x) * 0.5)
  nd <- d[c(1, 7, 13, 50), ]
  fits <- suppressMessages(list(
    poisson = illume::ilm_model(yc ~ x + (1 | g), data = d, family = "poisson", verbose = FALSE),
    nbinom = illume::ilm_model(yn ~ x + (1 | g), data = d, family = "nbinom", verbose = FALSE),
    zip = illume::ilm_model(yz ~ x + (1 | g), ziformula = ~ 1, data = d, family = "poisson",
                            verbose = FALSE),
    gaussian_disp = illume::ilm_model(yg ~ x + (1 | g), dispformula = ~ x, data = d,
                                      family = "gaussian", verbose = FALSE)))
  ys <- list(poisson = nd$yc, nbinom = nd$yn, zip = nd$yz, gaussian_disp = nd$yg)
  for (k in names(fits)) {
    p <- suppressWarnings(ilm_predict(fits[[k]], nd, target = "response", groups = "new",
                                      nsim = 300))
    check_blocks(p, ys[[k]])
  }
})

test_that("small blocks and one-draw blocks agree", {
  x <- suppressWarnings(illume::ilm_model(yc ~ x, data = data.frame(x = 1:30 / 10,
                          yc = rpois(30, exp(1:30 / 20))), family = "poisson", verbose = FALSE))
  p <- ilm_predict(x, data.frame(x = c(0.5, 2)), target = "response", nsim = 200)
  expect_gt(length(lmn_draw_blocks(p, 2)), 0L)
  expect_length(lmn_draw_blocks(p, lmn_block), 200L)       # a full block: one draw at a time
  expect_length(lmn_draw_blocks(p, 2), 1L)                 # a small one: every draw at once
})
