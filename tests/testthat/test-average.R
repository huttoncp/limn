## Weighted averages of predicted means, checked draw by draw against the
## weighted average done by hand, and against ilm_aggregate(fun = "mean").

avg_fit <- function(seed = 1) {
  set.seed(seed); G <- 24; n <- 8
  d <- data.frame(g = factor(rep(sprintf("g%02d", 1:G), each = n)), x = rnorm(G * n),
                  age = rep(c("young", "old"), G * n / 2))
  d$y <- rpois(nrow(d), exp(0.4 + 0.3 * d$x + 0.5 * (d$age == "old") + rnorm(G, 0, 0.5)[d$g]))
  suppressMessages(suppressWarnings(illume::ilm_model(y ~ x + age + (1 | g), data = d,
                                                      family = "poisson", verbose = FALSE)))
}
## g = NA: a population average does not use the unit, but the rows need its column
avg_cells <- data.frame(region = rep(c("north", "south"), each = 2), age = rep(c("young", "old"), 2),
                        x = c(-0.5, -0.5, 0.5, 0.5), size = c(300, 700, 500, 500), g = NA)

test_that("a weighted average is the weighted mean of each draw's means", {
  p <- ilm_predict(avg_fit(), avg_cells, target = "mean", groups = "population", nsim = 200)
  a <- ilm_average(p, by = "region", weights = "size")
  expect_identical(a$kind, "mean")
  expect_equal(a$rows$region, c("north", "south"))
  expect_equal(a$rows$weight, c(1000, 1000))
  expect_equal(a$rows$n_rows, c(2L, 2L))
  w <- avg_cells$size
  north <- (w[1] * p$draws[1, ] + w[2] * p$draws[2, ]) / (w[1] + w[2])
  south <- (w[3] * p$draws[3, ] + w[4] * p$draws[4, ]) / (w[3] + w[4])
  expect_equal(unname(a$draws[1, ]), unname(north), tolerance = 1e-12)
  expect_equal(unname(a$draws[2, ]), unname(south), tolerance = 1e-12)
  expect_equal(a$estimate[1], sum(w[1:2] * p$estimate[1:2]) / 1000, tolerance = 1e-12)
  ## weights given as a vector are the same as the column, and only their ratios matter
  expect_equal(ilm_average(p, by = "region", weights = avg_cells$size / 100)$draws, a$draws,
               tolerance = 1e-12)
})

test_that("equal weights give ilm_aggregate()'s average of means", {
  p <- ilm_predict(avg_fit(), avg_cells, target = "mean", groups = "population", nsim = 100)
  expect_equal(ilm_average(p)$draws, ilm_aggregate(p, fun = "mean")$draws, tolerance = 1e-12)
  expect_equal(ilm_average(p, by = "region")$estimate,
               ilm_aggregate(p, by = "region", fun = "mean")$estimate, tolerance = 1e-12)
})

test_that("rows of different scenarios are never averaged together", {
  p <- ilm_predict(avg_fit(), list(low = transform(avg_cells, x = -1), high = transform(avg_cells, x = 1)),
                   target = "mean", groups = "population", nsim = 100)
  a <- ilm_average(p, weights = "size")
  expect_equal(nrow(a$rows), 2L)
  expect_true(all(a$estimate[a$rows$.scenario == "high"] > a$estimate[a$rows$.scenario == "low"]))
})

test_that("bad requests stop with a reason", {
  f <- avg_fit()
  p <- ilm_predict(f, avg_cells, target = "mean", groups = "population", nsim = 50)
  expect_error(ilm_average(ilm_predict(f, avg_cells, target = "response", groups = "population",
                                       nsim = 50)), "prediction of the mean")
  expect_error(ilm_average(p, weights = "pop"), "column the prediction's rows do not have")
  expect_error(ilm_average(p, weights = 1:3), "one number per row")
  expect_error(ilm_average(p, weights = c(1, -1, 1, 1)), "non-negative")
  expect_error(ilm_average(p, by = "region", weights = c(0, 0, 1, 1)), "sum to zero")
  expect_error(ilm_average(p, by = "zone"), "do not have: zone")
})
