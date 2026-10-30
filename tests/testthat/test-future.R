hist_panel <- function() {
  set.seed(2)
  data.frame(site = rep(c("a", "b", "c"), times = c(4, 5, 3)),
             month = c(1:4, 1:5, 1:3),
             region = rep(c("north", "south", "north"), times = c(4, 5, 3)),
             rain = rnorm(12), y = rpois(12, 3))
}

test_that("each series gets the next h time points after its own last one", {
  d <- hist_panel()
  f <- suppressMessages(ilm_future(d, "month", "site", h = 2, response = "y"))
  expect_equal(nrow(f), 6L)
  expect_equal(f$site, rep(c("a", "b", "c"), each = 2))
  expect_equal(f$month, c(5, 6, 6, 7, 4, 5))
  expect_equal(f$horizon, rep(1:2, 3))
})

test_that("series-constant columns are copied; time-varying ones are left missing and named", {
  d <- hist_panel()
  expect_message(f <- ilm_future(d, "month", "site", h = 1, response = "y"),
                 "vary over time.*rain")
  expect_equal(f$region, c("north", "south", "north"))
  expect_true(all(is.na(f$rain)) && all(is.na(f$y)))
  expect_identical(attr(f, "unknown"), "rain")          # the response is not listed
})

test_that("monthly dates advance by calendar month, not by 30 days", {
  ## month-ends stay month-ends: 31 March is followed by 30 April, not 1 May
  d <- data.frame(day = as.Date(c("2026-01-31", "2026-02-28", "2026-03-31")),
                  y = 1:3)
  f <- ilm_future(d, "day", h = 2, response = "y")
  expect_equal(f$day, as.Date(c("2026-04-30", "2026-05-31")))
  ## any other day keeps its number, capped at the month's length
  m <- data.frame(day = as.Date(c("2025-11-15", "2025-12-15", "2026-01-15")), y = 1:3)
  expect_equal(ilm_future(m, "day", h = 2, response = "y")$day,
               as.Date(c("2026-02-15", "2026-03-15")))
  expect_equal(lmn_add_months(as.Date("2026-01-30"), 1:2),
               as.Date(c("2026-02-28", "2026-03-30")))
  expect_equal(lmn_add_months(as.Date("2027-11-30"), 3), as.Date("2028-02-29"))
  w <- data.frame(day = as.Date("2026-01-05") + 7 * 0:3, y = 1:4)
  expect_equal(ilm_future(w, "day", h = 1, response = "y")$day, as.Date("2026-02-02"))
  expect_equal(ilm_future(w, "day", h = 1, step = "quarter", response = "y")$day,
               as.Date("2026-04-26"))
})

test_that("irregular spacing needs an explicit step", {
  d <- data.frame(t = c(1, 2, 4, 7), y = 1:4)
  expect_error(ilm_future(d, "t", h = 1), "give `step`")
  expect_equal(ilm_future(d, "t", h = 2, step = 3, response = "y")$t, c(10, 13))
  dd <- data.frame(day = as.Date("2026-01-01") + c(0, 3, 10), y = 1:3)
  expect_error(ilm_future(dd, "day", h = 1), "give `step`")
})

test_that("bad requests stop with a reason", {
  d <- hist_panel()
  expect_error(ilm_future(d, "week", "site"), "not in `data`")
  expect_error(ilm_future(d, "month", "site", h = 0), "at least 1")
  expect_error(ilm_future(d, "region", "site"), "numeric, a Date or a date-time")
  d$horizon <- 1
  expect_error(ilm_future(d, "month", "site"), "already has a column called 'horizon'")
})
