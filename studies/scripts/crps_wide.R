## crps_wide: the discrete CRPS's wide-support path against the exact sum, and
## its speed against the path it replaces. The design is
## studies/designs/crps_wide.md, written before the method was built.
##
## Run from the package root, with LIMN_STUDY_LIB set to the pinned library
## (it needs tsibbledata for the PBS rows):
##   Rscript studies/scripts/crps_wide.R
source("studies/scripts/pinned_library.R")
tools::psnice(value = 19)
suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
version <- as.character(utils::packageVersion("limn"))
out_dir <- file.path("studies", "runs", version, "crps_wide")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
t_start <- Sys.time()

nb <- list(family = "nbinom", discrete = TRUE, lp_names = "mu", par_names = "logsize",
  d = function(y, lp, par, log = FALSE) stats::dnbinom(y, size = exp(par$logsize), mu = exp(lp$mu), log = log),
  p = function(q, lp, par) stats::pnbinom(q, size = exp(par$logsize), mu = exp(lp$mu)),
  q = function(prob, lp, par) stats::qnbinom(prob, size = exp(par$logsize), mu = exp(lp$mu)),
  r = function(lp, par) stats::rnbinom(length(lp$mu), size = exp(par$logsize), mu = exp(lp$mu)),
  mean = function(lp, par) exp(lp$mu))
po <- list(family = "poisson", discrete = TRUE, lp_names = "mu", par_names = character(0),
  d = function(y, lp, par, log = FALSE) stats::dpois(y, exp(lp$mu), log = log),
  p = function(q, lp, par) stats::ppois(q, exp(lp$mu)),
  q = function(prob, lp, par) stats::qpois(prob, exp(lp$mu)),
  r = function(lp, par) stats::rpois(length(lp$mu), exp(lp$mu)),
  mean = function(lp, par) exp(lp$mu))
S <- 200L

## a mixture of S draws for rows with log-means m, each draw shifted by `shift`
mix <- function(dist, m, shift, logsize = NULL) {
  lp <- list(mu = outer(m, shift, `+`))
  lmn_new_pred(dist, lp, par = if (is.null(logsize)) list() else list(logsize = rep(logsize, length(shift))))
}
## the support the method uses, from the mixture's 1e-6 and 1 - 1e-6 quantiles
support <- function(p) {
  n <- lmn_pred_n(p); i <- seq_len(n)
  qlo <- lmn_draw_quantiles(p, i, 1e-6); qhi <- lmn_draw_quantiles(p, i, 1 - 1e-6)
  klo <- lmn_mix_qint(p, i, 1e-6, apply(qlo, 1, min), apply(qlo, 1, max))
  khi <- lmn_mix_qint(p, i, 1 - 1e-6, apply(qhi, 1, min), apply(qhi, 1, max))
  data.frame(klo = klo, khi = khi, width = khi - klo + 1)
}
## score rows by both paths: wide (the method) and exact (every integer)
score <- function(case, p, y) {
  w <- lmn_crps_discrete(p, y)
  ex <- as.numeric(lmn_crps_discrete(p, y, max_support = 1e9))
  sp <- support(p)
  data.frame(case = case, width = sp$width, y = y, exact = ex, estimate = as.numeric(w),
             abs_error = abs(as.numeric(w) - ex), error_estimate = attr(w, "error_estimate"),
             bracket = attr(w, "bound"), flagged = attr(w, "flagged"), step = attr(w, "step"),
             wide = sp$width > 2000)
}
yq <- function(p, probs) vapply(seq_len(lmn_pred_n(p)), function(i)
  lmn_pred_quantile(lmn_subset_pred(p, i), probs[(i - 1) %% length(probs) + 1]), 0)

set.seed(20260929)
rows <- list()
## a row's log-mean by bisection so that its support has width W, for the
## draws' shifts and size given (the addendum)
mean_for_width <- function(dist, W, shift, logsize = NULL) {
  lo <- log(5); hi <- log(5e6)
  for (it in 1:60) {
    mid <- (lo + hi) / 2
    if (support(mix(dist, mid, shift, logsize))$width < W) lo <- mid else hi <- mid
  }
  hi
}
## (A) typical: NB draws, sizes 5 to 50, log-means spread SD 0.05, each row's
## mean chosen for a target width of 1,000 to 10,000 (the addendum)
dropped_A <- 0L
for (sz in c(5, 20, 50)) {
  shiftA <- stats::rnorm(S, 0, 0.05)
  m <- vapply(c(1000, 2500, 4000, 6000, 8000, 10000), mean_for_width, 0, dist = nb, shift = shiftA,
              logsize = log(sz))
  p <- mix(nb, m, shiftA, log(sz))
  rows[[length(rows) + 1L]] <- score(sprintf("A: nbinom size %g", sz), p, yq(p, c(0.1, 0.5, 0.9)))
}
## (B) draws that differ widely: bimodal (m and 2m), and heavy-tailed log-means
m <- log(seq(1500, 6000, length.out = 10))
pb <- mix(nb, m, c(rep(0, S / 2), rep(log(2), S / 2)), log(30))
rows[[length(rows) + 1L]] <- score("B: bimodal draws", pb, yq(pb, c(0.1, 0.5, 0.9)))
shiftH <- 0.3 * stats::rt(S, 3)
mH <- vapply(seq(1000, 10000, length.out = 10), mean_for_width, 0, dist = nb, shift = shiftH, logsize = log(30))
ph <- mix(nb, mH, shiftH, log(30))
keep <- support(ph)$width <= 12000
dropped_B <- sum(!keep)
if (any(keep)) {
  ph <- lmn_subset_pred(ph, which(keep))
  rows[[length(rows) + 1L]] <- score("B: heavy-tailed draws", ph, yq(ph, c(0.1, 0.5, 0.9)))
}
## (C) y far in the tails, outside the support, and at its ends
m <- log(c(3000, 8000, 15000, 25000))
pc <- mix(nb, m, stats::rnorm(S, 0, 0.05), log(20))
sp <- support(pc)
yc <- c(yq(pc, 1e-5), yq(pc, 1 - 1e-5), sp$klo - 50, sp$khi + 50, sp$klo, sp$khi)
pc6 <- lmn_subset_pred(pc, rep(seq_len(4), 6))
rows[[length(rows) + 1L]] <- score("C: y in the far tails, outside and at the ends", pc6, yc)
## (D) steep F: near-Poisson NB (size 1e4) and Poisson itself, widths 2,000 to 5,000
m <- log(seq(40000, 250000, length.out = 8))
pd <- mix(nb, m, stats::rnorm(S, 0, 0.002), log(1e4))
rows[[length(rows) + 1L]] <- score("D: near-Poisson (size 1e4)", pd, yq(pd, c(0.1, 0.5, 0.9)))
pp <- mix(po, m, stats::rnorm(S, 0, 0.002))
rows[[length(rows) + 1L]] <- score("D: Poisson", pp, yq(pp, c(0.1, 0.5, 0.9)))
## (E) the threshold: means tuned to widths of 1,990 to 2,010
shiftE <- stats::rnorm(S, 0, 0.02)
width_at <- function(lm) support(mix(nb, lm, shiftE, log(20)))$width
target <- 1990:2010
lmE <- vapply(target, function(W) {
  lo <- log(50); hi <- log(20000)
  for (it in 1:60) { mid <- (lo + hi) / 2; if (width_at(mid) < W) lo <- mid else hi <- mid }
  hi
}, 0)
pe <- mix(nb, lmE, shiftE, log(20))
rows[[length(rows) + 1L]] <- score("E: widths either side of 2,000", pe, yq(pe, c(0.1, 0.5, 0.9)))
res <- do.call(rbind, rows)

## ---- speed on the dry run's PBS rows: the current path against the old -----------
old_env <- new.env(parent = asNamespace("limn"))
old_src <- system("git show main:R/scores.R", intern = TRUE)
eval(parse(text = old_src), envir = old_env)
stopifnot(identical(formals(old_env$lmn_crps_discrete)$max_support, 1e4),
          identical(formals(lmn_crps_discrete)$max_support, 2000),
          exists("lmn_crps_bracket", envir = old_env, inherits = FALSE))
pbs <- as.data.frame(tsibbledata::PBS)
pbs$series <- paste(pbs$Concession, pbs$Type, pbs$ATC1, pbs$ATC2, sep = "|")
months <- sort(unique(pbs$Month)); pbs$t <- match(pbs$Month, months)
full <- names(which(tapply(pbs$t, pbs$series, length) == 204))
speed <- list()
for (s in full[c(10, 40, 120)]) {
  d <- pbs[pbs$series == s, ]; d$series <- factor(d$series)
  d$s12 <- sin(2*pi*d$t/12); d$c12 <- cos(2*pi*d$t/12); d$s6 <- sin(2*pi*d$t/6); d$c6 <- cos(2*pi*d$t/6)
  tr <- d[d$t <= 192, ]; te <- d[d$t > 192, ]
  f <- suppressWarnings(suppressMessages(illume::ilm_model(Scripts ~ s12 + c12 + s6 + c6, data = tr,
         family = "nbinom", ar = illume::ilm_rw1(~ t | series), verbose = FALSE)))
  fc <- suppressWarnings(ilm_forecast(f, te[c("series", "t", "s12", "c12", "s6", "c6")], nsim = 1000, seed = 1))
  t0 <- proc.time()[["elapsed"]]; new <- lmn_crps_discrete(fc, te$Scripts); t_new <- proc.time()[["elapsed"]] - t0
  t0 <- proc.time()[["elapsed"]]; old <- old_env$lmn_crps_discrete(fc, te$Scripts); t_old <- proc.time()[["elapsed"]] - t0
  speed[[length(speed) + 1L]] <- data.frame(series = s, mean_scripts = mean(te$Scripts), rows = nrow(te),
    median_width = stats::median(support(fc)$width), seconds_old = t_old, seconds_new = t_new,
    max_abs_diff = max(abs(as.numeric(new) - as.numeric(old))),
    within = all(abs(as.numeric(new) - as.numeric(old)) <= attr(old, "bound") + attr(new, "error_estimate") + 1e-9),
    max_flagged = sum(attr(new, "flagged")))
}
speed <- do.call(rbind, speed)
minutes <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
saveRDS(list(rows = res, speed = speed), file.path(out_dir, "crps_wide.rds"))

## ---- the pass rules and the write-up ----------------------------------------------
wide <- res[res$wide, ]
rule1 <- all(res$abs_error <= 1e-6 * res$exact)
rule2 <- all(wide$abs_error <= wide$error_estimate)
under <- wide[wide$abs_error > wide$error_estimate, ]
e_rows <- res[res$case == "E: widths either side of 2,000", ]
rule3 <- all(e_rows$abs_error <= 1e-6 * e_rows$exact) && any(e_rows$width <= 2000) && any(e_rows$width > 2000)
inside <- mean(wide$abs_error <= wide$bracket)
summ <- do.call(rbind, lapply(split(res, res$case), function(z) data.frame(
  case = z$case[1], rows = nrow(z), wide_rows = sum(z$wide), widths = sprintf("%d to %d", min(z$width), max(z$width)),
  max_rel_error = max(z$abs_error / z$exact), max_error_over_estimate = if (any(z$wide))
    max(z$abs_error[z$wide] / pmax(z$error_estimate[z$wide], 1e-300)) else NA,
  flagged = sum(z$flagged))))
write.csv(summ, file.path(out_dir, "cases.csv"), row.names = FALSE)
write.csv(speed, file.path(out_dir, "speed.csv"), row.names = FALSE)
f3 <- function(v) ifelse(is.na(v), "", formatC(v, format = "e", digits = 1))
md <- c(
  "# The discrete CRPS on wide supports", "",
  sprintf("## limn %s (%s; illume %s)", version, format(Sys.Date()), utils::packageVersion("illume")), "",
  "Generated by `studies/scripts/crps_wide.R`, to the design in `studies/designs/crps_wide.md`. Do not edit by hand.", "",
  lmn_study_header(), "",
  sprintf("Run on limn %s: %.1f minutes on 1 core.", system("git rev-parse --short HEAD", intern = TRUE), minutes), "",
  "### Pass rules (fixed in the design)", "",
  sprintf("1. |estimate - exact| <= 1e-6 x exact in every row: **%s**.", if (rule1) "pass" else "FAIL"),
  sprintf("2. |estimate - exact| <= error_estimate in every wide-path row: **%s** (%d of %d rows).",
          if (rule2) "pass" else "FAIL", sum(wide$abs_error <= wide$error_estimate), nrow(wide)),
  sprintf("3. Rows either side of 2,000 agree with the exact sum: **%s** (widths %d to %d).",
          if (rule3) "pass" else "FAIL", min(e_rows$width), max(e_rows$width)),
  sprintf("4. The error lies within the rigorous bracket in %.3f of wide-path rows (reported, not required).", inside), "",
  if (nrow(under)) c("Rows where the estimate undershoots:", "",
    sprintf("- %s: width %d, y %g, error %s against estimate %s", under$case, under$width, under$y,
            f3(under$abs_error), f3(under$error_estimate)), "") else character(0),
  sprintf("Rows left out for widths over 12,000 (the exact sum's limit): %d in (B, heavy-tailed).", dropped_B), "",
  "### By case", "",
  "| case | rows | wide-path rows | widths | max relative error | max error / estimate | flagged |",
  "|---|---|---|---|---|---|---|",
  sprintf("| %s | %d | %d | %s | %s | %s | %d |", summ$case, summ$rows, summ$wide_rows, summ$widths,
          f3(summ$max_rel_error), f3(summ$max_error_over_estimate), summ$flagged), "",
  "### Speed on the dry run's PBS rows (1,000 draws, 12 rows a series)", "",
  "| series | mean scripts | median width | seconds, old path | seconds, new | max difference | within old bound + new estimate | flagged |",
  "|---|---|---|---|---|---|---|---|",
  sprintf("| %s | %.0f | %.0f | %.1f | %.1f | %s | %s | %d |", speed$series, speed$mean_scripts, speed$median_width,
          speed$seconds_old, speed$seconds_new, f3(speed$max_abs_diff), speed$within, speed$max_flagged))
writeLines(md, file.path("studies", "findings", "crps_wide.md"))
cat(md, sep = "\n")
