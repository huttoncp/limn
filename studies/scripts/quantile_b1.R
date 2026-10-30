## quantile_b1: the capped Illinois root-finder for the mixture quantile, checked
## as registered in studies/designs/quantile_b1.md (with its addendum).
##
## Run from the package root, with LIMN_STUDY_LIB set to a pinned library:
##   Rscript studies/scripts/quantile_b1.R [pilot]
## "pilot" runs only the pilot: the gaussian calibration at 2,000 x 1,000, both
## ways, and one call at S = 4,000.
source("studies/scripts/pinned_library.R")
tools::psnice(value = 19)
Sys.setenv(OMP_NUM_THREADS = "1")
pilot <- identical(commandArgs(TRUE)[1], "pilot")
suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
ns <- asNamespace("limn")
version <- as.character(utils::packageVersion("limn"))
out_dir <- file.path("studies", "runs", version, "quantile_b1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
t_start <- Sys.time()
quiet <- function(e) suppressMessages(suppressWarnings(e))
stamp <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

## ---- the reference: main 3841634's bracket check and bisection, verbatim ---------
ref_widen_bracket <- function(pred, i, prob, lo, hi, discrete) {
  Fx <- function(rows, x) lmn_mix(pred, rows, matrix(x), "p")[, 1]
  lo_ok <- function(rows) { v <- Fx(rows, if (discrete) lo[rows_k(rows)] - 1 else lo[rows_k(rows)])
    if (discrete) v < prob else v <= prob }
  rows_k <- function(rows) match(rows, i)
  hi_ok <- function(rows) Fx(rows, hi[rows_k(rows)]) >= prob
  bad_lo <- !lo_ok(i); bad_hi <- !hi_ok(i)
  bad_lo[is.na(bad_lo)] <- TRUE; bad_hi[is.na(bad_hi)] <- TRUE
  step <- pmax(hi - lo, 1e-8 * pmax(1, abs(lo), abs(hi)), if (discrete) 1 else 0)
  for (it in seq_len(1100L)) {
    if (!any(bad_lo | bad_hi)) break
    if (any(bad_lo)) {
      lo[bad_lo] <- lo[bad_lo] - step[bad_lo]
      if (discrete) lo[bad_lo] <- floor(lo[bad_lo])
      ok <- lo_ok(i[bad_lo]); ok[is.na(ok)] <- FALSE; bad_lo[bad_lo] <- !ok
    }
    if (any(bad_hi)) {
      hi[bad_hi] <- hi[bad_hi] + step[bad_hi]
      if (discrete) hi[bad_hi] <- ceiling(hi[bad_hi])
      ok <- hi_ok(i[bad_hi]); ok[is.na(ok)] <- FALSE; bad_hi[bad_hi] <- !ok
    }
    step <- 2 * step
  }
  stuck <- bad_lo | bad_hi | !is.finite(lo) | !is.finite(hi)
  if (any(stuck)) {
    message("a quantile could not be bracketed for ", sum(stuck), " row(s): ",
            "their predictive CDF is not a number there, so the quantile is NA.")
    lo[stuck] <- NA_real_; hi[stuck] <- NA_real_
  }
  list(lo = lo, hi = hi)
}
ref_pred_quantile <- function(pred, prob) {
  if (lmn_is_empirical(pred))
    return(lmn_rescale(pred, apply(pred$draws, 1L, stats::quantile, probs = prob, names = FALSE,
                                   type = if (isTRUE(pred$dist$discrete)) 1L else 7L)))
  n <- nrow(pred$lp[[1]])
  out <- numeric(n)
  for (i in lmn_row_chunks(n)) {
    Q <- lmn_draw_quantiles(pred, i, prob)
    lo <- apply(Q, 1L, min); hi <- apply(Q, 1L, max)
    if (pred$dist$discrete) {
      out[i] <- lmn_mix_qint(pred, i, prob, lo, hi)
    } else {
      b <- ref_widen_bracket(pred, i, prob, lo, hi, FALSE); lo <- b$lo; hi <- b$hi
      for (it in seq_len(200L)) {
        open <- !is.na(lo) & (hi - lo) > 1e-12 * pmax(1, abs(lo), abs(hi))
        if (!any(open)) break
        mid <- (lo + hi) / 2
        Fm <- lmn_mix(pred, i[open], matrix(mid[open]), "p")[, 1]
        up <- Fm >= prob
        hi[open][up] <- mid[open][up]
        lo[open][!up] <- mid[open][!up]
      }
      out[i] <- (lo + hi) / 2
    }
  }
  lmn_rescale(pred, out)
}
## evaluated in limn's namespace, so both methods' lmn_mix() calls are traced alike
environment(ref_widen_bracket) <- ns
environment(ref_pred_quantile) <- ns
new_pred_quantile <- get("lmn_pred_quantile", envir = ns)

## the mixture-CDF passes a call takes
passes <- 0L
trace("lmn_mix", quote(passes <<- passes + 1L), print = FALSE, where = ns)
counted <- function(f, ...) { passes <<- 0L; v <- f(...); list(value = v, passes = passes) }

## ---- families and mixtures, as in crps_continuous ------------------------------------
fams <- list(
  gaussian = list(p = function(q, m, s, lt) stats::pnorm(q, m, s, lower.tail = lt), q = function(p, m, s) stats::qnorm(p, m, s)),
  lognormal = list(p = function(q, m, s, lt) stats::plnorm(q, m, s, lower.tail = lt), q = function(p, m, s) stats::qlnorm(p, m, s)),
  weibull = list(p = function(q, m, s, lt) stats::pweibull(q, 1 / s, exp(m), lower.tail = lt),
                 q = function(p, m, s) stats::qweibull(p, 1 / s, exp(m))),
  loglogistic = list(p = function(q, m, s, lt) stats::plogis((log(q) - m) / s, lower.tail = lt),
                     q = function(p, m, s) exp(m + s * stats::qlogis(p))),
  beta = list(p = function(q, m, s, lt) stats::pbeta(q, stats::plogis(m) / s, (1 - stats::plogis(m)) / s, lower.tail = lt),
              q = function(p, m, s) stats::qbeta(p, stats::plogis(m) / s, (1 - stats::plogis(m)) / s)))
mkdist <- function(f) {
  F <- fams[[f]]
  list(family = f, discrete = FALSE, lp_names = "mu", par_names = "logsd",
       d = function(y, lp, par, log = FALSE) stop("unused"),
       p = function(q, lp, par, lower.tail = TRUE) F$p(q, lp$mu, exp(par$logsd), lower.tail),
       q = function(prob, lp, par) F$q(prob, lp$mu, exp(par$logsd)),
       r = function(lp, par) stop("unused"), mean = function(lp, par) stop("unused"))
}
mkpred <- function(f, m, s, w) lmn_new_pred(mkdist(f), list(mu = matrix(m, 1)), par = list(logsd = log(s)), weights = w)
lu <- function(n, a, b) exp(stats::runif(n, log(a), log(b)))
rand_mix <- function(f) {
  S <- round(lu(1, 2, 1000)); w <- stats::rexp(S); w <- w / sum(w)
  b_lim <- switch(f, gaussian = c(1e-3, 10), lognormal = c(0.05, 3), weibull = c(0.1, 3.4),
                  loglogistic = c(0.05, 1.95), beta = c(1 / 1000, 1 / 0.6))
  t_max <- if (f == "gaussian") 30 else 3
  m <- stats::rnorm(1) + lu(1, 0.01, t_max) * stats::rnorm(S)
  a <- stats::runif(1, 0, log(1000) / 2)
  s <- pmin(pmax(lu(1, b_lim[1], b_lim[2]) * exp(stats::runif(S, -a, a)), b_lim[1]), b_lim[2])
  list(m = m, s = s, w = w)
}
shapes <- function(f, S) {
  base <- switch(f, gaussian = 1, lognormal = 0.5, weibull = 0.5, loglogistic = 0.3, beta = 0.05)
  tiny <- switch(f, gaussian = 1e-4, lognormal = 1e-4, weibull = 1e-3, loglogistic = 1e-3, beta = 1e-5)
  sep <- if (f == "gaussian") 20 else 20 * base
  out <- list(
    baseline = list(m = stats::rnorm(S, 0, 0.3 * base), s = rep(base, S)),
    a_bimodal = list(m = c(rep(0, S %/% 2), rep(sep, S - S %/% 2)), s = rep(base, S)),
    b_scale_spread = list(m = stats::rnorm(S, 0, 0.5 * base),
                          s = pmin(base * 10^stats::runif(S, -1.5, 1.5), if (f == "loglogistic") 1.95 else if (f == "beta") 1.6 else if (f == "lognormal") 3 else Inf)),
    c_tail = list(m = stats::rnorm(S, 0, 0.3 * base), s = rep(base, S)),
    d1_degenerate = list(m = stats::rnorm(S, 0, tiny), s = rep(tiny, S)),
    d2_spiky = list(m = stats::rnorm(S, 0, if (f == "gaussian") 1 else 0.3), s = rep(tiny, S)))
  if (f != "gaussian") {
    add <- switch(f,
      weibull = list(weibull_shape0.3 = list(m = stats::rnorm(S, 0, 0.1), s = rep(1 / 0.3, S))),
      beta = list(beta_0.3_0.3 = list(m = stats::rnorm(S, 0, 0.05), s = rep(1 / 0.6, S))),
      loglogistic = lapply(c(s0.8 = 0.8, s1.2 = 1.2, s1.6 = 1.6, s1.9 = 1.9),
                           function(sv) list(m = stats::rnorm(S, 0, 0.1), s = rep(sv, S))),
      lognormal = list(sdlog3 = list(m = stats::rnorm(S, 0, 0.1), s = rep(3, S)),
                       sdlog15_spread = list(m = stats::rnorm(S, 0, 0.25), s = 0.5 * 10^stats::runif(S, -1.5, 1.5))))
    mixed <- lapply(add, function(a) list(m = a$m, s = c(a$s[seq_len(S %/% 2)], rep(base, S - S %/% 2))))
    names(mixed) <- paste0(names(add), "_mixed")
    out <- c(out, add, mixed)
  }
  lapply(out, function(x) { x$w <- rep(1 / S, S); x })
}
probs <- c(1e-6, 0.001, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.999, 1 - 1e-6)

## one mixture at every probability, both ways: agreement, passes, the bracket property
compare_one <- function(f, nm, mx) {
  pr <- mkpred(f, mx$m, mx$s, mx$w)
  do.call(rbind, lapply(probs, function(p) {
    a <- quiet(counted(ref_pred_quantile, pr, p)); b <- quiet(counted(new_pred_quantile, pr, p))
    qa <- a$value; qb <- b$value
    tol <- 1e-12 * max(1, abs(qa), na.rm = TRUE)
    Fm <- function(x) lmn_mix(pr, 1L, matrix(x), "p")[1, 1]
    passes_saved <- passes
    lowF <- if (is.finite(qb)) Fm(qb - 2 * tol) else NA; highF <- if (is.finite(qb)) Fm(qb + 2 * tol) else NA
    flat <- is.finite(qb) && abs(Fm(qb) - p) < 1e-15 && abs(highF - lowF) < 1e-15
    passes <<- passes_saved
    data.frame(family = f, case = nm, S = length(mx$m), prob = p,
               agree = if (is.na(qa) && is.na(qb)) TRUE else if (xor(is.na(qa), is.na(qb))) FALSE
                       else abs(qb - qa) <= 2e-12 * max(1, abs(qa)),
               rel_diff = if (is.na(qa) || is.na(qb)) NA else abs(qb - qa) / max(1, abs(qa)),
               passes_ref = a$passes, passes_new = b$passes,
               bracket_ok = is.finite(qb) && lowF < p && p <= highF, flat = flat)
  }))
}

## ---- the profiling data -------------------------------------------------------------
set.seed(1)
G <- 50; n <- 2000
d <- data.frame(g = factor(sample(sprintf("g%02d", 1:G), n, TRUE)), x = rnorm(n))
d$yg <- 1 + 0.5 * d$x + rnorm(G, 0, 0.5)[as.integer(d$g)] + rnorm(n)
fg <- quiet(illume::ilm_model(yg ~ x + (1 | g), data = d, verbose = FALSE))
pg <- quiet(ilm_predict(fg, d, target = "response", groups = "fitted", nsim = 1000))
## run f() with the package's lmn_pred_quantile set to `which` (today's or B1's)
with_quantile <- function(which, f) {
  assignInNamespace("lmn_pred_quantile", if (which == "today") ref_pred_quantile else new_pred_quantile, ns = "limn")
  on.exit(assignInNamespace("lmn_pred_quantile", new_pred_quantile, ns = "limn"))
  t <- system.time(v <- quiet(f()))[["elapsed"]]
  list(value = v, seconds = t)
}

if (pilot) {
  stamp("pilot: gaussian calibration at 2,000 x 1,000, both ways")
  cal_new <- with_quantile("B1", function() illume::ilm_calibration(pg, d$yg))
  cal_old <- with_quantile("today", function() illume::ilm_calibration(pg, d$yg))
  cat(sprintf("calibration: today %.1f s, B1 %.1f s (%.2fx); widths max rel diff %.1e; verdict identical %s\n",
              cal_old$seconds, cal_new$seconds, cal_old$seconds / cal_new$seconds,
              max(abs(cal_new$value$coverage$width - cal_old$value$coverage$width) / cal_old$value$coverage$width),
              identical(cal_new$value$verdict, cal_old$value$verdict)))
  stamp("pilot: one call at S = 4,000")
  set.seed(2950)
  mx <- shapes("gaussian", 4000)$baseline
  t1 <- system.time(r1 <- compare_one("gaussian", "baseline", mx))[["elapsed"]]
  cat(sprintf("S = 4,000, 13 probabilities both ways: %.1f s; passes today median %d, B1 median %d; all agree %s\n",
              t1, as.integer(stats::median(r1$passes_ref)), as.integer(stats::median(r1$passes_new)), all(r1$agree)))
  cat("PILOT DONE\n")
  quit(save = "no")
}

## ---- checks 1 to 3: agreement, the cap, the bracket property ---------------------------
res <- list()
rows <- list()
for (f in names(fams)) {
  stamp("checks 1-3:", f)
  set.seed(2950)
  for (S in c(100, 1000, 4000)) {
    sh <- shapes(f, S)
    for (nm in names(sh)) rows[[length(rows) + 1L]] <- compare_one(f, nm, sh[[nm]])
  }
}
for (f in names(fams)) {
  stamp("checks 1-3, random:", f)
  set.seed(2951)
  for (j in seq_len(200)) rows[[length(rows) + 1L]] <- compare_one(f, "random", rand_mix(f))
}
c13 <- do.call(rbind, rows)
res$c13 <- c13

## ---- check 4: downstream, on the profiling data -----------------------------------------
stamp("check 4 and 6")
tt <- rweibull(n, 1.5, exp(1 + 0.4 * d$x)); cc <- rexp(n, 0.15)
d$time <- pmin(tt, cc); d$event <- as.integer(tt <= cc)
fw <- quiet(illume::ilm_model(time ~ x, data = d, family = "weibull", censor = illume::ilm_surv(d$time, d$event), verbose = FALSE))
pw <- quiet(ilm_predict(fw, d, target = "response", nsim = 1000))
arm <- function(f) list(today = with_quantile("today", f), B1 = with_quantile("B1", f))
a_gcal <- arm(function() illume::ilm_calibration(pg, d$yg))
a_gq <- arm(function() ilm_scores(pg, d$yg, rules = c("pinball", "interval")))
a_wcal <- arm(function() illume::ilm_calibration(pw, tt))
reldiff <- function(a, b) max(abs(a - b) / pmax(abs(b), .Machine$double.xmin), na.rm = TRUE)
c4 <- data.frame(
  result = c("gaussian calibration widths", "gaussian calibration coverage", "gaussian verdict text",
             "gaussian pinball and interval scores", "weibull calibration widths", "weibull calibration coverage",
             "weibull verdict text"),
  value = c(sprintf("%.1e", reldiff(a_gcal$B1$value$coverage$width, a_gcal$today$value$coverage$width)),
            sprintf("%.1e", reldiff(a_gcal$B1$value$coverage$coverage, a_gcal$today$value$coverage$coverage)),
            identical(a_gcal$B1$value$verdict, a_gcal$today$value$verdict),
            sprintf("%.1e", reldiff(as.matrix(a_gq$B1$value[attr(a_gq$B1$value, "rules")]),
                                    as.matrix(a_gq$today$value[attr(a_gq$today$value, "rules")]))),
            sprintf("%.1e", reldiff(a_wcal$B1$value$coverage$width, a_wcal$today$value$coverage$width)),
            sprintf("%.1e", reldiff(a_wcal$B1$value$coverage$coverage, a_wcal$today$value$coverage$coverage)),
            identical(a_wcal$B1$value$verdict, a_wcal$today$value$verdict)))
c4$pass <- c(reldiff(a_gcal$B1$value$coverage$width, a_gcal$today$value$coverage$width) <= 1e-10,
             reldiff(a_gcal$B1$value$coverage$coverage, a_gcal$today$value$coverage$coverage) <= 1e-10,
             identical(a_gcal$B1$value$verdict, a_gcal$today$value$verdict),
             reldiff(as.matrix(a_gq$B1$value[attr(a_gq$B1$value, "rules")]),
                     as.matrix(a_gq$today$value[attr(a_gq$today$value, "rules")])) <= 1e-10,
             reldiff(a_wcal$B1$value$coverage$width, a_wcal$today$value$coverage$width) <= 1e-10,
             reldiff(a_wcal$B1$value$coverage$coverage, a_wcal$today$value$coverage$coverage) <= 1e-10,
             identical(a_wcal$B1$value$verdict, a_wcal$today$value$verdict))
cost <- data.frame(arm = c("gaussian calibration", "gaussian pinball + interval", "weibull calibration"),
                   today = c(a_gcal$today$seconds, a_gq$today$seconds, a_wcal$today$seconds),
                   B1 = c(a_gcal$B1$seconds, a_gq$B1$seconds, a_wcal$B1$seconds))
cost$ratio <- cost$B1 / cost$today
res$c4 <- c4; res$cost <- cost
untrace("lmn_mix", where = ns)
minutes <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
saveRDS(res, file.path(out_dir, "results.rds"))

## ---- findings -------------------------------------------------------------------------
yes <- function(ok) if (isTRUE(ok)) "yes" else "NO"
nf <- !c13$flat
checks <- data.frame(
  check = c("1 agreement with the bisection, every row (max rel diff)", "1 rows that disagree",
            "2 passes at most the bisection's + 2, every call (worst excess)",
            "3 bracket property, every non-flat row (failures)", "3 flat-at-p rows (reported apart)",
            "4 downstream results unchanged (all rows of the table below)"),
  value = c(sprintf("%.1e", max(c13$rel_diff, na.rm = TRUE)), sum(!c13$agree),
            max(c13$passes_new - c13$passes_ref), sum(!c13$bracket_ok[nf]), sum(c13$flat),
            sprintf("%d of %d", sum(c4$pass), nrow(c4))),
  rule = c("<= 2e-12", "0", "<= 2", "0", "reported", "all"),
  pass = c(yes(all(c13$agree)), yes(sum(!c13$agree) == 0), yes(max(c13$passes_new - c13$passes_ref) <= 2),
           yes(all(c13$bracket_ok[nf])), "", yes(all(c4$pass))))
write.csv(checks, file.path(out_dir, "checks.csv"), row.names = FALSE)
pass_tab <- stats::aggregate(cbind(passes_ref, passes_new) ~ family + case, data = c13, FUN = stats::median)
write.csv(pass_tab, file.path(out_dir, "passes.csv"), row.names = FALSE)
write.csv(cost, file.path(out_dir, "cost.csv"), row.names = FALSE)
merge_ok <- all(c13$agree) && max(c13$passes_new - c13$passes_ref) <= 2 && all(c13$bracket_ok[nf]) && all(c4$pass)
md <- c("# A faster quantile root-finder, capped: checks", "",
        sprintf("## limn %s (%s; illume %s)", version, format(Sys.Date()), utils::packageVersion("illume")), "",
        "Generated by `studies/scripts/quantile_b1.R`, to the design in `studies/designs/quantile_b1.md` (with its addendum). Do not edit by hand.", "",
        lmn_study_header(), "",
        sprintf("Run on limn %s: %.1f minutes on 1 core. %d quantile calls compared, both ways.",
                system("git rev-parse --short HEAD", intern = TRUE), minutes, nrow(c13)), "",
        "| check | result | rule | pass |", "|---|---|---|---|",
        sprintf("| %s | %s | %s | %s |", checks$check, checks$value, checks$rule, checks$pass), "",
        "Check 5, the test suite, is run separately and reported with the push request.", "",
        sprintf("**Decision rule** (checks 1, 2, 4 and 5, and 3 on every non-flat row): %s.",
                if (merge_ok) "checks 1 to 4 pass; B1 merges if the suite passes" else "a check FAILED; B1 does not merge"), "",
        "### Check 4: downstream results, today against B1", "",
        "| result | difference | pass |", "|---|---|---|", sprintf("| %s | %s | %s |", c4$result, c4$value, ifelse(c4$pass, "yes", "NO")), "",
        "### Passes per call, median by family and shape (today, B1)", "",
        "| family | shape | today | B1 |", "|---|---|---|---|",
        sprintf("| %s | %s | %.0f | %.0f |", pass_tab$family, pass_tab$case, pass_tab$passes_ref, pass_tab$passes_new), "",
        "**Why the gain shrinks with S.** The cap allows each row 2 passes beyond the bisection's count, and a row may take a secant step only while, were that step to gain nothing, bisection could still finish within the budget. The initial bracket is the extreme draws' own quantiles, and with more draws it spans far more than the mixture's own scale. The first Illinois steps from so wide a bracket often shrink it by less than half, and two such steps spend the slack before the secant's superlinear phase starts; the row then bisects. The pilot measured this: at S = 4,000 on the gaussian baseline shape, today took a median 43 passes and B1 44, within the cap but with no gain, while the calibration at S = 1,000 was 2.07 times faster. A better starting bracket would leave the slack for the superlinear phase; it is not part of this registration.", "",
        "### Check 6: cost (seconds, one run each, the two arms back to back in one R session; reported, not a gate)", "",
        "| arm | today | B1 | B1 / today |", "|---|---|---|---|",
        sprintf("| %s | %.1f | %.1f | %.2f |", cost$arm, cost$today, cost$B1, cost$ratio))
writeLines(md, file.path("studies", "findings", "quantile_b1.md"))
cat(md, sep = "\n")
