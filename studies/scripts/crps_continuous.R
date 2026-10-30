## crps_continuous: the continuous CRPS, exact for gaussian mixtures, checked as
## registered in studies/designs/crps_continuous.md (with amendments 1 and 2).
##
## Run from the package root, with LIMN_STUDY_LIB set to a pinned library:
##   Rscript studies/scripts/crps_continuous.R [reps_slow] [reps_fast]
## reps_slow: check 6's repetitions for arms over 5 minutes; reps_fast: for the rest.
source("studies/scripts/pinned_library.R")
tools::psnice(value = 19)
Sys.setenv(OMP_NUM_THREADS = "1")
args <- commandArgs(TRUE)
reps_slow <- if (length(args) >= 1) as.integer(args[1]) else 1L
reps_fast <- if (length(args) >= 2) as.integer(args[2]) else 3L
suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
version <- as.character(utils::packageVersion("limn"))
out_dir <- file.path("studies", "runs", version, "crps_continuous")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
t_start <- Sys.time()
quiet <- function(e) suppressMessages(suppressWarnings(e))
stamp <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

## ---- families in illume's parameterisations, with location e and scale s ------
fams <- list(
  gaussian = list(p = function(q, m, s, lt) stats::pnorm(q, m, s, lower.tail = lt),
                  q = function(p, m, s) stats::qnorm(p, m, s), axis = "raw"),
  lognormal = list(p = function(q, m, s, lt) stats::plnorm(q, m, s, lower.tail = lt),
                   q = function(p, m, s) stats::qlnorm(p, m, s), axis = "log",
                   lp = function(u, m, s, lt) stats::pnorm((u - m) / s, lower.tail = lt, log.p = TRUE)),
  weibull = list(p = function(q, m, s, lt) stats::pweibull(q, 1 / s, exp(m), lower.tail = lt),
                 q = function(p, m, s) stats::qweibull(p, 1 / s, exp(m)), axis = "log",
                 ## z = (log y - e)/s: F = 1 - exp(-e^z), S = exp(-e^z)
                 lp = function(u, m, s, lt) { z <- (u - m) / s
                   if (lt) log(-expm1(-exp(z))) else -exp(z) }),
  loglogistic = list(p = function(q, m, s, lt) stats::plogis((log(q) - m) / s, lower.tail = lt),
                     q = function(p, m, s) exp(m + s * stats::qlogis(p)), axis = "log",
                     lp = function(u, m, s, lt) stats::plogis((u - m) / s, lower.tail = lt, log.p = TRUE)),
  ## beta: mean plogis(e), precision 1 / s
  beta = list(p = function(q, m, s, lt) stats::pbeta(q, stats::plogis(m) / s, (1 - stats::plogis(m)) / s, lower.tail = lt),
              q = function(p, m, s) stats::qbeta(p, stats::plogis(m) / s, (1 - stats::plogis(m)) / s), axis = "logit",
              lp = function(u, m, s, lt) stats::pbeta(stats::plogis(u), stats::plogis(m) / s, (1 - stats::plogis(m)) / s,
                                                      lower.tail = lt, log.p = TRUE)))
mkdist <- function(f) {
  F <- fams[[f]]
  list(family = f, discrete = FALSE, lp_names = "mu", par_names = "logsd",
       d = function(y, lp, par, log = FALSE) stop("unused"),
       p = function(q, lp, par, lower.tail = TRUE) F$p(q, lp$mu, exp(par$logsd), lower.tail),
       q = function(prob, lp, par) F$q(prob, lp$mu, exp(par$logsd)),
       r = function(lp, par) stop("unused"), mean = function(lp, par) stop("unused"))
}
## one row's mixture as a prediction
mkpred <- function(f, m, s, w) lmn_new_pred(mkdist(f), list(mu = matrix(m, 1)), par = list(logsd = log(s)), weights = w)

## ---- references -----------------------------------------------------------------
lse <- function(L) { mx <- max(L); if (!is.finite(mx)) return(mx); mx + log(sum(exp(L - mx))) }
## integrate() over the whole support, the integrand in log space from the
## family's own tails; NA (a failed case) on any warning or error
ref_axis <- function(f, m, s, w, y) {
  F <- fams[[f]]
  if (f == "gaussian") {
    lF <- function(z, lt) vapply(z, function(zz) lse(log(w) + stats::pnorm(zz, m, s, lower.tail = lt, log.p = TRUE)), 0)
    g1 <- function(z) exp(2 * lF(z, TRUE)); g2 <- function(z) exp(2 * lF(z, FALSE)); uy <- y; lJ <- NULL
  } else {
    lJ <- if (F$axis == "logit") function(u) stats::plogis(u, log.p = TRUE) + stats::plogis(-u, log.p = TRUE) else function(u) u
    lF <- function(u, lt) vapply(u, function(uu) lse(log(w) + F$lp(uu, m, s, lt)), 0)
    g1 <- function(u) exp(lJ(u) + 2 * lF(u, TRUE)); g2 <- function(u) exp(lJ(u) + 2 * lF(u, FALSE))
    uy <- if (F$axis == "logit") stats::qlogis(y) else log(y)
  }
  tryCatch(stats::integrate(g1, -Inf, uy, rel.tol = 1e-13, subdivisions = 20000L)$value +
             stats::integrate(g2, uy, Inf, rel.tol = 1e-13, subdivisions = 20000L)$value,
           warning = function(w) NA_real_, error = function(e) NA_real_)
}
## the raw-axis reference (condition 3)
ref_raw <- function(f, m, s, w, y) {
  F <- fams[[f]]
  Fm <- function(z, lt) vapply(z, function(zz) sum(w * F$p(zz, m, s, lt)), 0)
  top <- if (f == "beta") 1 else Inf
  tryCatch(stats::integrate(function(z) Fm(z, TRUE)^2, 0, y, rel.tol = 1e-13, subdivisions = 20000L)$value +
             stats::integrate(function(z) Fm(z, FALSE)^2, y, top, rel.tol = 1e-13, subdivisions = 20000L)$value,
           warning = function(w) NA_real_, error = function(e) NA_real_)
}
## the split reference (amendment 4): integrate() over each interval between the
## components' own quantiles (and y), -Inf to Inf, so no narrow or far mass is
## missed; up to `nb` breakpoints. NA (a failed case) on any warning or error.
ref_split <- function(f, m, s, w, y, nb = 1000L) {
  F <- fams[[f]]; ax <- F$axis
  tr <- switch(ax, raw = identity, log = log, logit = stats::qlogis)
  pr <- c(1e-14, 1e-8, 1e-4, 0.01, 0.25, 0.5, 0.75, 0.99, 1 - 1e-4, 1 - 1e-8, 1 - 1e-14)
  bk <- suppressWarnings(tr(as.vector(outer(pr, seq_along(m), function(p, k) F$q(p, m[k], s[k])))))
  bk <- sort(unique(bk[is.finite(bk)]))
  if (length(bk) > nb) bk <- unique(stats::quantile(bk, seq(0, 1, length.out = nb), names = FALSE))
  uy <- tr(y); bk <- sort(unique(c(bk, if (is.finite(uy)) uy)))
  if (ax == "raw") {
    lF <- function(z, lt) vapply(z, function(zz) lse(log(w) + stats::pnorm(zz, m, s, lower.tail = lt, log.p = TRUE)), 0)
    lJ <- function(u) 0
  } else {
    lJ <- if (ax == "logit") function(u) stats::plogis(u, log.p = TRUE) + stats::plogis(-u, log.p = TRUE) else function(u) u
    lF <- function(u, lt) vapply(u, function(uu) lse(log(w) + F$lp(uu, m, s, lt)), 0)
  }
  ends <- c(-Inf, bk, Inf)
  tryCatch({
    tot <- 0
    for (k in seq_len(length(ends) - 1L)) {
      a <- ends[k]; b <- ends[k + 1L]
      below <- if (is.finite(uy)) b <= uy else (uy > 0)
      tot <- tot + stats::integrate(function(u) exp(lJ(u) + 2 * lF(u, below)), a, b,
                                    rel.tol = 1e-13, subdivisions = 5000L)$value
    }
    tot
  }, warning = function(e) NA_real_, error = function(e) NA_real_)
}
exact_g <- function(m, s, w, y) as.numeric(lmn_crps_gaussian(mkpred("gaussian", m, s, w), y, tau = NA))
rel <- function(a, b) abs(a - b) / abs(b)

## ---- mixtures -----------------------------------------------------------------
lu <- function(n, a, b) exp(stats::runif(n, log(a), log(b)))
## a random mixture as registered; `t_max` and the scale bounds by family
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
## an outcome: the mixture's median plus k quartile spreads, on the family's axis
pick_y <- function(f, mx, k) {
  pr <- mkpred(f, mx$m, mx$s, mx$w)
  q <- function(p) lmn_pred_quantile(pr, p)
  ax <- fams[[f]]$axis
  tr <- switch(ax, raw = identity, log = log, logit = stats::qlogis)
  it <- switch(ax, raw = identity, log = exp, logit = stats::plogis)
  it(tr(q(0.5)) + k * (tr(q(0.75)) - tr(q(0.25))))
}
## the six shapes (and the added ones), on each family's axis, with S draws
shapes <- function(f, S) {
  base <- switch(f, gaussian = 1, lognormal = 0.5, weibull = 0.5, loglogistic = 0.3, beta = 0.05)
  tiny <- switch(f, gaussian = 1e-4, lognormal = 1e-4, weibull = 1e-3, loglogistic = 1e-3, beta = 1e-5)
  sep <- if (f == "gaussian") 20 else 20 * base
  u <- function(...) rep(1 / S, S)
  out <- list(
    baseline = list(m = stats::rnorm(S, 0, 0.3 * base), s = rep(base, S)),
    a_bimodal = list(m = c(rep(0, S %/% 2), rep(sep, S - S %/% 2)), s = rep(base, S)),
    b_scale_spread = list(m = stats::rnorm(S, 0, 0.5 * base),
                          s = pmin(base * 10^stats::runif(S, -1.5, 1.5), if (f == "loglogistic") 1.95 else if (f == "beta") 1.6 else if (f == "lognormal") 3 else Inf)),
    c_tail = list(m = stats::rnorm(S, 0, 0.3 * base), s = rep(base, S), k = c(10, -10, 12, -12)),
    d1_degenerate = list(m = stats::rnorm(S, 0, tiny), s = rep(tiny, S)),
    d2_spiky = list(m = stats::rnorm(S, 0, if (f == "gaussian") 1 else 0.3), s = rep(tiny, S)))
  if (f != "gaussian") {
    add <- switch(f,
      weibull = list(weibull_shape0.3 = list(m = stats::rnorm(S, 0, 0.1), s = rep(1 / 0.3, S))),
      beta = list(beta_0.3_0.3 = list(m = stats::rnorm(S, 0, 0.05), s = rep(1 / 0.6, S))),
      loglogistic = lapply(c(s0.8 = 0.8, s1.2 = 1.2, s1.6 = 1.6, s1.9 = 1.9),
                           function(sv) list(m = stats::rnorm(S, 0, 0.1), s = rep(sv, S))),
      lognormal = list(sdlog3 = list(m = stats::rnorm(S, 0, 0.1), s = rep(3, S)),
                       ## the unclipped scale spread, sdlog to about 16 (amendment 4)
                       sdlog15_spread = list(m = stats::rnorm(S, 0, 0.25), s = 0.5 * 10^stats::runif(S, -1.5, 1.5))))
    ## each added shape alone, and mixed half and half with components of another scale
    mixed <- lapply(add, function(a) list(m = a$m, s = c(a$s[seq_len(S %/% 2)], rep(base, S - S %/% 2))))
    names(mixed) <- paste0(names(add), "_mixed")
    out <- c(out, add, mixed)
  }
  lapply(out, function(x) { x$w <- rep(1 / S, S); x })
}
row_ys <- function(f, mx, nr) {
  ks <- if (!is.null(mx$k)) rep(mx$k, length.out = nr) else rep(c(-3, -1, 0, 0.5, 1, 3, 6, -6), length.out = nr)
  vapply(ks, function(k) pick_y(f, mx, k), 0)
}

res <- list()
## ---- check 1: the exact form -----------------------------------------------------
stamp("check 1")
set.seed(2850)
c1 <- do.call(rbind, lapply(seq_len(200), function(j) {
  mx <- rand_mix("gaussian"); y <- pick_y("gaussian", mx, stats::runif(1, -12, 12))
  e <- exact_g(mx$m, mx$s, mx$w, y)
  r <- scoringRules::crps_mixnorm(y, matrix(mx$m, 1), matrix(mx$s, 1), matrix(mx$w, 1))
  data.frame(S = length(mx$m), rel = rel(e, r))
}))
set.seed(2855)
sh100 <- shapes("gaussian", 100)
c1b <- do.call(rbind, lapply(names(sh100), function(nm) {
  mx <- sh100[[nm]]; ys <- row_ys("gaussian", mx, 4)
  do.call(rbind, lapply(ys, function(y) {
    rf <- ref_split("gaussian", mx$m, mx$s, mx$w, y, nb = 1e6)
    data.frame(case = nm, y = y, rel = if (is.na(rf)) NA else rel(exact_g(mx$m, mx$s, mx$w, y), rf), ref_failed = is.na(rf))
  }))
}))
res$check1 <- list(mixnorm = c1, integrate = c1b)

## ---- check 2: S = 100 (above) and S = 4,000 against scoringRules ------------------
stamp("check 2")
set.seed(2856)
sh4k <- shapes("gaussian", 4000)
c2 <- do.call(rbind, lapply(names(sh4k), function(nm) {
  mx <- sh4k[[nm]]; ys <- row_ys("gaussian", mx, 4)
  do.call(rbind, lapply(ys, function(y) {
    r <- scoringRules::crps_mixnorm(y, matrix(mx$m, 1), matrix(mx$s, 1), matrix(mx$w, 1))
    data.frame(case = nm, y = y, rel = rel(exact_g(mx$m, mx$s, mx$w, y), r))
  }))
}))
res$check2 <- c2

## ---- check 3: tau, calibrated then held out ----------------------------------------
stamp("check 3")
rrow <- function(m, s, w, y) {
  pr <- mkpred("gaussian", m, s, w)
  lo <- min(stats::qnorm(1e-10, m, s)); hi <- max(stats::qnorm(1 - 1e-10, m, s))
  width <- max(max(min(y, hi), lo) - lo, hi - min(max(y, lo), hi)) / 16
  list(r = width / min(s), err = rel(as.numeric(lmn_crps_uniform(pr, y)), exact_g(m, s, w, y)))
}
set.seed(2851)
grid <- expand.grid(sep = c(0, 5, 20, 50), ratio = c(1, 10, 100, 1000), k = c(0, 5, 10), S = c(100, 1000))
cal <- do.call(rbind, lapply(seq_len(nrow(grid)), function(g) {
  G <- grid[g, ]
  do.call(rbind, lapply(1:4, function(j) {
    S <- G$S
    m <- stats::rnorm(S, 0, 0.3) + c(rep(0, S %/% 2), rep(G$sep, S - S %/% 2))
    s <- exp(stats::runif(S, -log(G$ratio) / 2, log(G$ratio) / 2))
    w <- rep(1 / S, S); y <- stats::median(m) + G$k * sign(stats::runif(1) - 0.5)
    v <- rrow(m, s, w, y)
    data.frame(G, r = v$r, err = v$err)
  }))
}))
taus <- 2^-(0:6)
pass <- vapply(taus, function(t) all(cal$err[cal$r <= t] <= 1e-13), TRUE)
tau <- if (any(pass)) max(taus[pass]) else NA_real_
## held out
set.seed(2852)
ho <- c(lapply(shapes("gaussian", 1000), function(mx) list(mx = mx, ys = row_ys("gaussian", mx, 12))),
        lapply(seq_len(200), function(j) { mx <- rand_mix("gaussian"); list(mx = mx, ys = pick_y("gaussian", mx, stats::runif(1, -12, 12))) }))
c3 <- do.call(rbind, lapply(seq_along(ho), function(h) {
  mx <- ho[[h]]$mx
  do.call(rbind, lapply(ho[[h]]$ys, function(y) {
    pr <- mkpred("gaussian", mx$m, mx$s, mx$w)
    v <- lmn_crps_gaussian(pr, y, tau = tau)
    data.frame(case = if (h <= 6) names(ho)[h] else "random", fast = attr(v, "fast") > 0,
               rel = rel(as.numeric(v), exact_g(mx$m, mx$s, mx$w, y)))
  }))
}))
res$check3 <- list(calibration = cal, tau = tau, held_out = c3)

## ---- check 4: the other families ---------------------------------------------------
stamp("check 4")
raw_variant <- function(f, m, s, w, y, tau_c) {
  ## today's raw axis with tau_c panels sized by the narrowest raw-axis component
  pr <- mkpred(f, m, s, w); F <- fams[[f]]
  lo <- min(F$q(1e-10, m, s)); hi <- max(F$q(1 - 1e-10, m, s))
  sc <- min((F$q(0.75, m, s) - F$q(0.25, m, s))) / 1.349
  np <- min(4096L, max(16L, ceiling((hi - lo) / (tau_c * sc))))
  as.numeric(lmn_crps_uniform(pr, y, panels = np))
}
fam_others <- c("lognormal", "weibull", "loglogistic", "beta")
c4 <- list()
for (f in fam_others) {
  stamp(" ", f, "calibration")
  set.seed(2853)
  gcal <- expand.grid(sep = c(0, 2, 5), ratio = c(1, 10, 100, 1000), k = c(-3, 0, 3), S = c(100, 1000))
  rows <- lapply(seq_len(nrow(gcal)), function(g) {
    G <- gcal[g, ]; S <- G$S
    b <- switch(f, lognormal = 0.5, weibull = 0.5, loglogistic = 0.3, beta = 0.05)
    lim <- switch(f, lognormal = c(0.05, 3), weibull = c(0.1, 3.4), loglogistic = c(0.05, 1.95), beta = c(1e-3, 1.6))
    m <- stats::rnorm(S, 0, 0.3 * b) + c(rep(0, S %/% 2), rep(G$sep * b, S - S %/% 2))
    s <- pmin(pmax(b * exp(stats::runif(S, -log(G$ratio) / 2, log(G$ratio) / 2)), lim[1]), lim[2])
    mx <- list(m = m, s = s, w = rep(1 / S, S)); y <- pick_y(f, mx, G$k)
    list(mx = mx, y = y, ref = ref_split(f, m, s, mx$w, y))
  })
  ## the shapes check 4 holds out join the calibration, at S = 100 (amendment 4)
  sh_c <- shapes(f, 100)
  sh_c <- sh_c[!names(sh_c) %in% c("baseline", "a_bimodal", "b_scale_spread", "c_tail", "d1_degenerate", "d2_spiky")]
  for (nm in names(sh_c)) {
    mx <- sh_c[[nm]]
    for (y in row_ys(f, mx, 3)) rows[[length(rows) + 1L]] <- list(mx = mx, y = y, ref = ref_split(f, mx$m, mx$s, mx$w, y))
  }
  ok <- !vapply(rows, function(r) is.na(r$ref), TRUE)
  tau_c <- NA_real_
  for (t in 2^-(0:6)) {
    errs <- vapply(rows[ok], function(r) rel(as.numeric(quiet(lmn_crps_axis(mkpred(f, r$mx$m, r$mx$s, r$mx$w), r$y, tau_c = t))), r$ref), 0)
    if (all(errs <= 1e-11, na.rm = FALSE)) { tau_c <- t; break }
  }
  stamp(" ", f, "tau_c", tau_c)
  tc_use <- if (is.na(tau_c)) 2^-6 else tau_c
  ## held out: the shapes at S = 100 and 1,000, and 200 random mixtures
  nshape <- 2L * length(shapes(f, 10))         # before the seed, so it draws nothing from it
  set.seed(2858)
  hos_s <- unlist(lapply(c(100, 1000), function(S) { sh <- shapes(f, S); names(sh) <- paste0(names(sh), "_S", S)
                                                    lapply(sh, function(mx) list(mx = mx, ys = row_ys(f, mx, 4))) }), recursive = FALSE)
  set.seed(2859)
  hos <- c(hos_s, lapply(seq_len(200), function(j) { mx <- rand_mix(f); list(mx = mx, ys = pick_y(f, mx, stats::runif(1, -12, 12))) }))
  names(hos)[names(hos) == ""] <- "random"
  rr <- do.call(rbind, lapply(seq_along(hos), function(h) {
    mx <- hos[[h]]$mx
    do.call(rbind, lapply(hos[[h]]$ys, function(y) {
      rf <- ref_split(f, mx$m, mx$s, mx$w, y)
      v <- quiet(lmn_crps_axis(mkpred(f, mx$m, mx$s, mx$w), y, tau_c = tc_use))
      raw <- if (h <= nshape) raw_variant(f, mx$m, mx$s, mx$w, y, tc_use) else NA_real_
      data.frame(family = f, case = names(hos)[h], y = y, ref_failed = is.na(rf),
                 value = as.numeric(v), rel = if (is.na(rf)) NA else rel(as.numeric(v), rf),
                 method_failed = !is.null(attr(v, "failed")), integrated = attr(v, "integrated"),
                 tail_bound = attr(v, "tail_bound"), raw_rel = if (is.na(rf)) NA else rel(raw, rf))
    }))
  }))
  ## the raw-axis reference cross-check: 20 rows where the raw axis converges
  set.seed(2857)
  rawx <- do.call(rbind, lapply(seq_len(20), function(j) {
    repeat { mx <- rand_mix(f); if (f != "loglogistic" || max(mx$s) <= 1.2) break }
    y <- pick_y(f, mx, stats::runif(1, -3, 3))
    rf <- ref_raw(f, mx$m, mx$s, mx$w, y)
    v <- quiet(lmn_crps_axis(mkpred(f, mx$m, mx$s, mx$w), y, tau_c = tc_use))
    data.frame(family = f, ref_failed = is.na(rf), rel = if (is.na(rf)) NA else rel(as.numeric(v), rf))
  }))
  c4[[f]] <- list(tau_c = tau_c, held_out = rr, raw_reference = rawx)
}
## the infinite case
inf_rows <- vapply(c(2, 2.5), function(sv) {
  v <- quiet(lmn_crps_continuous(mkpred("loglogistic", c(0, 0.5), c(0.5, sv), c(0.5, 0.5)), 1))
  identical(as.numeric(v), Inf)
}, TRUE)
res$check4 <- c4; res$check4_inf <- inf_rows

## ---- check 6: cost at the profiling sizes --------------------------------------------
stamp("check 6")
set.seed(1)
G <- 50; n <- 2000
d <- data.frame(g = factor(sample(sprintf("g%02d", 1:G), n, TRUE)), x = rnorm(n))
d$yg <- 1 + 0.5 * d$x + rnorm(G, 0, 0.5)[as.integer(d$g)] + rnorm(n)
fg <- quiet(illume::ilm_model(yg ~ x + (1 | g), data = d, verbose = FALSE))
pg <- quiet(ilm_predict(fg, d, target = "response", groups = "fitted", nsim = 1000))
S <- 40; Tn <- 48; h <- 12
ts <- data.frame(site = factor(rep(sprintf("s%02d", 1:S), each = Tn)), t = rep(1:Tn, S)); ts$x <- rnorm(nrow(ts))
ts$y <- log1p(rpois(nrow(ts), exp(1 + 0.2 * ts$x + ave(rnorm(nrow(ts), 0, 0.15), ts$site, FUN = cumsum))))
ftg <- quiet(illume::ilm_model(y ~ x, data = ts, ar = illume::ilm_rw1(~ t | site), verbose = FALSE))
fut <- quiet(ilm_future(ts, time = "t", group = "site", h = h, response = "y")); fut$x <- 0
fcg <- quiet(ilm_forecast(ftg, fut, nsim = 1000))
tt <- rweibull(n, 1.5, exp(1 + 0.4 * d$x)); cc <- rexp(n, 0.15)
d$time <- pmin(tt, cc); d$event <- as.integer(tt <= cc)
fw <- quiet(illume::ilm_model(time ~ x, data = d, family = "weibull", censor = illume::ilm_surv(d$time, d$event), verbose = FALSE))
pw <- quiet(ilm_predict(fw, d, target = "response", nsim = 1000))
timeit <- function(expr, reps) {
  e <- substitute(expr); env <- parent.frame()
  stats::median(vapply(seq_len(reps), function(k) system.time(quiet(eval(e, env)))[["elapsed"]], 0))
}
yfc <- rep(1, nrow(fcg$rows))
tc_w <- if (is.na(c4$weibull$tau_c)) 2^-6 else c4$weibull$tau_c
## the slow arms once (amendment 3), each rerun twice more, taking the median,
## when its ratio to today lands between 1.6x and 2.4x
one <- function(expr) { e <- substitute(expr); env <- parent.frame()
  function() system.time(quiet(eval(e, env)))[["elapsed"]] }
arms <- list(g_today = one(lmn_crps_uniform(pg, d$yg)), g_exact = one(lmn_crps_gaussian(pg, d$yg, tau = NA)),
             g_hybrid = one(lmn_crps_gaussian(pg, d$yg, tau = tau)),
             w_today = one(lmn_crps_uniform(pw, tt)), w_axis = one(lmn_crps_axis(pw, tt, tau_c = tc_w)))
runs <- lapply(arms, function(f) vapply(seq_len(reps_slow), function(k) f(), 0))
guard <- function(arm, base) {
  ratio <- stats::median(runs[[arm]]) / stats::median(runs[[base]])
  if (ratio >= 1.6 && ratio <= 2.4 && length(runs[[arm]]) < 3L)
    runs[[arm]] <<- c(runs[[arm]], arms[[arm]](), arms[[arm]]())
}
for (a2 in c("g_exact", "g_hybrid")) guard(a2, "g_today")
guard("w_axis", "w_today")
med <- vapply(runs, stats::median, 0)
cost <- data.frame(
  arm = c("gaussian today", "gaussian exact", "gaussian hybrid", "forecast today", "forecast exact", "forecast hybrid",
          "weibull today", "weibull axis"),
  runs = c(lengths(runs)[1:3], reps_fast, reps_fast, reps_fast, lengths(runs)[4:5]),
  seconds = c(med[["g_today"]], med[["g_exact"]], med[["g_hybrid"]],
              timeit(lmn_crps_uniform(fcg, yfc), reps_fast), timeit(lmn_crps_gaussian(fcg, yfc, tau = NA), reps_fast),
              timeit(lmn_crps_gaussian(fcg, yfc, tau = tau), reps_fast),
              med[["w_today"]], med[["w_axis"]]))
cost$ratio_to_today <- cost$seconds / cost$seconds[c(1, 1, 1, 4, 4, 4, 7, 7)]
fast_share <- attr(lmn_crps_gaussian(pg, d$yg, tau = tau), "fast") / n
res$check6 <- list(cost = cost, fast_share = fast_share, reps = c(slow = reps_slow, fast = reps_fast))
minutes <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
saveRDS(res, file.path(out_dir, "results.rds"))

## ---- findings ---------------------------------------------------------------------------
fmt <- function(v) if (all(is.na(v))) "NA" else sprintf("%.1e", max(v, na.rm = TRUE))
ho4 <- do.call(rbind, lapply(c4, `[[`, "held_out"))
rx4 <- do.call(rbind, lapply(c4, `[[`, "raw_reference"))
rule <- function(ok) if (isTRUE(ok)) "yes" else "NO"
checks <- data.frame(
  check = c("1 exact vs crps_mixnorm, 200 random (max rel)",
            "1 the split reference vs the exact form, 24 gaussian rows (reference validation, max rel)",
            "1 reference failures", "2 exact vs crps_mixnorm, six cases S = 4,000 (max rel)",
            "3 tau", "3 held out: fast-path rows (max rel vs exact)",
            paste("4", fam_others, "tau_c"), "4 held out, all families (max rel vs integrate)",
            "4 reference failures", "4 method failures (NA rows)",
            "4 sdlog15_spread, S = 100 and 1,000 (max rel)",
            "4 raw-axis reference cross-check (max rel)", "4 raw-axis reference failures", "4 Inf for s = 2 and 2.5"),
  value = c(fmt(c1$rel), fmt(c1b$rel), sum(c1b$ref_failed), fmt(c2$rel),
            format(tau), fmt(c3$rel[c3$fast]),
            vapply(fam_others, function(f) format(c4[[f]]$tau_c), ""), fmt(ho4$rel),
            sum(ho4$ref_failed), sum(ho4$method_failed), fmt(ho4$rel[grepl("^sdlog15_spread_S", ho4$case)]), fmt(rx4$rel), sum(rx4$ref_failed), paste(inf_rows, collapse = ", ")),
  rule = c("<= 1e-12", "<= 1e-12", "0", "<= 1e-12", "reported", "<= 1e-12",
           rep("reported", 4), "<= 1e-10", "0", "reported", "<= 1e-10", "<= 1e-10", "0", "TRUE, TRUE"),
  pass = c(rule(max(c1$rel) <= 1e-12), rule(max(c1b$rel, na.rm = TRUE) <= 1e-12), rule(sum(c1b$ref_failed) == 0),
           rule(max(c2$rel) <= 1e-12), "", rule(!any(c3$fast) || max(c3$rel[c3$fast]) <= 1e-12),
           rep("", 4), rule(max(ho4$rel, na.rm = TRUE) <= 1e-10), rule(sum(ho4$ref_failed) == 0), "",
           rule(max(ho4$rel[grepl("^sdlog15_spread_S", ho4$case)]) <= 1e-10),
           rule(max(rx4$rel, na.rm = TRUE) <= 1e-10), rule(sum(rx4$ref_failed) == 0), rule(all(inf_rows))))
write.csv(checks, file.path(out_dir, "checks.csv"), row.names = FALSE)
by_case <- stats::aggregate(cbind(rel, raw_rel) ~ family + case, data = transform(ho4, case = sub("_S[0-9]+$", "", case)),
                            FUN = function(v) max(v, na.rm = TRUE), na.action = stats::na.pass)
write.csv(by_case, file.path(out_dir, "by_case.csv"), row.names = FALSE)
write.csv(cost, file.path(out_dir, "cost.csv"), row.names = FALSE)
md <- c("# The continuous CRPS, exact for gaussian mixtures: checks", "",
        sprintf("## limn %s (%s; illume %s; scoringRules %s)", version, format(Sys.Date()),
                utils::packageVersion("illume"), utils::packageVersion("scoringRules")), "",
        "Generated by `studies/scripts/crps_continuous.R`, to the design in `studies/designs/crps_continuous.md` (with its amendments). Do not edit by hand.", "",
        lmn_study_header(), "",
        sprintf("Run on limn %s: %.1f minutes on 1 core. Check 6 repetitions: %d for the 2,000-row arms, %d for the rest.",
                system("git rev-parse --short HEAD", intern = TRUE), minutes, reps_slow, reps_fast), "",
        "## Run 4, the registered run after amendment 4", "",
        "| check | result | rule | pass |", "|---|---|---|---|",
        sprintf("| %s | %s | %s | %s |", checks$check, checks$value, checks$rule, checks$pass), "",
        "### Check 4 by family and shape: worst relative error, the method and the raw-axis variant", "",
        "| family | shape | method | raw axis |", "|---|---|---|---|",
        sprintf("| %s | %s | %.1e | %.1e |", by_case$family, by_case$case, by_case$rel, by_case$raw_rel), "",
        "### Check 6: cost (seconds, median of the repetitions)", "",
        "| arm | runs | seconds | ratio to today |", "|---|---|---|---|",
        sprintf("| %s | %d | %.1f | %.2f |", cost$arm, cost$runs, cost$seconds, cost$ratio_to_today), "",
        sprintf("Share of the 2,000 profiling rows on the gaussian fast path: %.2f.", fast_share), "",
        "## The four runs, and why there were four", "",
        "Nothing earlier is replaced: run 3's findings are kept as generated in `studies/findings/crps_continuous_run3.md`, with its tables in `studies/runs/0.0.0.9000/crps_continuous_run3/`, and each amendment is in the design file, written before the run it governs.", "",
        "- **Run 1** (5d6ee62) stopped in check 3: the script's `row_ys()` indexed eight outcome offsets with twelve rows, giving NA outcomes, and the gaussian fast-path indicator failed on NA. Fixed in the script, and the CRPS now returns NA for an NA outcome (151d5df). No results.",
        "- **Run 2** (151d5df) stopped in check 4's held-out beta rows: R's `qbeta()` returns quantiles off by up to 0.5 in probability at tiny shape parameters, so the method's range and panel scale were wrong. Fixed in the method: a row trusts its quantiles only when F at lo and S at hi say so, and otherwise is integrated over the whole axis (ae5900b). Checks 1 to 3 and the calibration had run.",
        "- **Run 3** (42d4f4d) completed and FAILED four rules, as registered: check 1's integrate() comparison (Inf on four near-degenerate rows, two reference failures), and check 4's held-out rule (max Inf) with two reference failures and one method NA. Weibull shape 0.3 was 1.2e-10 against 1e-10. The diagnosis (amendment 4) found the registered reference silently wrong on most failing rows, and three method faults: a lognormal with sdlog 14 to 16 was 4.1% and 20% low (the upper tail cut at the 1 - 1e-16 quantile), one spiky lognormal row was NA (the whole-axis fallback failed), and two rows were marginal at tau_c = 1.",
        "- **Run 4** (this run) follows amendment 4: the split reference, validated against the exact form on the gaussian rows; the method's upper tail for lognormal and weibull integrated until negligible; split integration as the method's fallback; the added shapes in the calibration and fresh held-out seeds; and the unclipped sdlog-15 case kept as a named shape that must pass.")
writeLines(md, file.path("studies", "findings", "crps_continuous.md"))
cat(md, sep = "\n")
