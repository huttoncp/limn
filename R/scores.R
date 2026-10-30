## ---- scoring a predictive distribution ---------------------------------------
##
## Every score is a loss, so lower is better, and every score is computed from
## the predictive mixture itself: its density, CDF, quantiles or mean. None
## simulates a response.
##
## The workers behind ilm_scores.ilm_pred() (R/ilm_scores.R), limn's method for
## illume's generic.

## Rules and what each needs.
##   log       -log f(y), the negative log predictive density or mass
##   crps      continuous ranked probability score; for integer outcomes the
##             sum over the support, which for ordered categories is the
##             (unnormalised) ranked probability score
##   brier     (P(Y = 1) - y)^2, binary outcomes
##   pinball   quantile loss at each of `probs`
##   interval  the interval score for the central `level` interval
##             (Gneiting and Raftery 2007, section 6.2, equation (43))
##   squared   (y - predictive mean)^2
##   absolute  |y - predictive median|, the point forecast absolute error suits
##   mase      |y - median| / q, q the series' in-sample mean absolute naive
##             error at lag `period` (Hyndman & Koehler 2006; fpp3 5.8): the
##             mean over rows is the MASE
##   rmsse     (y - mean)^2 / q2, q2 the mean squared naive error: the square
##             root of the mean over rows is the RMSSE
lmn_rules <- c("log", "crps", "brier", "pinball", "interval", "squared", "absolute",
               "mase", "rmsse")

#' @keywords internal
#' @noRd
lmn_pred_scores <- function(pred, y, rules = c("log", "crps"),
                            probs = c(0.1, 0.5, 0.9), level = 0.9, naive = NULL) {
  lmn_check_response_pred(pred)
  n <- lmn_pred_n(pred)
  y <- lmn_check_outcome(pred, y, n)
  bad <- setdiff(rules, lmn_rules)
  if (length(bad))
    stop("unknown rule(s): ", paste(bad, collapse = ", "), ". The rules are ",
         paste(lmn_rules, collapse = ", "), ".", call. = FALSE)
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1)
    stop("`level` must be one number between 0 and 1.", call. = FALSE)
  if (!is.numeric(probs) || any(probs <= 0 | probs >= 1))
    stop("`probs` must be probabilities strictly between 0 and 1.", call. = FALSE)

  scaled <- intersect(rules, c("mase", "rmsse"))
  if (length(scaled)) sc <- lmn_row_scales(pred, naive)
  ok <- !is.na(y)
  out <- list()
  sub <- lmn_subset_pred(pred, which(ok))
  yo <- y[ok]
  put <- function(v) { full <- rep(NA_real_, n); full[ok] <- v; full }
  crps_bound <- rep(0, n)                                  # the discrete CRPS's error estimate
  crps_bracket <- rep(0, n)                                # and its rigorous bracket width
  ## an average is scored through its total: the log score on the total's own
  ## values (plus the density's change of scale for a continuous outcome), the
  ## CRPS scaled back
  avg_sc <- if (is.null(sub$scale)) rep(1, length(yo)) else sub$scale
  tot <- sub; tot$scale <- NULL; yt <- yo / avg_sc
  if (sub$dist$discrete) yt <- round(yt)                   # the total is a whole number
  for (r in rules) {
    if (r == "pinball") {
      for (p in probs) {
        q <- lmn_pred_quantile(sub, p)
        out[[paste0("pinball_", format(p))]] <- put((yo - q) * (p - (yo < q)))
      }
      next
    }
    out[[if (r == "interval") paste0("interval_", format(level)) else r]] <-
      put(switch(r,
        log = if (identical(sub$kind, "sample")) lmn_no_log_score() else
          lmn_score_log(tot, yt) + (if (sub$dist$discrete) 0 else log(avg_sc)),
        crps = { v <- if (identical(sub$kind, "sample")) lmn_crps_sample(tot, yt)
                   else if (sub$dist$discrete) lmn_crps_discrete(tot, yt)
                   else lmn_crps_continuous(tot, yt)
                 if (!is.null(attr(v, "error_estimate"))) {
                   crps_bound <- put(avg_sc * attr(v, "error_estimate"))
                   crps_bracket <- put(avg_sc * attr(v, "bound"))
                 }
                 avg_sc * as.numeric(v) },
        brier = lmn_score_brier(sub, yo),
        interval = lmn_score_interval(sub, yo, level),
        squared = (yo - lmn_pred_mean(sub))^2,
        absolute = abs(yo - lmn_pred_quantile(sub, 0.5)),
        mase = abs(yo - lmn_pred_quantile(sub, 0.5)) / sc$abs[ok],
        rmsse = (yo - lmn_pred_mean(sub))^2 / sc$sq[ok]))
  }
  res <- cbind(pred$rows, as.data.frame(out, check.names = FALSE))
  structure(res, class = c("ilm_pred_scores", "data.frame"),
            rules = names(out), orientation = "lower is better",
            crps_bound = if (any(crps_bound > 0, na.rm = TRUE)) crps_bound,
            crps_bracket = if (any(crps_bracket > 0, na.rm = TRUE)) crps_bracket)
}

#' @keywords internal
#' @noRd
lmn_check_outcome <- function(pred, y, n) {
  if (!is.numeric(y) && !is.logical(y))
    stop("`y` must be the observed outcomes, as numbers.", call. = FALSE)
  y <- as.numeric(y)
  if (length(y) != n)
    stop("`y` has ", length(y), " values for ", n, " predicted rows.", call. = FALSE)
  yt <- if (is.null(pred$scale)) y else y / pred$scale    # an average's total
  if (pred$dist$discrete && any(!is.na(yt) & abs(yt - round(yt)) > 1e-8))
    stop(if (is.null(pred$scale)) "the outcome is discrete, so `y` must hold whole numbers."
         else "the outcome is an average of whole numbers, so `y` times the number of rows averaged must be whole.",
         call. = FALSE)
  y
}

## The scale of the scaled errors: each series' in-sample naive-forecast error,
## at lag `period` (1 for the naive forecast, the season's length for the
## seasonal naive), from the history the model was trained on.
#' @keywords internal
#' @noRd
lmn_naive_scale <- function(history, time, group = NULL, response, period = 1L) {
  tm <- lmn_check_column(history, time, "time")
  y <- as.numeric(lmn_check_column(history, response, "response"))
  g <- if (is.null(group)) rep("1", nrow(history)) else
    as.character(lmn_check_column(history, group, "group"))
  period <- lmn_check_count(period, "period")
  keys <- unique(g)
  res <- lapply(keys, function(k) {
    ix <- which(g == k); yy <- y[ix][order(tm[ix])]
    if (length(yy) <= period)
      stop("series '", k, "' has ", length(yy), " observation(s), too few for a ",
           "naive error at lag ", period, ".", call. = FALSE)
    e <- diff(yy, lag = period)
    c(abs = mean(abs(e)), sq = mean(e^2))
  })
  out <- data.frame(key = keys, abs = vapply(res, `[[`, 0, "abs"),
                    sq = vapply(res, `[[`, 0, "sq"))
  if (any(out$abs == 0))
    stop("series ", paste0("'", out$key[out$abs == 0], "'", collapse = ", "),
         " never change in the history, so a scaled error is undefined.", call. = FALSE)
  structure(out, group = group, period = period)
}

## Each predicted row's scale, matched through its series.
#' @keywords internal
#' @noRd
lmn_row_scales <- function(pred, naive) {
  if (is.null(naive) || !is.data.frame(naive) || !all(c("key", "abs", "sq") %in% names(naive)))
    stop("the scaled errors (mase, rmsse) need `naive`, the in-sample naive ",
         "scale of each series.", call. = FALSE)
  group <- attr(naive, "group")
  key <- if (is.null(group)) rep("1", nrow(pred$rows)) else {
    if (!group %in% names(pred$rows))
      stop("the predicted rows have no '", group, "' column to match the series by.",
           call. = FALSE)
    as.character(pred$rows[[group]])
  }
  i <- match(key, naive$key)
  if (anyNA(i))
    stop("no in-sample scale for series ", paste0("'", unique(key[is.na(i)]), "'", collapse = ", "),
         ": a scaled error needs the series' own history.", call. = FALSE)
  list(abs = naive$abs[i], sq = naive$sq[i])
}

## The same predictive distribution restricted to rows i.
#' @keywords internal
#' @noRd
lmn_subset_pred <- function(pred, i) {
  n <- lmn_pred_n(pred)
  if (length(i) == n) return(pred)
  if (!is.null(pred$scale)) pred$scale <- pred$scale[i]
  if (lmn_is_empirical(pred)) {
    pred$draws <- pred$draws[i, , drop = FALSE]
    pred$estimate <- pred$estimate[i]
  } else {
    pred$lp <- lapply(pred$lp, function(M) M[i, , drop = FALSE])
    pred$par <- lapply(pred$par, function(v) if (is.matrix(v)) v[i, , drop = FALSE] else v)
  }
  pred$rows <- pred$rows[i, , drop = FALSE]
  if (!is.null(pred$origin)) pred$origin <- pred$origin[i]
  if (!is.null(pred$new_unit)) pred$new_unit <- pred$new_unit[i]
  pred
}

## Scores, PITs and calibration compare outcomes with a predictive distribution
## of the RESPONSE; a prediction of the mean has no outcome to be compared with.
#' @keywords internal
#' @noRd
lmn_check_response_pred <- function(pred) {
  if (!inherits(pred, "ilm_pred"))
    stop("`pred` must be a predictive distribution (class \"ilm_pred\").",
         call. = FALSE)
  if (lmn_is_draws(pred))
    stop("this is a prediction of the ", if (pred$kind == "mean") "mean" else
           "linear predictor", ", not of an observation, so there is no outcome to ",
         "score it against. Predict with target = \"response\".", call. = FALSE)
  invisible(pred)
}

## A distribution held as samples has no density to take the log of.
#' @keywords internal
#' @noRd
lmn_no_log_score <- function() {
  stop("the log score needs the predictive density, and a distribution held as ",
       "samples (such as a total summed from simulated outcomes) has none. Score it ",
       "with crps, pinball or interval, which the samples give directly.", call. = FALSE)
}

## CRPS from samples, by the fair (unbiased) estimator: the mean distance from
## the outcome, less the mean distance between two different samples over two
## (Ferro 2014; Zamo and Naveau 2018). Over sorted samples, the sum over pairs
## of |x_i - x_j| is 2 * sum_k (2k - M - 1) x_(k).
#' @keywords internal
#' @noRd
lmn_crps_sample <- function(pred, y) {
  X <- pred$draws; M <- ncol(X)
  vapply(seq_len(nrow(X)), function(i) {
    x <- sort(X[i, ])
    mean(abs(x - y[i])) - sum((2 * seq_len(M) - M - 1) * x) / (M * (M - 1))
  }, 0)
}

## -log of the mixture density, on the log scale throughout: log-sum-exp over
## the draws, so a density too small for double precision still scores finitely.
#' @keywords internal
#' @noRd
lmn_score_log <- function(pred, y) {
  n <- length(y); S <- length(pred$weights)
  out <- numeric(n)
  lw <- log(pred$weights)
  for (i in lmn_row_chunks(n, S)) {
    ld <- matrix(NA_real_, length(i), S)
    for (ss in lmn_draw_blocks(pred, length(i))) {
      a <- if (length(pred$pvec)) lmn_draw_args(pred, i, ss) else lmn_block_args(pred, i, ss)
      ld[, ss] <- pred$dist$d(rep(y[i], length(ss)), a$lp, a$par, log = TRUE) +
        rep(lw[ss], each = length(i))
    }
    mx <- apply(ld, 1L, max)
    fin <- is.finite(mx)
    out[i] <- Inf
    out[i][fin] <- -(mx[fin] + log(rowSums(exp(ld[fin, , drop = FALSE] - mx[fin]))))
  }
  out
}

## CRPS for a continuous outcome: the integral of (F(z) - 1{z >= y})^2 over z
## (studies/designs/crps_continuous.md).
##   gaussian     the exact form of a gaussian mixture (Grimit et al. 2006),
##                with the uniform quadrature kept only for rows whose panels
##                are fine against the narrowest component (r <= tau);
##   lognormal, weibull, loglogistic
##                quadrature on log z, beta on logit z, with panels sized by the
##                narrowest component and the upper tail from the survival
##                function, never from 1 - F;
##   any other    an error: there is no fallback to a method known to be wrong
##                (the uniform quadrature was up to 100% off on positive
##                families), and a test holds illume's families to these.
## Rows over the panel cap are integrated by integrate(); a failure there is NA,
## named in a message and counted in attribute "failed", never a number.
#' @keywords internal
#' @noRd
lmn_crps_continuous <- function(pred, y, tau = lmn_crps_tau, tau_c = lmn_crps_tau_c, cap = 4096L) {
  fam <- pred$dist$family
  if (!isTRUE(fam %in% lmn_crps_families))
    stop("limn has no CRPS for a continuous '", if (is.null(fam)) "unnamed" else fam,
         "' predictive distribution. Score it with rules = \"log\", or with \"pinball\" ",
         "or \"interval\", which need only its quantiles.", call. = FALSE)
  out <- if (identical(fam, "gaussian")) lmn_crps_gaussian(pred, y, tau)
    else lmn_crps_axis(pred, y, if (length(tau_c) > 1L) tau_c[[fam]] else tau_c, cap)
  failed <- attr(out, "failed")
  if (length(failed))
    message("the CRPS could not be computed for ", length(failed), " row(s) (",
            paste(utils::head(failed, 10L), collapse = ", "), if (length(failed) > 10L) ", ..." else "",
            "): integrate() failed there, so they are NA.")
  out
}

## The continuous families with a CRPS method, and how illume's other families
## are scored (test-crps-continuous.R holds every family of illume's
## ilm_family() to one of these, so one added upstream breaks a test first).
lmn_crps_families <- c("gaussian", "lognormal", "weibull", "loglogistic", "beta")
lmn_discrete_families <- c("binomial", "poisson", "nbinom", "multinomial", "ordinal",
                           "ordinal_probit", "ordinal_cloglog")
lmn_no_dist_families <- c("rp", "rp_odds", "rp_normal")

## The calibrated thresholds (studies/findings/crps_continuous.md, run 4): tau
## for the gaussian fast path, and tau_c, the panel width against the
## narrowest component, by family. NA for tau would mean no fast path.
lmn_crps_tau <- 0.5
## the fallback's breakpoints at most (between the components' own quantiles)
lmn_crps_breaks <- 120L
lmn_crps_tau_c <- c(lognormal = 1, weibull = 0.5, loglogistic = 1, beta = 1)

## E|N(mu, v)|; v = 0 gives |mu|
#' @keywords internal
#' @noRd
lmn_eabs_norm <- function(mu, v) {
  s <- sqrt(v)
  out <- s * 2 * stats::dnorm(mu / s) + mu * (2 * stats::pnorm(mu / s) - 1)
  z <- v == 0
  out[z] <- abs(mu[z])
  out
}

## The gaussian mixture's CRPS, exact:
##   sum_s w_s A(m_s - y, s_s^2) - sum_{j<k} w_j w_k A(m_j - m_k, s_j^2 + s_k^2)
##   - (1/2) sum_j w_j^2 A(0, 2 s_j^2),
## each component's mean and sd read from its own quantiles, so totals and
## scaled predictions are covered. Rows whose uniform panels are fine against
## the narrowest component (r <= tau) keep the quadrature.
#' @keywords internal
#' @noRd
lmn_crps_gaussian <- function(pred, y, tau = NA_real_, eps = 1e-10, panels = 16L) {
  n <- length(y); w <- pred$weights; S <- length(w)
  out <- numeric(n)
  fast <- logical(n)
  for (i in lmn_row_chunks(n, 4L)) {
    m <- lmn_draw_quantiles(pred, i, 0.5)
    sd <- lmn_draw_quantiles(pred, i, stats::pnorm(1)) - m
    if (!is.na(tau)) {
      lo <- apply(lmn_draw_quantiles(pred, i, eps), 1L, min)
      hi <- apply(lmn_draw_quantiles(pred, i, 1 - eps), 1L, max)
      yi <- y[i]
      width <- pmax(pmax(pmin(yi, hi), lo) - lo, hi - pmin(pmax(yi, lo), hi)) / panels
      ratio <- width / apply(sd, 1L, min)
      fast[i] <- !is.na(ratio) & ratio <= tau
    }
    for (k in seq_along(i)) {
      r <- i[k]
      if (fast[r]) next
      if (is.na(y[r])) { out[r] <- NA_real_; next }
      mk <- m[k, ]; v <- sd[k, ]^2
      t1 <- sum(w * lmn_eabs_norm(mk - y[r], v))
      ## the pairs j < k once each: for each j, the components after it as one
      ## contiguous vector (no index vectors, which cost six times as much)
      t2 <- 0
      for (j in seq_len(S - 1L)) {
        k <- (j + 1L):S
        t2 <- t2 + w[j] * sum(w[k] * lmn_eabs_norm(mk[j] - mk[k], v[j] + v[k]))
      }
      t2 <- t2 + 0.5 * sum(w^2 * lmn_eabs_norm(0, 2 * v))
      out[r] <- t1 - t2
    }
  }
  if (any(fast)) out[fast] <- lmn_crps_uniform(lmn_subset_pred(pred, which(fast)), y[fast])
  attr(out, "fast") <- sum(fast)
  out
}

## Quadrature on the family's natural axis: u = log z for the positive
## families, u = logit z for beta, with dz = e^u du or z(1 - z) du. On each side
## of y, panels of 8 Gauss-Legendre nodes no wider than tau_c times the
## narrowest component's scale there (its IQR on the axis over 1.349). Below y
## the integrand is F^2, above it the survival function squared, from the
## family's upper tail. A loglogistic is integrated past its 1 - 1e-10
## quantile until a bound on the rest (from 1 - F_k <= (z / c_k)^(-1/s_k), by
## convexity over the mixture) is 1e-13 of the integral; any component with
## s >= 2 makes the CRPS infinite. Rows over `cap` panels a side are
## integrated by integrate().
#' @keywords internal
#' @noRd
lmn_crps_axis <- function(pred, y, tau_c, cap = 4096L, eps = 1e-10) {
  fam <- pred$dist$family
  logit <- identical(fam, "beta"); llog <- identical(fam, "loglogistic")
  tr <- if (logit) stats::qlogis else log
  ## log of the Jacobian dz/du at u
  logJ <- if (logit) function(u) stats::plogis(u, log.p = TRUE) + stats::plogis(-u, log.p = TRUE)
    else function(u) u
  inv <- if (logit) stats::plogis else exp
  w <- pred$weights; lw <- log(w)
  gl <- lmn_gauss_legendre(8L)
  n <- length(y); out <- numeric(n); bound <- numeric(n)
  failed <- integer(0); integrated <- 0L; infinite <- integer(0)
  for (r in seq_len(n)) {
    yr <- y[r]
    if (is.na(yr)) { out[r] <- NA_real_; next }
    q <- function(p) as.vector(lmn_draw_quantiles(pred, r, p))
    q25 <- tr(q(0.25)); q75 <- tr(q(0.75))
    if (llog) {
      e <- tr(q(0.5)); s <- (q75 - q25) / (2 * stats::qlogis(0.75))
      if (any(s[w > 0] >= 2)) { out[r] <- Inf; infinite <- c(infinite, r); next }
      ## log F or log S of the mixture at u, in log space, from each
      ## component's (e, s) read off its own quantiles
      lmix <- function(u, upper) {
        L <- outer(u, e, "-") / rep(s, each = length(u))
        L <- stats::plogis(if (upper) -L else L, log.p = TRUE) + rep(lw, each = length(u))
        mx <- apply(L, 1L, max); mx[!is.finite(mx)] <- 0
        mx + log(rowSums(exp(L - mx)))
      }
    } else {
      lmix <- function(u, upper) {
        z <- inv(u)
        log(lmn_mix(pred, r, matrix(z, 1L), if (upper) "s" else "p")[1, ])
      }
    }
    q_top <- q(if (llog || logit) 1 - eps else 1 - 1e-16)
    lo_z <- min(q(eps)); hi_z <- max(q_top)
    if (logit) { lo_z <- max(lo_z, .Machine$double.xmin); hi_z <- min(hi_z, 1 - .Machine$double.eps) }
    lo <- tr(lo_z); hi <- tr(hi_z)
    uy <- if (yr <= 0) -Inf else if (logit && yr >= 1) Inf else tr(yr)
    sc <- min((q75 - q25)[w > 0]) / 1.349
    ## the quantiles behind the range and the panels must be right: a family's
    ## own quantile function can fail (R's qbeta at shape parameters of 1e4 and
    ## more, for one). A row where they are not goes to integrate() over the
    ## whole axis, which needs no quantiles.
    trusted <- is.finite(sc) && sc > 0 && is.finite(lo) && is.finite(hi) && hi > lo &&
      isTRUE(exp(lmix(lo, FALSE)) <= 1e-6) && isTRUE(exp(lmix(hi, TRUE)) <= 1e-6)
    ## the integral of exp(logJ(u) + 2 lmix(u)) over [a, b] by panels, or NA when
    ## the panels would exceed the cap
    piece <- function(a, b, upper, scale) {
      if (!(b > a)) return(0)
      np <- max(1L, ceiling((b - a) / (tau_c * scale)))
      if (!is.finite(np) || np > cap) return(NA_real_)
      h <- (b - a) / np
      mid <- a + (seq_len(np) - 0.5) * h
      u <- as.vector(outer(gl$x * h / 2, mid, "+")); wt <- rep(gl$w * h / 2, np)
      sum(wt * exp(logJ(u) + 2 * lmix(u, upper)))
    }
    v <- if (trusted) piece(lo, max(min(uy, hi), lo), FALSE, sc) + piece(min(max(uy, lo), hi), hi, TRUE, sc)
      else NA_real_
    if (llog && !is.na(v)) {
      ## past hi, to the point where the rest's bound is 1e-13 of the integral
      lb <- function(u) {
        k <- (2 / s) * e + (1 - 2 / s) * u - log(2 / s - 1) + lw
        mx <- max(k); mx + log(sum(exp(k - mx)))
      }
      target <- log(1e-13 * max(v, .Machine$double.xmin))
      top <- hi; step <- 1
      while (lb(top) > target && top < hi + 1e5) { top <- top + step; step <- step * 2 }
      bot <- hi
      for (it in 1:60) { mid <- (top + bot) / 2; if (lb(mid) > target) bot <- mid else top <- mid }
      ## the tail's panels go by the components that matter there
      Sk <- stats::plogis(-(hi - e) / s, log.p = TRUE) + lw
      keep <- Sk >= max(Sk) - log(1e30)
      v <- v + piece(max(hi, uy), top, TRUE, min((q75 - q25)[keep]) / 1.349)
      bound[r] <- exp(lb(top))
    } else if (!logit && !is.na(v)) {
      ## lognormal, weibull: past hi (the 1 - 1e-16 quantile), until the
      ## integrand on the log axis is below 1e-16 of the integral and falling.
      ## A wide lognormal's e^u S(u)^2 can peak beyond hi (amendment 4).
      g <- function(u) logJ(u) + 2 * lmix(u, TRUE)
      sck <- (q75 - q25) / 1.349
      keep <- (hi - tr(q_top)) / sck <= 10              # the components that reach hi
      sc_t <- suppressWarnings(min(sck[keep & w > 0]))
      if (!is.finite(sc_t) || sc_t <= 0) sc_t <- sc
      target <- log(1e-16 * max(v, .Machine$double.xmin))
      top <- max(hi, uy); step <- sc_t
      done <- function(u) {
        gu <- g(u)
        if (isTRUE(gu == -Inf)) return(TRUE)
        isTRUE(is.finite(gu) && gu < target && g(u + sc_t) <= gu)
      }
      while (!done(top) && top < hi + 1e5) { top <- top + step; step <- step * 2 }
      v <- v + piece(max(hi, uy), top, TRUE, sc_t)
    }
    if (is.na(v)) {
      ## over the cap, or quantiles not to be trusted: integrate() over the whole
      ## axis, split at the components' own quantiles so no narrow or far mass is
      ## missed (amendment 4). Any warning or error is NA, never a number.
      integrated <- integrated + 1L
      g1 <- function(u) exp(logJ(u) + 2 * lmix(u, FALSE))
      g2 <- function(u) exp(logJ(u) + 2 * lmix(u, TRUE))
      pr <- c(1e-14, 1e-8, 1e-4, 0.01, 0.25, 0.5, 0.75, 0.99, 1 - 1e-4, 1 - 1e-8, 1 - 1e-14)
      bk <- suppressWarnings(tr(as.vector(lmn_draw_quantiles(pred, r, pr[1]))))
      for (pp in pr[-1]) bk <- c(bk, suppressWarnings(tr(as.vector(lmn_draw_quantiles(pred, r, pp)))))
      bk <- sort(unique(bk[is.finite(bk)]))
      if (length(bk) > lmn_crps_breaks) bk <- unique(stats::quantile(bk, seq(0, 1, length.out = lmn_crps_breaks), names = FALSE))
      bk <- sort(unique(c(bk, if (is.finite(uy)) uy)))
      ends <- c(-Inf, bk, Inf)
      v <- tryCatch({
        tot <- 0
        for (k in seq_len(length(ends) - 1L)) {
          a <- ends[k]; b <- ends[k + 1L]
          f <- if (b <= uy) g1 else g2
          tot <- tot + stats::integrate(f, a, b, rel.tol = 1e-12, subdivisions = 5000L)$value
        }
        tot
      }, warning = function(w) NA_real_, error = function(e) NA_real_)
      lo_z <- 0; hi_z <- if (logit) 1 else Inf
      if (is.na(v)) failed <- c(failed, r)
    }
    ## outside [lo, hi] the indicator's stretches, exactly, on the raw axis
    out[r] <- v + max(0, lo_z - yr) + max(0, yr - hi_z)
  }
  if (length(infinite))
    message("the CRPS is infinite for ", length(infinite), " row(s): a loglogistic component ",
            "with scale s >= 2 has a tail too heavy for (1 - F)^2 to be integrable.")
  out[failed] <- NA_real_
  attr(out, "failed") <- if (length(failed)) failed
  attr(out, "integrated") <- integrated
  attr(out, "tail_bound") <- bound
  out
}

## The uniform quadrature, as before 1 October: the integral of
## (F(z) - 1{z >= y})^2 over z.
## Outside [lo, hi], the smallest and largest of the draws' 1e-10 and 1 - 1e-10
## quantiles, F is 0 or 1 to that precision, so those stretches are added
## exactly. Inside, the integral is split at y, where the indicator jumps, and
## each smooth piece is integrated by composite Gauss-Legendre: 16 panels of 8
## nodes.
#' @keywords internal
#' @noRd
lmn_crps_uniform <- function(pred, y, eps = 1e-10, panels = 16L, nodes = 8L) {
  n <- length(y)
  gl <- lmn_gauss_legendre(nodes)
  L <- panels * nodes
  out <- numeric(n)
  for (i in lmn_row_chunks(n, 2L * L)) {
    lo <- apply(lmn_draw_quantiles(pred, i, eps), 1L, min)
    hi <- apply(lmn_draw_quantiles(pred, i, 1 - eps), 1L, max)
    yi <- y[i]
    a1 <- lo; b1 <- pmax(pmin(yi, hi), lo)           # below y: F^2
    a2 <- pmin(pmax(yi, lo), hi); b2 <- hi           # above y: (1 - F)^2
    g1 <- lmn_gl_grid(a1, b1, panels, gl)
    g2 <- lmn_gl_grid(a2, b2, panels, gl)
    F1 <- lmn_mix(pred, i, g1$x, "p")
    F2 <- lmn_mix(pred, i, g2$x, "p")
    out[i] <- rowSums(g1$w * F1^2) + rowSums(g2$w * (1 - F2)^2) +
      pmax(0, lo - yi) + pmax(0, yi - hi)
  }
  out
}

## Composite Gauss-Legendre points and weights on [a, b], row by row: an
## m x (panels * nodes) matrix of each.
#' @keywords internal
#' @noRd
lmn_gl_grid <- function(a, b, panels, gl) {
  m <- length(a); M <- length(gl$x)
  h <- (b - a) / panels
  x <- matrix(0, m, panels * M); w <- x
  for (j in seq_len(panels)) {
    mid <- a + (j - 0.5) * h
    cols <- (j - 1L) * M + seq_len(M)
    x[, cols] <- mid + outer(h / 2, gl$x)
    w[, cols] <- outer(h / 2, gl$w)
  }
  list(x = x, w = w)
}

## CRPS for an integer outcome: the sum over the support of
## (F(k) - 1{k >= y})^2.
#' @keywords internal
#' @noRd
## The support runs from the mixture's own eps quantile to its 1 - eps
## quantile (eps = 1e-6), found by bisection inside the draws' quantiles.
## Beyond it the mixture's F is below eps or above 1 - eps, so the terms left
## out, F(k)^2 or (1 - F(k))^2, are each at most eps^2, and the CRPS is exact to
## rounding. An outcome outside the support is different: the points between
## the support's end and the outcome are counted as whole ones, which is off by
## up to 2 x eps each. Those rows take the union of the draws' 1e-12 and
## 1 - 1e-12 quantiles (eps_outside) instead, as every row did before;
## attribute "recomputed" counts them. Up to `max_support` values (2,000) a
## row is summed exactly. A wider row is summed by
## lmn_crps_wide(): exact near y, Euler-Maclaurin with Simpson elsewhere, to an
## error estimate of at most 1e-6 of the CRPS (studies/designs/crps_wide.md).
## The attributes give each row's "error_estimate", its rigorous "bound" (the
## monotone bracket's width, loose, for the record), "flagged" where the
## estimate exceeds 1e-6 of the CRPS, and the final "step"; all 0 (or FALSE)
## for rows summed exactly.
lmn_crps_discrete <- function(pred, y, eps = 1e-6, eps_outside = 1e-12, max_support = 2000) {
  n <- length(y)
  out <- numeric(n); bound <- numeric(n); err_est <- numeric(n); flagged <- logical(n)
  step <- rep(1, n); recomputed <- 0L
  for (i in lmn_row_chunks(n)) {
    qlo <- lmn_draw_quantiles(pred, i, eps); qhi <- lmn_draw_quantiles(pred, i, 1 - eps)
    klo <- lmn_mix_qint(pred, i, eps, apply(qlo, 1L, min), apply(qlo, 1L, max))
    khi <- lmn_mix_qint(pred, i, 1 - eps, apply(qhi, 1L, min), apply(qhi, 1L, max))
    ## outcomes outside the support: those rows' bounds as before
    far <- which(y[i] < klo | y[i] > khi)
    if (length(far) && eps_outside < eps) {
      klo[far] <- apply(lmn_draw_quantiles(pred, i[far], eps_outside), 1L, min)
      khi[far] <- apply(lmn_draw_quantiles(pred, i[far], 1 - eps_outside), 1L, max)
      recomputed <- recomputed + length(far)
    }
    wide <- (khi - klo + 1) > max_support
    if (any(wide)) {
      for (r in which(wide)) {
        b <- lmn_crps_wide(pred, i[r], y[i[r]], klo[r], khi[r])
        out[i[r]] <- b$estimate; bound[i[r]] <- b$bound; err_est[i[r]] <- b$error_estimate
        flagged[i[r]] <- b$flagged; step[i[r]] <- b$step
      }
    }
    exact <- which(!wide)
    if (!length(exact)) next
    ## rows of like support widths together, so that a narrow row is not padded
    ## to the widest's; each group in pieces small enough to evaluate at once
    width <- khi[exact] - klo[exact] + 1
    for (grp in lmn_width_groups(width)) {
      g <- exact[grp]
      L <- max(width[grp])
      for (jj in lmn_row_chunks(length(g), L)) {
        j <- g[jj]; ii <- i[j]
        Fk <- lmn_mix_cdf_run(pred, ii, klo[j], L)
        K <- klo[j] + matrix(0:(L - 1L), length(j), L, byrow = TRUE)
        inside <- K <= khi[j]
        yi <- y[ii]
        out[ii] <- rowSums(inside * (Fk - (K >= yi))^2) +
          pmax(0, klo[j] - yi) + pmax(0, yi - khi[j] - 1)
      }
    }
  }
  structure(out, error_estimate = err_est, bound = bound, flagged = flagged, step = step,
            recomputed = recomputed)
}

## The mixture's CDF at klo, klo + 1, ..., klo + L - 1 for rows i: its value at
## klo, then the mass function summed along the row. The mixture is linear in
## its draws, so this is exact, and a mass function costs far less than a
## distribution function (an incomplete beta or gamma) for large counts.
#' @keywords internal
#' @noRd
lmn_mix_cdf_run <- function(pred, i, klo, L) {
  F0 <- lmn_mix(pred, i, matrix(klo), "p")[, 1L]
  if (L == 1L) return(matrix(F0))
  K <- klo + matrix(seq_len(L - 1L), length(i), L - 1L, byrow = TRUE)
  Fk <- cbind(F0, lmn_mix(pred, i, K, "d"), deparse.level = 0)
  for (k in 2:L) Fk[, k] <- Fk[, k] + Fk[, k - 1L]  # running sums along each row
  pmin(Fk, 1)
}

## Groups of rows whose widths are within a factor of 1.5 of each other.
#' @keywords internal
#' @noRd
lmn_width_groups <- function(width) {
  o <- order(width)
  grp <- integer(length(width)); g <- 1L; lo <- width[o[1L]]
  for (k in o) {
    if (width[k] > 1.5 * lo + 5) { g <- g + 1L; lo <- width[k] }
    grp[k] <- g
  }
  split(seq_along(width), grp)
}

## One row's discrete CRPS over a support wider than `max_support` values
## (studies/designs/crps_wide.md).
## The integers within `h` of y, where 1{k >= y} jumps, are summed exactly.
## The rest is two smooth monotone pieces, F(k)^2 below and (1 - F(k))^2
## above, each summed by Euler-Maclaurin: composite Simpson on integer nodes
## a, a + d, ..., e, plus (g(a) + g(e)) / 2 and (G'(e) - G'(a)) / 12 from
## central differences, with the integers past e summed exactly. The error
## estimate is Richardson's, |S_d - S_2d| / 15 per piece, times `safety`. The
## step halves, from a power of 2 near width / 400, until the estimate is at
## most `rel_tol` of the CRPS, or the pieces are summed exactly. The rigorous
## bracket of a monotone piece's block sums, d times its rise, is returned for
## the record.
#' @keywords internal
#' @noRd
lmn_crps_wide <- function(pred, i, y, klo, khi, h = 50, rel_tol = 1e-6, safety = 10) {
  cache <- new.env(parent = emptyenv())
  Fat <- function(k) {                              # the mixture's CDF at integers, cached
    key <- as.character(k)
    have <- vapply(key, exists, TRUE, envir = cache, inherits = FALSE)
    if (any(!have)) {
      kn <- k[!have]
      v <- lmn_mix(pred, i, matrix(kn, 1L), "p")[1, ]
      for (j in seq_along(kn)) assign(as.character(kn[j]), v[j], envir = cache)
    }
    vapply(key, get, 0, envir = cache, inherits = FALSE)
  }
  exact_run <- function(a, b, g) {                  # sum of g(F(k)) over a..b by the pmf run
    if (b < a) return(0)
    sum(g(lmn_mix_cdf_run(pred, i, a, b - a + 1)[1, ]))
  }
  gb <- function(F) F^2; ga <- function(F) (1 - F)^2
  near_lo <- max(klo, y - h); near_hi <- min(khi, y + h)
  near <- if (near_hi >= near_lo) {
    Fk <- lmn_mix_cdf_run(pred, i, near_lo, near_hi - near_lo + 1)[1, ]
    k <- near_lo:near_hi
    sum((Fk - (k >= y))^2)
  } else 0
  pieces <- list(list(a = klo, b = min(khi, y - h - 1), g = gb),
                 list(a = max(klo, y + h + 1), b = khi, g = ga))
  simpson <- function(v, d) d / 3 * (v[1] + v[length(v)] + 4 * sum(v[seq(2, length(v) - 1, 2)]) +
                                        2 * sum(v[seq(3, length(v) - 2, 2)]))
  piece_sum <- function(p, d) {
    n <- p$b - p$a
    if (n < 0) return(list(sum = 0, err = 0, bracket = 0))
    m <- if (d >= 2) floor(n / (4 * d)) * 4 else 0
    if (m < 8) return(list(sum = exact_run(p$a, p$b, p$g), err = 0, bracket = 0))
    e <- p$a + m * d
    nodes <- p$a + d * (0:m)
    gv <- p$g(Fat(nodes))
    s_d <- simpson(gv, d); s_2d <- simpson(gv[seq(1, m + 1, 2)], 2 * d)
    gd <- p$g(Fat(c(p$a - 1, p$a + 1, e - 1, e + 1)))
    em <- (gv[1] + gv[m + 1]) / 2 + ((gd[4] - gd[3]) / 2 - (gd[2] - gd[1]) / 2) / 12
    list(sum = s_d + em + exact_run(e + 1, p$b, p$g), err = abs(s_d - s_2d) / 15,
         bracket = d * abs(gv[m + 1] - gv[1]))
  }
  width <- khi - klo + 1
  d <- 2^ceiling(log2(max(2, width / 400)))
  outside <- max(0, klo - y) + max(0, y - khi - 1)
  repeat {
    r <- lapply(pieces, piece_sum, d = d)
    total <- near + outside + sum(vapply(r, `[[`, 0, "sum"))
    est <- safety * sum(vapply(r, `[[`, 0, "err"))
    if (est <= rel_tol * total || d < 2) break
    d <- d / 2
  }
  list(estimate = total, error_estimate = est, bound = sum(vapply(r, `[[`, 0, "bracket")),
       flagged = est > rel_tol * total, step = d)
}

## Brier score for a binary outcome: (P(Y = 1) - y)^2.
#' @keywords internal
#' @noRd
lmn_score_brier <- function(pred, y) {
  if (!pred$dist$discrete || any(!y %in% c(0, 1)))
    stop("the Brier score is for a binary outcome coded 0/1.", call. = FALSE)
  F0 <- lmn_pred_cdf(pred, rep(0, length(y)))
  F1 <- lmn_pred_cdf(pred, rep(1, length(y)))
  if (any(abs(F1 - 1) > 1e-10))
    stop("the Brier score is for a binary outcome, but the predictive ",
         "distribution puts mass above 1.", call. = FALSE)
  (1 - F0 - y)^2
}

## The interval score for the central `level` interval [l, u] (Gneiting and
## Raftery 2007, section 6.2, equation (43)): its width, plus 2 / alpha times
## any distance by which y falls outside it.
#' @keywords internal
#' @noRd
lmn_score_interval <- function(pred, y, level) {
  alpha <- 1 - level
  l <- lmn_pred_quantile(pred, alpha / 2)
  u <- lmn_pred_quantile(pred, 1 - alpha / 2)
  (u - l) + (2 / alpha) * (l - y) * (y < l) + (2 / alpha) * (y - u) * (y > u)
}
