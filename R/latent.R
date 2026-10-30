## ---- a correlation over time, carried to new times --------------------------------
##
## A row's latent value, draw by draw, for AR(1), CAR(1) and random-walk terms.
## The process is Markov, so given the fitted cells each stretch of new times
## depends only on the cells around it:
##   - after a series' last cell, a forecast: forward from that cell;
##   - between two cells, a bridge between them;
##   - before its first cell, backward from that cell;
##   - a new series (or any under "new"), a fresh path: from the stationary
##     distribution for AR(1) and CAR(1), from zero at its first time for a walk.
## Times in one stretch are drawn in sequence, each given the one before, so a
## series' draws are one path, not independent rows. Rows naming the same
## series and time share one value, across scenarios as well.

## Each draw's latent parameters, as vectors over draws.
#' @keywords internal
#' @noRd
lmn_latent_par <- function(lat, S, step = NULL) {
  per_draw <- function(a) {
    if (is.null(a)) return(NULL)
    if (length(dim(a)) == 3L && any(dim(a)[1:2] != 1L))
      stop("forecasting a correlation over time with more than one latent value per ",
           "cell (a multinomial outcome) is not in limn yet.", call. = FALSE)
    a <- as.numeric(a)
    if (length(a) == 1L) rep(a, S) else a
  }
  list(type = lat$type, sig2 = per_draw(lat$Sigma), rho = per_draw(lat$rho),
       q = per_draw(lat$var_per_time), step = if (is.null(step)) 1 else as.numeric(step))
}

## The correlation across gaps `dt` (a vector), one column per draw: rho per grid
## step for AR(1), rho per unit of time for CAR(1).
#' @keywords internal
#' @noRd
lmn_latent_corr <- function(dt, par) {
  h <- if (identical(par$type, "ar1")) round(dt / par$step) else dt
  outer(h, par$rho, function(h, r) r^h)
}

## The distribution of a latent value given its neighbours, one row per value
## and one column per draw. `kind`:
##   "start"  - no neighbour: the stationary distribution (AR(1), CAR(1));
##   "step"   - one neighbour L at distance dt_from, before or after (a walk
##              from an anchor has L = 0);
##   "bridge" - L at dt_from before and R at dt_to after.
#' @keywords internal
#' @noRd
lmn_latent_step <- function(kind, L, R, dt_from, dt_to, par) {
  n <- length(dt_from)
  if (identical(par$type, "rw1")) {
    q <- matrix(par$q, n, length(par$q), byrow = TRUE)
    if (kind == "bridge") {
      w <- dt_from / (dt_from + dt_to)
      return(list(mean = L + w * (R - L), var = q * (dt_from * dt_to / (dt_from + dt_to))))
    }
    return(list(mean = L, var = q * dt_from))
  }
  s2 <- matrix(par$sig2, n, length(par$sig2), byrow = TRUE)
  if (kind == "start") return(list(mean = s2 * 0, var = s2))
  r1 <- lmn_latent_corr(dt_from, par)
  if (kind == "step") return(list(mean = r1 * L, var = s2 * (1 - r1^2)))
  r2 <- lmn_latent_corr(dt_to, par)
  den <- 1 - r1^2 * r2^2
  list(mean = (r1 * (1 - r2^2) * L + r2 * (1 - r1^2) * R) / den,
       var = s2 * (1 - r1^2) * (1 - r2^2) / den)
}

## The latent value of every row, one column per draw. `fresh` marks rows whose
## series is drawn afresh; the others continue their fitted cells. Returns the
## values, and which rows were drawn rather than read off a cell. With
## `noise = FALSE` every step takes its conditional mean: the cells carried
## forward with no new shocks, which for a Gaussian latent is the conditional
## mean of the new times given the cells.
#' @keywords internal
#' @noRd
lmn_latent_paths <- function(ar, cells, a_cell, par, fresh, seed, noise = TRUE) {
  n <- nrow(ar); S <- ncol(a_cell)
  tm <- as.numeric(ar$time); grp <- as.character(ar$group)
  ctime <- as.numeric(cells$time); cgrp <- as.character(cells$group)
  out <- matrix(0, n, S)
  at_cell <- !fresh & !is.na(ar$cell)
  if (any(at_cell)) out[at_cell, ] <- a_cell[ar$cell[at_cell], , drop = FALSE]
  todo <- which(!at_cell)
  if (!length(todo)) return(list(values = out, drawn = rep(FALSE, n)))

  ## one value per series and time, in stretches that depend on nothing else
  key <- paste(grp, tm, sep = "\r")
  first <- !duplicated(key[todo])
  k_row <- todo[first]                             # a row standing for each value
  K <- length(k_row)
  k_of_row <- match(key, key[k_row])               # each row's value
  kt <- tm[k_row]; kg <- grp[k_row]; kfresh <- fresh[k_row]
  prev <- ar$prev_cell[k_row]; nxt <- ar$next_cell[k_row]
  stretch <- ifelse(kfresh, paste("new", kg),
                    paste(kg, ifelse(is.na(prev), "", prev), ifelse(is.na(nxt), "", nxt)))
  backward <- !kfresh & is.na(prev)
  ord <- order(stretch, ifelse(backward, -kt, kt))
  pos <- stats::ave(seq_along(ord), stretch[ord], FUN = seq_along)[order(ord)]

  from_kind <- character(K); from_idx <- integer(K); dt_from <- numeric(K)
  to_idx <- rep(NA_integer_, K); dt_to <- rep(NA_real_, K)
  prev_key <- integer(K)                           # the value before, in its stretch
  prev_key[ord] <- c(NA, ord[-K])
  same <- c(FALSE, stretch[ord][-1] == stretch[ord][-K])
  prev_key[ord][!same] <- NA
  anchor_t <- tapply(ctime, cgrp, min)
  fresh_t <- if (any(kfresh)) tapply(kt[kfresh], kg[kfresh], min) else numeric(0)
  for (k in seq_len(K)) {
    if (pos[k] > 1L) {
      from_kind[k] <- "key"; from_idx[k] <- prev_key[k]; dt_from[k] <- abs(kt[k] - kt[prev_key[k]])
    } else if (kfresh[k]) {
      if (identical(par$type, "rw1")) {
        t0 <- min(c(if (kg[k] %in% names(anchor_t)) anchor_t[[kg[k]]], fresh_t[[kg[k]]]))
        from_kind[k] <- "anchor"; dt_from[k] <- kt[k] - t0
      } else from_kind[k] <- "start"
    } else {
      c0 <- if (backward[k]) nxt[k] else prev[k]
      from_kind[k] <- "cell"; from_idx[k] <- c0; dt_from[k] <- abs(kt[k] - ctime[c0])
    }
    if (!kfresh[k] && !is.na(prev[k]) && !is.na(nxt[k])) {
      to_idx[k] <- nxt[k]; dt_to[k] <- ctime[nxt[k]] - kt[k]
    }
  }

  Z <- if (noise) lmn_with_seed(seed, matrix(stats::rnorm(K * S), K, S)) else matrix(0, K, S)
  val <- matrix(0, K, S)
  for (p in seq_len(max(pos))) {
    for (kind in c("start", "step", "bridge")) {
      ks <- which(pos == p & switch(kind, start = from_kind == "start",
                                    step = from_kind != "start" & is.na(to_idx),
                                    bridge = !is.na(to_idx)))
      if (!length(ks)) next
      L <- matrix(0, length(ks), S)
      fk <- from_kind[ks]
      if (any(fk == "cell")) L[fk == "cell", ] <- a_cell[from_idx[ks][fk == "cell"], , drop = FALSE]
      if (any(fk == "key")) L[fk == "key", ] <- val[from_idx[ks][fk == "key"], , drop = FALSE]
      R <- if (kind == "bridge") a_cell[to_idx[ks], , drop = FALSE] else NULL
      st <- lmn_latent_step(kind, L, R, dt_from[ks], dt_to[ks], par)
      val[ks, ] <- st$mean + sqrt(pmax(st$var, 0)) * Z[ks, , drop = FALSE]
    }
  }
  out[todo, ] <- val[k_of_row[todo], , drop = FALSE]
  drawn <- rep(FALSE, n); drawn[todo] <- TRUE
  list(values = out, drawn = drawn)
}

## The latent's variance at each row, for a population mean to integrate over:
## the stationary variance for AR(1) and CAR(1); for a walk, the variance built
## up since the series' first cell (or its first time here, for a new series).
#' @keywords internal
#' @noRd
lmn_latent_var <- function(ar, cells, par) {
  n <- nrow(ar); S <- length(if (identical(par$type, "rw1")) par$q else par$sig2)
  if (!identical(par$type, "rw1")) return(matrix(par$sig2, n, S, byrow = TRUE))
  grp <- as.character(ar$group); tm <- as.numeric(ar$time)
  t0 <- tapply(as.numeric(cells$time), as.character(cells$group), min)
  first <- ifelse(grp %in% names(t0), t0[grp], stats::ave(tm, grp, FUN = min))
  outer(pmax(0, tm - first), par$q)
}
