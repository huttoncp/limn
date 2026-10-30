## p10: does REML close the few-groups shortfall in new-group coverage?
##
## P3 found that with 8 groups, joint draws (M3) cover a new group's row 0.879
## (gaussian) and 0.877 (Poisson) at nominal 0.90. Maas and Hox (2005, pp.
## 89-90) found the standard errors of the group-level variances too small with
## 50 groups or fewer, in linear models fitted by restricted maximum
## likelihood. That maximum likelihood also estimates the group-level SD too
## low with so few groups is limn's own measurement, made by this script: the
## mean estimate was 0.673 (gaussian) and 0.720 (Poisson) under ML, and 0.739
## and 0.784 under REML, against a true SD of 0.8. illume's `reml = TRUE`
## (exact for gaussian, approximate otherwise) is the cheapest candidate fix,
## so this measures it on P3's own design and seeds: the ML arm should
## reproduce P3's M3 cells.
##
## Coverage is judged by the randomised PIT, as in P3. Both arms see the same
## data and the same standard-normal draws, so their difference is paired.
##
## A prototype, kept as it was run on 24 September 2026 with illume 0.0.8.9000
## (1000 replicates per family); the findings are in p10_few_groups_reml.md.
## It reads the fit's joint precision and parameter vectors by name, as that
## version laid them out, so a later illume may need changes. Run from the
## package root: Rscript studies/prototypes/p10_few_groups_reml.R [R] [cores]
args <- commandArgs(TRUE)
R     <- if (length(args) >= 1) as.integer(args[1]) else 300L
ncore <- if (length(args) >= 2) as.integer(args[2]) else 2L
res_dir <- file.path("studies", "prototypes", "results")

worker_setup <- function() {
  suppressPackageStartupMessages({ library(illume); library(Matrix) })
  TRUE
}

## ---- shared pieces, as in P3 ------------------------------------------------
mu_joint <- function(f) {
  rn <- rownames(f$jointPrecision); fx <- f$opt$par; rd <- f$sdr$par.random
  mu <- numeric(length(rn))
  for (nm in unique(rn))
    mu[rn == nm] <- if (nm %in% names(rd)) rd[names(rd) == nm] else fx[names(fx) == nm]
  mu
}
## N(0, Q^-1) from given standard-normal draws Z, holding flat directions
rmvn_prec_z <- function(Q, Z) {
  e <- eigen(as.matrix(Q), symmetric = TRUE)
  keep <- e$values > 1e-8 * max(e$values)
  list(D = e$vectors[, keep, drop = FALSE] %*% (Z[seq_len(sum(keep)), , drop = FALSE] /
                                                  sqrt(e$values[keep])),
       held = sum(!keep))
}
gh <- function(K) {
  i <- seq_len(K - 1); b <- sqrt(i / 2)
  J <- matrix(0, K, K); J[cbind(i, i + 1)] <- b; J[cbind(i + 1, i)] <- b
  e <- eigen(J, symmetric = TRUE); list(x = e$values * sqrt(2), w = e$vectors[1, ]^2)
}
pit <- function(y, fam, eta_s, sd_s = NULL, re_sd_s = NULL, G20, v) {
  if (fam == "gaussian") {
    tot <- if (is.null(re_sd_s)) sd_s else sqrt(sd_s^2 + re_sd_s^2)
    return(mean(pnorm(y, eta_s, tot)))
  }
  cdf <- function(q) {
    if (q < 0) return(0)
    if (is.null(re_sd_s)) return(mean(ppois(q, exp(eta_s))))
    mean(vapply(seq_along(eta_s), function(s)
      sum(G20$w * ppois(q, exp(eta_s[s] + re_sd_s[s] * G20$x))), 0))
  }
  lo <- cdf(y - 1); hi <- cdf(y)
  lo + v * (hi - lo)                      # v: the row's own uniform, shared by both arms
}

rep_reml <- function(rep, fam, G = 8L, n = 6L, sb = 0.8, S = 1000L) {
  ## P3's seed and data-generating steps, unchanged
  set.seed(700000 + rep + (if (fam == "poisson") 50000 else 0) + G * 1000)
  b0 <- if (fam == "gaussian") 1 else log(3); b1 <- 0.5
  d <- data.frame(g = factor(rep(seq_len(G), each = n)), x = rnorm(G * n))
  b <- rnorm(G, 0, sb)
  gen <- function(eta) if (fam == "gaussian") rnorm(length(eta), eta, 1) else rpois(length(eta), exp(eta))
  d$y <- gen(b0 + b1 * d$x + b[d$g])
  ts <- data.frame(g = factor(seq_len(G), levels = levels(d$g)), x = rnorm(G))
  ts$y <- gen(b0 + b1 * ts$x + b)
  nn <- 5L; tn <- data.frame(x = rnorm(nn)); tn$y <- gen(b0 + b1 * tn$x + rnorm(nn, 0, sb))
  ## common random numbers for both arms
  set.seed(800000 + rep + (if (fam == "poisson") 50000 else 0))
  Z <- matrix(rnorm(64L * S), 64L, S); v_seen <- runif(G); v_new <- runif(nn)
  G20 <- gh(20)
  out <- list()
  for (arm in c("ML", "REML")) {
    f <- tryCatch(suppressWarnings(suppressMessages(
           ilm_model(y ~ x + (1 | g), data = d, family = fam, joint = TRUE,
                     reml = (arm == "REML"), verbose = FALSE))),
         error = function(e) NULL)
    if (is.null(f) || is.null(f$jointPrecision)) return(NULL)
    Q <- f$jointPrecision; rn <- rownames(Q); mu <- mu_joint(f)
    ib <- which(rn == "beta"); it <- which(rn == "theta"); iu <- which(rn == "bvec")
    idp <- which(rn == "logdisp")
    jd <- rmvn_prec_z(Q, Z); D <- jd$D + mu
    sd_of <- if (length(idp)) exp(D[idp, ]) else NULL
    re_sd <- exp(D[it, ])
    u_seen <- vapply(seq_len(G), function(i) {
      eta <- as.numeric(c(1, ts$x[i]) %*% D[ib, , drop = FALSE]) + D[iu[i], ]
      pit(ts$y[i], fam, eta, sd_of, NULL, G20, v_seen[i])
    }, 0)
    u_new <- vapply(seq_len(nn), function(i) {
      eta <- as.numeric(c(1, tn$x[i]) %*% D[ib, , drop = FALSE])
      pit(tn$y[i], fam, eta, sd_of, re_sd, G20, v_new[i])
    }, 0)
    out[[arm]] <- data.frame(rep = rep, fam = fam, arm = arm,
                             target = c(rep("seen", G), rep("new", nn)),
                             u = c(u_seen, u_new), held = jd$held,
                             sb_hat = sqrt(f$Sigma$g[1, 1]))
  }
  do.call(rbind, out)
}

## ---- run --------------------------------------------------------------------
library(parallel)
cl <- makeCluster(ncore)
invisible(clusterCall(cl, worker_setup))
clusterExport(cl, c("mu_joint", "rmvn_prec_z", "gh", "pit", "rep_reml"))
t0 <- proc.time()[3]
jobs <- expand.grid(rep = seq_len(R), fam = c("gaussian", "poisson"), stringsAsFactors = FALSE)
res <- parLapplyLB(cl, seq_len(nrow(jobs)), function(j, jobs)
  tryCatch(rep_reml(jobs$rep[j], jobs$fam[j]), error = function(e) NULL), jobs = jobs)
t1 <- proc.time()[3]
stopCluster(cl)
ok <- vapply(res, Negate(is.null), TRUE)
z <- do.call(rbind, res[ok])

z$c90 <- z$u > 0.05 & z$u < 0.95; z$c50 <- z$u > 0.25 & z$u < 0.75
## replicate-level means first: rows within a replicate share a fit
agg <- aggregate(cbind(c90, c50) ~ rep + fam + arm + target, data = z, FUN = mean)
mse <- function(v) c(mean = mean(v), se = sd(v) / sqrt(length(v)))
tab <- do.call(data.frame, aggregate(cbind(c90, c50) ~ fam + arm + target, data = agg, FUN = mse))
## paired REML - ML difference in 90% coverage
w <- reshape(agg[c("rep", "fam", "target", "arm", "c90")], idvar = c("rep", "fam", "target"),
             timevar = "arm", direction = "wide")
w$diff <- w$c90.REML - w$c90.ML
dtab <- do.call(data.frame, aggregate(diff ~ fam + target, data = w, FUN = mse))
sbt <- do.call(data.frame, aggregate(sb_hat ~ fam + arm, data = unique(z[c("rep", "fam", "arm", "sb_hat")]),
                                     FUN = function(v) c(mean = mean(v), median = median(v))))

cat(sprintf("%d of %d replicates usable (both arms), %.1f min\n", sum(ok), length(ok), (t1 - t0) / 60))
cat("\n== coverage of central 90% and 50% intervals, 8 groups of 6, joint draws (M3) ==\n")
print(format(tab, digits = 3), row.names = FALSE)
cat("\n== paired REML - ML difference in 90% coverage ==\n")
print(format(dtab, digits = 3), row.names = FALSE)
cat("\n== estimated group SD (truth 0.8) ==\n")
print(format(sbt, digits = 3), row.names = FALSE)
saveRDS(list(z = z, tab = tab, diff = dtab, sb = sbt, minutes = (t1 - t0) / 60),
        file.path(res_dir, "p10_few_groups_reml.rds"))
