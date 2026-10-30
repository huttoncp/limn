# Design: the continuous CRPS, exact for gaussian mixtures (crps_continuous)

Written before the method is built and before any validation run. The
direction was settled on 1 October 2026.

**Why:** an edge-case study
measured today's `lmn_crps_continuous()` against the exact gaussian-mixture
CRPS, S = 1,000 draws, 12 rows per case:

| mixture | today's worst relative error |
|---|---|
| (b) component sds from 0.01 to 10 across draws | 7.7e-3 |
| (d2) sd 1e-4, means spread by 1 (a staircase CDF) | 9.8e-4 |
| (a) bimodal, half the means 20 sd away | 4.0e-13 |
| baseline, tail (c), tight degenerate (d1) | at most 4.7e-11 |

Today the integral runs over [lo, hi], the extreme 1e-10 and 1 - 1e-10
quantiles of the draws, split only at y. Each piece gets 16 equal panels of
8-node Gauss-Legendre. The panel width is set by the widest component, while
the integrand's structure is as narrow as the smallest sd, so (b) and (d2)
are under-resolved. Even 64 panels leaves 1.6e-3 on (b).

Families that reach this code (illume's `ilm_dist()`, discrete = FALSE):
gaussian, beta, weibull, lognormal and loglogistic. Gaussian totals from
`ilm_aggregate()` reach it too. Samples have their own CRPS (the fair form),
and nothing else in limn integrates with uniform panels: quantiles are found by
bisection on F, the log score and the PIT are exact, and the discrete CRPS has
its own error estimate.

## The method

**Gaussian mixtures: the exact closed form** (Grimit et al. 2006; as in
`scoringRules::crps_mixnorm`). With weights w_s, means m_s and sds s_s,

  CRPS = sum_s w_s A(m_s - y, s_s^2)
         - sum_{j<k} w_j w_k A(m_j - m_k, s_j^2 + s_k^2)
         - (1/2) sum_j w_j^2 A(0, 2 s_j^2),

where A(mu, v) = E|N(mu, v)| = sqrt(v) 2 phi(mu / sqrt(v)) + mu (2 Phi(mu / sqrt(v)) - 1).

The pairs are used once by symmetry, S(S - 1)/2 of them. Rows and pairs are
chunked so memory stays below about 200 MB.
The draws' weights are used as they are, so covariate-draw blocks and
`ilm_aggregate()` totals are covered.

**Gaussian fast path, only where it is proven.** Today's 16 x 8 quadrature is
kept for a row when the indicator

  r = (largest panel width) / (smallest component sd)

is at most a threshold tau. Otherwise the row takes the exact form. tau is set
by the calibration in check 3 and confirmed on held-out cases. If no tau keeps
every calibration case at or below 1e-13, there is no fast path: exact only.

**Beta, weibull, lognormal, loglogistic: quadrature with resolved panels, on
the family's natural axis.** These have no closed form for a mixture's pair
term. (Amended 1 October, before any run: the axes, the tails and the infinite
case below.)
- **The axis.** The positive families are integrated on u = log x, with
  dx = e^u du. They are multiplicative, and a raw-axis scale misreads a skewed
  component: at sdlog 3 a lognormal's IQR/1.349 is 5.5, while F climbs from 0
  to 0.25 inside [0, 0.13]. Beta is integrated on u = logit x, with
  dx = x(1 - x) du, which also tames a density that is infinite at 0 or 1.
- **The panels.** On that axis, the pieces either side of y get the smallest
  panel count that makes the panel width at most tau_c times the smallest
  component scale, capped at 4,096 panels a piece. The scale is the
  component's IQR on the axis over 1.349, from its own quantiles: sdlog for a
  lognormal, 2.197 s / 1.349 for a loglogistic, and so on, for any family.
  A row over the cap is integrated by `stats::integrate()` at rel.tol 1e-10,
  over the whole support, and counted in attribute "integrated".
- **The lower stretch**, from the support's lower end to lo (each component's
  1e-10 quantile), has F at most 1e-10. So it adds at most 1e-20 times its
  length on the raw axis, which is dropped, and the bound is kept.
- **The upper stretch beyond hi.** Today it is dropped. For a loglogistic it
  is not small. In illume's form, z = (log y - e) / s, the upper tail is
  1 - F ~ (x / e^e)^(-1/s), so the stretch beyond hi = e^e 10^(10 s) is
  about 1e-20 hi / (2/s - 1).
  - Measured against integration to Inf (one component, e = 0, y = 1), the formula matches the loss at every s:
    1.6e-10 relative at s = 1, 1.6e-8 at 1.2, 1.5e-5 at 1.5, 6.6e-4 at 1.67,
    1.2e-2 at 1.8, 17% at 1.92 and 33% at 1.95.
  - So the upper piece runs on the log axis past hi, to the point X where
    the tail bound
      sum_k w_k c_k^(2/s_k) X^(1 - 2/s_k) / (2/s_k - 1),
    from 1 - F <= sum_k w_k (x / c_k)^(-1/s_k), is at most 1e-13 times the
    integral so far. On the log axis X stays manageable: log X is about 570
    at s = 1.9 and 1,170 at s = 1.95 (one component, e = 0, y = 1). The bound is kept as the row's error estimate.
  - For weibull and lognormal the upper piece runs to the 1 - 1e-16 quantile,
    past which (1 - F)^2 is below 1e-32.
- **An infinite CRPS.** A loglogistic component with s >= 2 makes the CRPS
  infinite, since the integral of (1 - F)^2 diverges. Such a row returns Inf,
  with one message per call naming the rows. It is never a finite number.
- tau_c is set by calibration, as tau is (check 4).

## Validation (pre-registered)

The script is `studies/scripts/crps_continuous.R`, run on the pinned library
mrp-442639f (scoringRules 1.1.3), 1 core.

**Random mixtures** (checks 1, 3 and 4), fixed in advance by their seeds.
Each mixture is drawn as follows:
- S: log-uniform on [2, 1000], rounded.
- The weights: independent Exp(1), normalised, so unequal.
- Gaussian: means m_s = mu0 + t * N(0, 1), with mu0 ~ N(0, 1) and t log-uniform
  on [0.01, 30]. sds s_s = b * exp(U(-a, a)), with b log-uniform on [1e-3, 10]
  and a uniform on [0, log(1000) / 2], so the spread is up to 1000-fold within
  a mixture.
- The outcome: the mixture's median plus k times its quartile spread, with k
  uniform on [-12, 12].
- The other families, in illume's forms:
  - lognormal: e_s as for the gaussian means, with t up to 3, and sdlog log-uniform on [0.05, 3];
  - weibull: z = (log y - e)/s Gumbel, s log-uniform on [0.1, 3.4], so shape 0.3 to 10;
  - loglogistic: s log-uniform on [0.05, 1.95];
  - beta: mean plogis(e) and precision phi log-uniform on [0.6, 1000], which covers beta(0.3, 0.3).
  - Within a mixture their scale parameters spread as the gaussian sds do, within the same bounds.

1. **The exact form.**
   - Against `scoringRules::crps_mixnorm` on 200 random gaussian mixtures
     (seed 2850). Rule: at most 1e-12 relative.
   - Against `integrate()` at rel.tol 1e-13 over the whole line on the six
     cases (baseline, a, b, c, d1, d2) at S = 100. Rule: at most 1e-10
     relative, which is integrate's own limit on d1.
2. **S = 100 and S = 4,000** on the six cases: the exact form against
   `integrate()` (S = 100) and against scoringRules (S = 4,000). Same rules.
3. **Calibrating tau, then holding it out.**
   - Calibration set: gaussian mixtures over a grid of bimodal separation
     {0, 5, 20, 50} sd, sd ratio across draws {1, 10, 100, 1000}, the outcome at
     {0, 5, 10} sd, S in {100, 1000}, 4 rows each, seed 2851. For each row,
     record r and the quadrature's relative error against the exact form.
   - tau is the largest r such that every calibration row with r <= tau has
     error at most 1e-13, ten times inside the target.
   - Held out: the six cases plus 200 random mixtures as above, seed 2852. Rule: every
     row that takes the fast path is at most 1e-12 relative from the exact form.
4. **The other families.**
   - Beta, weibull, lognormal and loglogistic mixtures, each in the six shapes
     carried to that family's axis: a spread of locations (bimodal), a
     1000-fold spread of scale parameters, a tail outcome, near-degenerate, and
     spiky. S = 100 and 1,000.
   - Added shapes: densities infinite at a bound (weibull shape 0.3, beta(0.3,
     0.3)), and heavy tails (loglogistic s = 0.8, 1.2, 1.6 and 1.9, lognormal
     sdlog 3), each alone and mixed with components of other scales.
   - tau_c is calibrated on a grid like check 3's on these families (seed
     2853), then held out on the shapes above plus 200 random mixtures per
     family (seed 2854).
   - The reference is `integrate()` at rel.tol 1e-13 over the WHOLE support,
     (0, Inf) or (0, 1), on the same axis transform, so truncation error is
     measured, not hidden.
   - Rule: at most 1e-10 relative on every held-out row. The "integrated"
     count and the largest tail bound are reported.
   - A loglogistic component with s = 2 and with s = 2.5: the row returns Inf,
     with the message. Rule: Inf, never finite.
   - A raw-axis variant (today's axis with tau_c panels) is run beside it on
     the same shapes, so the log and logit axes are shown to matter or not.

5. **limn's existing CRPS tests** (against scoringRules on gaussian
   predictions) pass unchanged, with today's tolerances.
6. **Cost**, reported against a 2x rule. At the profiling sizes (2,000 rows x
   1,000 draws gaussian; the 480-row gaussian forecast; a weibull AFT
   prediction of 2,000 x 1,000), median of 3 on one core:
   - today's code;
   - exact only;
   - the hybrid.

   If exact only is more than 2x slower than today on the profiling shape,
   and the hybrid does not bring it back within 2x, it goes to the maintainer before
   the merge.

## Amendment 2 (1 October, before any run; the conditions set on review)

1. **The upper tail from the survival function, never 1 - F.**
   - Each component's survival S_k(x) comes from its family's own upper tail:
     illume's exported `p(q, lp, par, lower.tail = FALSE)`, so it is accurate
     down to about 1e-308.
   - They are mixed as S_mix = sum_k w_k S_k, and the log axis's integrand is
     exp(u + 2 log S_mix(u)). So e^u is never formed past u = 709, and 1 - F
     is never formed at all.
   - illume's p() has no log.p argument, so log S_k is not available directly.
     A true log-sum-exp is therefore not possible through the exported API.
     Over the range integrated, S_k stays above 1e-308 (at s = 1.95 and
     u = 1,170, S is about 1e-261), so the linear mixture then the log is
     exact to rounding there.
   - A component whose S_k underflows to 0 contributes 0, which is correct to
     within e^u 1e-616.
   - The same holds for the lower piece (F from the lower tail), for beta's
     1 - F near 1 on the logit axis, and for the reference integrands.
2. **Failure is never a silent number.**
   - When integrate() fails or warns, whether as the cap fallback in the
     method or as a reference in the study, the method returns NA for that
     row, with one message per call naming the rows, counted in attribute
     "failed".
   - The study reports a reference failure as a FAILED case, not a skipped
     one.
3. **Raw-axis reference rows as a cross-check.**
   - 20 rows per family where the raw axis converges, which is everything but
     loglogistic with s > 1.2.
   - The reference integrate() runs on the raw axis there, to the same 1e-10
     rule.

**Calibration rules, made explicit.**
- tau: every gaussian calibration row at or below it is within 1e-13 of the
  exact form.
- tau_c: every calibration row of the other families at or below it is within
  1e-11 of the reference, ten times inside check 4's 1e-10 rule.
- Both are searched over halvings, from 1 down to 1/64. The largest value
  that passes is taken. If none passes, the result is reported as is: no
  fast path for tau, the cap fallback for tau_c.

**Timing.** The largest pieces are piloted first: check 6 at 2,000 x 1,000,
and S = 4,000. If the pilot puts the whole run above twice the 1 core-hour
estimate, it stops and reports.

## Amendment 3 (1 October, after the timing pilot, before the run)

**The pilot**, 1 core:
- The exact gaussian form, with its pairs built as index vectors, was 5.9x
  slower than today on the profiling shape (100 rows x 1,000 draws: 71.9 s
  against 12.1 s).
- A lean loop, for each j the components k > j as one contiguous vector,
  brings it to parity (50 rows: 6.7 s against 7.0 s), agreeing to 5.6e-15.
  The method uses it.
- The realistic rows' indicator r was 0.46 to 0.62, and today's error there
  was at most 1.7e-15.
- The weibull axis method at tau_c = 0.25 was 2.9x today on a real AFT fit,
  before calibration.
- The projected total was about 2.6 core-hours, over twice the registered
  estimate, so the run was stopped and option B was chosen.

**Option B and its guard.**
- Check 6 runs ONCE for the four arms over 5 minutes: gaussian today, exact
  and hybrid, and weibull today and axis, all at 2,000 rows. The other arms
  keep the median of 3.
- If a single run's ratio to today lands between 1.6x and 2.4x (the 2x line,
  with noise), that arm runs twice more and the median is taken.
- The revised estimate is 1.6 core-hours, with a stop at 3.2 core-hours.

**The 2x rule covers the weibull arm.**
- If the axis method at the calibrated tau_c is still over 2x today on the
  real AFT fit, that goes to the maintainer with the options, such as a validated fast
  path for the non-gaussian families.
- No such fast path is built in advance.

**tau_c's search.** tau_c runs over 1, 1/2, 1/4, ..., 1/64, searched downward
from 1 and stopping at the first value that passes. So tau_c is capped at 1.
If none passes, tau_c is reported as not found, and the held-out check runs
at 1/64, with its rows over the cap integrated.

**No fallback to a method known to be wrong.**
- A continuous family without a CRPS method gets an informative error, not
  the uniform quadrature, which this work measured at up to 100% off on the
  positive families.
- The uniform quadrature stays only as the gaussian fast path.
- A test holds every family in illume's `ilm_family()` to one of three
  classes: a continuous family with a CRPS method, a discrete family, or a
  family `ilm_dist()` does not cover (the flexible survival families). So a
  family added upstream breaks limn's tests, not a user's scores.

## Amendment 4 (1 October, after run 3, before run 4; approved on review)

**Run 3** (on 42d4f4d; findings kept as `studies/findings/crps_continuous_run3.md`)
failed four rules: check 1's integrate() comparison, check 4's held-out rule,
and both checks' reference-failure counts. The failing rows were regenerated
from the study's seeds and scored against a
split reference: integrate() over the intervals between the components' own
quantiles.

1. **The registered reference was wrong.** integrate() over an infinite range
   silently misses narrow or far mass, returning about 0, or it fails:
   - gaussian near-degenerate and spiky rows: the exact form matched
     crps_mixnorm to 1e-15;
   - lognormal near-degenerate rows, and random rows with y near 1e-28;
   - beta rows at y = 1, where the registered reference was 10% to 260% off.
   On all of them the method matched the split reference.
2. **The method was wrong.**
   - (a) A lognormal with sdlog 14 to 16 (the scale-spread shape, which is not
     clipped for lognormal) was 4.1% and 20% low. The upper piece stopped at
     the 1 - 1e-16 quantile, but for so wide a lognormal the integrand
     e^u S(u)^2 peaks beyond it.
   - (b) One spiky lognormal row at S = 1,000 was NA: the whole-axis
     integrate() fallback failed.
   - (c) Two rows were marginal at tau_c = 1: weibull shape 0.3 (1.2e-10) and a
     tiny-shape beta (1.7e-10). The calibration grid had no such shapes.

**The changes.**
1. **The reference** for checks 1 and 4 is the split reference: integrate()
   at rel.tol 1e-13 over each interval between the components' quantiles at
   1e-14, 1e-8, 1e-4, 0.01, 0.25, 0.5, 0.75, 0.99 and their complements, plus
   y, from -Inf to Inf.
   - Check 1 (gaussian, S = 100) uses all the breakpoints.
   - Check 4 uses up to 1,000, a different split from the method's own
     fallback (120).
   - A reference that warns or fails is still a failed case.
2. **The reference is validated (the amendment's condition 1).** On every
   gaussian row given a split reference (check 1's 24 six-case rows at S = 100),
   it must agree with the exact form to 1e-12 or better. The worst is
   reported. The non-gaussian rows keep the raw-axis cross-check.
   - A pre-check gave at most 3.4e-15
     on those 24 rows, where the thinned version was 1e-7 on the staircase.
3. **The method's upper tail** for lognormal and weibull runs past hi until the
   integrand on the log axis is below 1e-16 of the integral and falling. The
   tail's panels are sized by the components that reach hi.
4. **The method's fallback** (over the cap, or untrusted quantiles) is split
   integration over the components' quantile breakpoints (up to 120), from -Inf
   to Inf, at rel.tol 1e-12. It is NA with a message on any warning or error.
   - A check: against an unthinned
     reference on spiky lognormal and weibull mixtures at S = 1,000, it was
     within 1.9e-14.
5. **tau_c's calibration** adds check 4's added shapes (weibull shape 0.3,
   beta(0.3, 0.3), loglogistic s = 0.8 to 1.9, sdlog 3, each alone and mixed)
   at S = 100, on seed 2853. The held-out set uses fresh seeds: 2858 for the
   shapes and 2859 for the random mixtures.
6. **The lognormal scale-spread shape** is clipped to the random-mixture bound
   (sdlog 3), as loglogistic and beta are. The unclipped case is kept as a
   named extra shape, `sdlog15_spread` (S = 100 and 1,000). It must pass the
   1e-10 rule (the amendment's condition 3).
7. **Run 4** goes from the top: the same checks and rules, and the same stop at
   3.2 core-hours. The findings report all four runs (condition 2), tau_c's new
   value, and check 6's weibull cost against the 2x rule (condition 4).

## What is reported

- Each check's worst case and its rule, pass or fail.
- tau and tau_c.
- The share of rows on the fast path at the profiling sizes.
- The costs.

NEWS gets a fix entry: what was wrong (equal panels sized by the widest
component), by how much (up to 0.77% relative on the sd-spread case, 0.1% on
the staircase case), and which families and cases.

Cost: about 1 core-hour, most of it in check 6 and the S = 4,000 runs.

2026-10-15: wording neutralised for publication; the design is unchanged.
