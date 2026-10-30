# Design: the discrete CRPS on wide supports (crps_wide)

Written before the method is built and before any validation run. It was
settled as option (i)(a) on 29 September 2026.

**Why:** the demand dry run found limn's exact discrete CRPS at 6 to 18 s a
row on PBS's large counts.
- Those rows' supports are 23,000 to 25,000 values wide. They take today's
  bracket path, which evaluates the mixture's CDF (over 1,000 draws) at about
  20,000 block ends a row.
- That path's rigorous bound is about 0.5 against CRPS values of 600 to 900:
  far more accuracy than a score needs, at great cost.

## The method

A row whose support, from the mixture's 1e-6 to 1 - 1e-6 quantiles, is wider
than **2,000 values** takes the wide path. That threshold replaces today's
1e4. Rows up to 2,000 wide are summed exactly, as now.

CRPS = sum over k of (F(k) - 1{k >= y})^2, plus the terms outside the support
as now. It is split into three parts:

1. **The neighbourhood of y**: the integers within h = 50 of y (101 values),
   clipped to the support, are summed exactly. That uses the running pmf sum,
   as today's exact path does. The jump in 1{k >= y} therefore never sits
   inside a Simpson panel.
2. **The smooth pieces** are F(k)^2 below the neighbourhood and
   (1 - F(k))^2 above it. Each piece, over the integers a..b, is summed by
   Euler-Maclaurin:
   sum_{k=a}^{e} g(k) = integral_a^e G + (g(a) + g(e))/2 + (G'(e) - G'(a))/12.
   - The integral is by composite Simpson on integer nodes a, a+d, ...,
     e = a + 2Nd, with N even.
   - G' comes from central differences, from F at a-1, a+1, e-1 and e+1.
   - The integers in (e, b], fewer than 2d of them, are summed exactly.
   - A piece shorter than 8 integers is summed exactly.
3. **The error estimate** is Richardson's: |S_d - S_2d| / 15 for each piece,
   the second Simpson sum on every other node, summed over both pieces and
   multiplied by a safety factor of 10.
   - d starts at ceiling(width / 400).
   - It halves, reusing its nodes, until the error estimate is at most 1e-6
     times the CRPS, or d = 1, where the sum is exact.
   - A row still above 1e-6 relative says so in the attribute "flagged".

**Recorded per row**, as attributes:
- `error_estimate` (the flag's basis);
- `bound`, the rigorous monotone-bracket width, d times each piece's rise in
  g. It is loose (O(d)), and kept for the record only;
- `flagged`, TRUE when error_estimate exceeds 1e-6 x the CRPS;
- `step`, the final d.

Rows summed exactly have error_estimate 0 and bound 0. The attribute that
`lmn_pred_scores` reports as `crps_bound` becomes the error estimate, and the
rigorous width joins it.

## Validation (pre-registered)

The exact sum, the pmf run over every integer of the support, is the
reference. It is computable up to widths of about 10,000 with 200 draws.

**Cases.** Each is a mixture built with `lmn_new_pred` from known draw
parameters. S is 200 draws unless stated. Each case has 20 rows, with y at
the mixture's 10%, 50% and 90% quantiles.

- **(A) Typical:** negative binomial draws, means 2,000 to 60,000 (spanning
  widths 1,000 to 10,000 by the size parameter), size 5 to 50, the draws'
  log-means spread by SD 0.05.
- **(B) Draws that differ widely:**
  - bimodal: half the draws at mean m, half at 2m;
  - heavy-tailed: log-means drawn from a t with 3 df, scale 0.3.
- **(C) y far in the tails:**
  - at the mixture's 1e-5 and 1 - 1e-5 quantiles;
  - outside the support on either side (the recomputed path);
  - exactly at the support's ends.
- **(D) Steep F:** near-Poisson negative binomial draws (size 1e4), and the
  Poisson family itself, with means chosen to give widths of 2,000 to 5,000.
- **(E) The threshold:** means tuned so the width is 1,990 to 2,010. Those
  rows are summed exactly (width 2,000 or less) or by the wide path (over
  2,000), each compared with the exact sum.

**Pass rules, fixed now:**
1. |estimate - exact| <= 1e-6 x exact, in every row.
2. |estimate - exact| <= error_estimate, in 100% of wide-path rows. Any row
   where the estimate undershoots fails the rule, and is reported.
3. In (E), rows either side of 2,000 agree with the exact sum within rule 1.
4. The rigorous `bound` is reported. The estimate lying inside the bracket is
   reported, not required.

**Speed:** the dry run's A2 PBS rows, 12 rows of one series (the 58,700-a-
month series) and 12 of two others. Current path against new, on the same
rows with the same draws (1,000), seconds per row. Their CRPS values are
compared as well: the new path must agree with the current one within the
current path's rigorous bound plus the new error_estimate.

The findings report every case, whichever way it goes. The runs use 1 to 2
cores on the pinned library `demand-442639f`, which holds tsibbledata for the
PBS rows.

## Addendum (29 September, after the first run stopped, before the rerun)

The first run stopped at a script fault before any row was compared. The
design's cases (A) and (B, heavy-tailed) gave rows mostly too wide for the
exact reference:
- (A)'s means of 2,000 to 60,000 kept only 6 of 21 rows within 12,000 values;
- the heavy-tailed draws (t3, scale 0.3) kept none, and an empty case stopped
  the script.

Both cases now choose each row's mean by bisection to a target width. That
is the width range the design names, which the fixed means could not give:
- **(A):** targets 1,000, 2,500, 4,000, 6,000, 8,000 and 10,000 for each size
  (5, 20 and 50): 18 rows.
- **(B, heavy-tailed):** targets 1,000 to 10,000 in 10 even steps, with the
  same t3 draws (scale 0.3).

Nothing else changes: the other cases, the pass rules and the speed
comparison. A case that still yields no rows is reported as empty instead of
stopping the script.
