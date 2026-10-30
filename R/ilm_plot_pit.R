#' Histogram of probability integral transforms
#'
#' Plots PITs from [ilm_pit()] as a histogram on the density scale, where a
#' correct predictive distribution gives a flat histogram at height 1. The
#' shaded band is where each bar falls with probability `level` if the PITs
#' really are uniform.
#'
#' @section Reading the plot:
#' * **Flat, within the band:** the outcomes fell where their predictive
#'   distributions said they would. The band holds for each bar on its own,
#'   so with many bars one may stray outside it by chance.
#' * **U shape**, high at both ends: too many outcomes fell in the tails. This
#'   often means predictive distributions too narrow.
#' * **Hump**, high in the middle: intervals too wide on average.
#' * **Triangle**, sloping one way: the predictions are biased, too high if the
#'   bars rise towards 0, too low if they rise towards 1.
#'
#' The shapes are those Gneiting, Balabdaoui and Raftery (2007, p. 253)
#' describe; which way a triangle points is limn's own reading of it. For a
#' verdict with what to do about it, use
#' [ilm_calibration()][ilm_calibration.ilm_pred].
#'
#' @param x PITs from [ilm_pit()], or any numbers between 0 and 1.
#' @param bins Number of bars.
#' @param level Coverage of the shaded band under uniformity.
#' @param ... Passed to [tinyplot::tinyplot()].
#' @return The data behind the plot, invisibly: each bar's limits, count and
#'   height, and the band.
#' @references Gneiting, T., Balabdaoui, F. and Raftery, A. E. (2007).
#'   Probabilistic forecasts, calibration and sharpness. *Journal of the Royal
#'   Statistical Society B* 69, 243-268. Section 3.1, p. 253.
#' @seealso [ilm_pit()], which computes the PITs;
#'   [ilm_calibration()][ilm_calibration.ilm_pred] for a verdict.
#' @examples
#' u <- runif(500)
#' ilm_plot_pit(u)
#' @export
ilm_plot_pit <- function(x, bins = 10L, level = 0.95, ...) {
  u <- as.numeric(x)
  u <- u[is.finite(u)]
  if (!length(u)) stop("`x` has no PITs to plot.", call. = FALSE)
  if (any(u < 0 | u > 1))
    stop("PITs lie between 0 and 1; `x` has values outside.", call. = FALSE)
  bins <- lmn_check_count(bins, "bins", 2L)
  n <- length(u)
  br <- seq(0, 1, length.out = bins + 1L)
  cnt <- tabulate(findInterval(u, br, rightmost.closed = TRUE, all.inside = TRUE), bins)
  band <- stats::qbinom(c((1 - level) / 2, (1 + level) / 2), n, 1 / bins) * bins / n
  d <- data.frame(lower = br[-(bins + 1L)], upper = br[-1L], count = cnt,
                  density = cnt * bins / n, band_lower = band[1], band_upper = band[2])
  top <- max(d$density, band[2]) * 1.08
  tinyplot::tinyplot(xmin = 0, xmax = 1, ymin = band[1], ymax = band[2],
                     type = tinyplot::type_rect(), fill = "grey88", col = NA,
                     xlim = c(0, 1), ylim = c(0, top), xlab = "PIT", ylab = "density",
                     ...)
  tinyplot::tinyplot_add(xmin = d$lower, xmax = d$upper, ymin = 0, ymax = d$density,
                         type = tinyplot::type_rect(), fill = "#2a6f9780", col = "#2a6f97")
  ## a horizontal line still needs a point to anchor the layer, or it is not drawn
  tinyplot::tinyplot_add(x = 0.5, y = 1, type = tinyplot::type_hline(h = 1),
                         lty = 2, col = "grey30")
  invisible(d)
}
