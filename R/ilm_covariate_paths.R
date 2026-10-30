#' Draws of a covariate's future, from its own model
#'
#' When a covariate is not known in advance, its uncertainty belongs in the
#' forecast. `ilm_covariate_paths()` simulates future paths of the covariate
#' from an illume model of it, and returns `newdata` stacked into one block of
#' rows per path, marked by a `.draw` column. [ilm_forecast()] and
#' [ilm_predict()] then mix over the paths; see *Covariates given as draws* in
#' [ilm_predict()].
#'
#' @section How the paths are drawn:
#' Each path is one draw from the covariate model's predictive distribution of
#' a new observation. That draw has three parts:
#' * a draw of the model's parameters, under `uncertainty`;
#' * each series' latent path continued over the new times, from its fitted
#'   state (`groups = "fitted"`);
#' * a new observation at each time.
#'
#' So a path is joint over its times. Together, the paths carry both the
#' covariate's future shocks and the uncertainty in its model's parameters.
#' Give `nsim` the same value as the forecast's own `nsim`: parameter draw s
#' of the forecast is then paired with path s.
#'
#' @section Several covariates:
#' Give a list of fits, one per covariate. Each covariate is drawn from its own
#' model, independently of the others. Covariates that move together would
#' need a joint model of them, which is not in limn yet.
#'
#' @param fit A fitted `"ilm_model"` of the covariate. Its response names the
#'   column of `newdata` to draw. Or a list of such fits, one per covariate.
#' @param newdata The rows to forecast, as for [ilm_forecast()], with the
#'   columns the covariate models need. The drawn columns may be missing, as
#'   [ilm_future()] leaves them.
#' @param nsim The number of paths.
#' @param uncertainty How the covariate model's own parameters enter:
#'   `"joint"` (the default), `"fixed"`, `"latent"` or `"none"`, as in
#'   [ilm_predict()].
#' @param seed Seed for the draws. The caller's random-number state is restored
#'   afterwards.
#' @return `newdata` repeated once per path, with a `.draw` column (1 to
#'   `nsim`) and the drawn columns filled in.
#' @seealso [ilm_as_draws()] for paths simulated by another tool;
#'   [ilm_future()] for the rows.
#' @examples
#' set.seed(1)
#' d <- data.frame(id = factor(rep(c("a", "b"), each = 30)), t = rep(1:30, 2))
#' d$x <- ave(rnorm(60), d$id, FUN = function(e) stats::filter(e, 0.8, "recursive"))
#' d$y <- 2 + 0.5 * d$x + rnorm(60, 0, 0.3)
#' fx <- illume::ilm_model(x ~ 1, data = d, family = "gaussian",
#'                         ar = illume::ilm_ar1(~ t | id), verbose = FALSE)
#' fy <- illume::ilm_model(y ~ x, data = d, family = "gaussian", verbose = FALSE)
#' fut <- ilm_future(d, time = "t", group = "id", h = 3, response = "y")
#' paths <- ilm_covariate_paths(fx, fut, nsim = 200)
#' head(paths)
#' fc <- ilm_forecast(fy, paths, nsim = 200)
#' as.data.frame(fc)
#' @export
ilm_covariate_paths <- function(fit, newdata, nsim = 1000L,
                                uncertainty = c("joint", "fixed", "latent", "none"),
                                seed = 1L) {
  fits <- if (inherits(fit, "ilm_model")) list(fit) else fit
  if (!is.list(fits) || !length(fits) || !all(vapply(fits, inherits, TRUE, "ilm_model")))
    stop("`fit` must be a fitted ilm_model of the covariate, or a list of them.", call. = FALSE)
  uncertainty <- match.arg(uncertainty)
  nsim <- lmn_check_count(nsim, "nsim")
  if (!is.data.frame(newdata))
    stop("`newdata` must be a data frame of the rows to forecast.", call. = FALSE)
  if (!is.null(newdata$.draw))
    stop("`newdata` already has a `.draw` column; give the rows once, without it.", call. = FALSE)
  vars <- vapply(fits, lmn_response_name, "")
  if (anyDuplicated(vars))
    stop("two of the fits model the same covariate, '", vars[duplicated(vars)][1], "'.",
         call. = FALSE)
  n <- nrow(newdata)
  out <- cbind(.draw = rep(seq_len(nsim), each = n),
               newdata[rep(seq_len(n), nsim), , drop = FALSE])
  for (k in seq_along(fits)) {
    p <- ilm_predict(fits[[k]], newdata, target = "response", groups = "fitted",
                     uncertainty = uncertainty, nsim = nsim, seed = lmn_sub_seed(seed, 100L * k))
    if (!is.null(p$dist$categories))
      stop("'", vars[k], "' is an ordered outcome; paths of categorical covariates are not ",
           "in limn yet.", call. = FALSE)
    out[[vars[k]]] <- as.vector(lmn_simulate_paths(p, nsim, lmn_sub_seed(seed, 100L * k + 1L)))
  }
  rownames(out) <- NULL
  out
}

## One simulated outcome per path, for every row: an n x nsim matrix. Path j
## uses parameter draw ((j - 1) mod S) + 1, so a prediction with fewer draws
## than paths (the estimates alone, say) still gives every path its own
## outcomes.
#' @keywords internal
#' @noRd
lmn_simulate_paths <- function(pred, nsim, seed) {
  n <- nrow(pred$lp[[1]]); S <- length(pred$weights)
  lmn_with_seed(seed, vapply(seq_len(nsim), function(j) {
    a <- lmn_draw_args(pred, seq_len(n), (j - 1L) %% S + 1L)
    as.numeric(pred$dist$r(a$lp, a$par))
  }, numeric(n)))
}

## The column a fit models: its response, which must be a plain column.
#' @keywords internal
#' @noRd
lmn_response_name <- function(fit) {
  f <- stats::formula(fit)
  if (length(f) < 3L) stop("the covariate model has no response.", call. = FALSE)
  lhs <- f[[2L]]
  if (!is.name(lhs))
    stop("the covariate model's response must be a column of the data, not an expression (",
         paste(deparse(lhs), collapse = ""), "); model the column itself.", call. = FALSE)
  as.character(lhs)
}

#' Covariate paths from another tool, as draws
#'
#' Turns simulated paths of covariates made elsewhere into `newdata`'s blocks
#' of draws, marked by a `.draw` column, which [ilm_forecast()] and
#' [ilm_predict()] mix over (see *Covariates given as draws* in
#' [ilm_predict()]). The paths might come from fable's `generate()`, a weather
#' ensemble, or a scenario generator. For a covariate modelled with illume,
#' [ilm_covariate_paths()] draws them directly.
#'
#' Every path must give a value for every row of `newdata`, matched on the
#' `by` columns, and paths are numbered in the order they first appear.
#' Columns of `newdata` that are not drawn are copied into every block.
#'
#' @param paths A data frame in long form: one row per path and row of
#'   `newdata`, with a path identifier, the `by` columns and the simulated
#'   values.
#' @param newdata The rows to forecast.
#' @param by The columns, in both, that identify a row of `newdata`, such as
#'   its series and time.
#' @param path The column of `paths` that identifies a path. The default,
#'   `".rep"`, is the name fable's `generate()` uses (fpp3, section 5.5).
#' @param values The drawn columns, as a named character vector mapping each
#'   column of `newdata` to its column in `paths`: for fable's `generate()`,
#'   whose simulated values are in `.sim`, `c(x = ".sim")`. By default every
#'   column of `paths` that is not in `by` or `path` and has the same name in
#'   `newdata`.
#' @return `newdata` repeated once per path, with a `.draw` column and the drawn
#'   columns filled in from the paths.
#' @references Hyndman, R. J. and Athanasopoulos, G. (2021). *Forecasting:
#'   Principles and Practice*, 3rd edition. OTexts. Section 5.5,
#'   "Distributional forecasts and prediction intervals": fable's `generate()`
#'   and its `.rep` and `.sim` columns (online edition,
#'   <https://otexts.com/fpp3/toolbox.html#sec-prediction-intervals>, read
#'   27 October 2026).
#' @seealso [ilm_covariate_paths()] for paths from an illume model; [ilm_forecast()],
#'   which mixes over the blocks.
#' @examples
#' fut <- data.frame(id = rep(c("a", "b"), each = 2), t = rep(31:32, 2), x = NA)
#' sims <- expand.grid(t = 31:32, id = c("a", "b"), .rep = 1:3)
#' sims$.sim <- rnorm(nrow(sims))
#' ilm_as_draws(sims, fut, by = c("id", "t"), values = c(x = ".sim"))
#' @export
ilm_as_draws <- function(paths, newdata, by, path = ".rep", values = NULL) {
  if (!is.data.frame(paths) || !is.data.frame(newdata))
    stop("`paths` and `newdata` must be data frames.", call. = FALSE)
  paths <- as.data.frame(paths); newdata <- as.data.frame(newdata)
  if (!is.null(newdata$.draw))
    stop("`newdata` already has a `.draw` column; give the rows once, without it.", call. = FALSE)
  miss <- c(setdiff(by, names(newdata)), setdiff(c(by, path), names(paths)))
  if (length(miss))
    stop("column ", paste0("'", unique(miss), "'", collapse = ", "), " is missing from ",
         if (all(miss %in% names(newdata))) "`paths`" else "the data", ".", call. = FALSE)
  if (is.null(values)) {
    v <- setdiff(intersect(names(paths), names(newdata)), c(by, path))
    values <- stats::setNames(v, v)
  }
  if (!length(values) || is.null(names(values)) || any(!nzchar(names(values))))
    stop("say which columns are drawn: `values = c(x = \".sim\")` maps newdata's x to ",
         "the paths' .sim.", call. = FALSE)
  if (!all(values %in% names(paths)))
    stop("`paths` has no column ", paste0("'", setdiff(values, names(paths)), "'", collapse = ", "),
         ".", call. = FALSE)
  key <- function(d) do.call(paste, c(lapply(d[by], as.character), sep = "\r"))
  k_nd <- key(newdata)
  if (anyDuplicated(k_nd))
    stop("`by` must identify each row of `newdata`; some rows share ",
         paste(by, collapse = ", "), ".", call. = FALSE)
  ids <- unique(paths[[path]])
  k_p <- paste(match(paths[[path]], ids), key(paths), sep = "\n")
  if (anyDuplicated(k_p))
    stop("a path gives more than one value for the same row of `newdata`.", call. = FALSE)
  want <- paste(rep(seq_along(ids), each = nrow(newdata)), rep(k_nd, length(ids)), sep = "\n")
  pos <- match(want, k_p)
  if (anyNA(pos)) {
    j <- which(is.na(pos))[1]
    stop("path ", format(ids[(j - 1L) %/% nrow(newdata) + 1L]), " has no value for row ",
         (j - 1L) %% nrow(newdata) + 1L, " of `newdata`; every path must cover every row.",
         call. = FALSE)
  }
  out <- cbind(.draw = rep(seq_along(ids), each = nrow(newdata)),
               newdata[rep(seq_len(nrow(newdata)), length(ids)), , drop = FALSE])
  for (v in names(values)) out[[v]] <- paths[[values[[v]]]][pos]
  rownames(out) <- NULL
  out
}
