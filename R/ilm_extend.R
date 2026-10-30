## ---- extension points: other ways of fitting illume's models -------------------------
##
## limn predicts, scores and cross-validates any fit that answers these two
## generics. An illume fit answers them through the methods here, which are
## what limn did before the generics existed. A package that fits illume's
## models another way (by sampling a posterior, say) registers its own
## methods, and every verb in limn then works on its fits. limn neither
## depends on nor suggests such a package: it owns the connection.

#' The parameter draws behind a prediction
#'
#' `ilm_predict_draws()` is the step of [ilm_predict()] that draws the
#' model's parameters. Every prediction, forecast, score, calibration and
#' cross-validation in limn is computed from these draws. You will not
#' normally call it yourself: it is a generic, so that a package fitting
#' illume's models another way can supply the draws from its own fits, and
#' limn then works on those fits unchanged.
#'
#' The method for illume's fits draws from the fit's normal approximation,
#' with [illume::ilm_draws()]. A method for another class of fit must return
#' the same shape: the class of fit must inherit from `"ilm_model"`, so that
#' illume's [illume::ilm_dist()], [illume::ilm_varcorr()],
#' [illume::ilm_matrices()] and [illume::ilm_cells()] accept it, and the draws
#' must be laid out as [illume::ilm_draws()] lays them out.
#'
#' @param fit A fitted model: an `"ilm_model"`, or a class that inherits from
#'   it.
#' @param nsim The number of draws to return.
#' @param seed The seed the draws are made from. The same seed must give the
#'   same draws, whatever else is running.
#' @param uncertainty Which parameters are drawn, as in [ilm_predict()]:
#'   `"joint"` (all of them), `"fixed"` (the variances held at their
#'   estimates), `"latent"` (every parameter held, so only the random effects
#'   and latent cells are drawn) or `"none"` (one draw, at the estimates). A
#'   method may refuse a value it has no meaning for, with an error that says
#'   which values it supports.
#' @param natural Whether each draw's variances are also needed on their
#'   natural scale. limn asks for them only when a unit's effects are drawn
#'   afresh or averaged over.
#' @param ... Passed on to methods.
#' @return An object of class `"ilm_draws"`, as [illume::ilm_draws()]
#'   returns:
#'   * `draws`: a matrix with one row per parameter, in illume's order and on
#'     its internal scale, and `nsim` columns (one for `uncertainty = "none"`);
#'   * `mode`: the estimates;
#'   * `map`: what each row of `draws` is, as [illume::ilm_draws()] gives it;
#'   * `natural`: when `natural = TRUE`, each draw's variances on their
#'     natural scale, which at the estimates equal [illume::ilm_varcorr()]'s.
#'
#'   [ilm_predict()] checks this shape, and stops with a message saying what
#'   is wrong when a method returns something else.
#' @eval lmn_diag_rd("draws_shape")
#' @seealso [ilm_cv_refit()], the other extension point.
#' @examples
#' set.seed(1)
#' d <- data.frame(g = factor(rep(letters[1:8], each = 5)), x = rnorm(40))
#' d$y <- rpois(40, exp(0.5 + 0.3 * d$x))
#' fit <- illume::ilm_model(y ~ x + (1 | g), data = d, family = "poisson",
#'                          verbose = FALSE)
#' dr <- ilm_predict_draws(fit, nsim = 10, seed = 1, uncertainty = "joint",
#'                         natural = FALSE)
#' dim(dr$draws)
#' @export
ilm_predict_draws <- function(fit, nsim, seed, uncertainty, natural, ...) {
  UseMethod("ilm_predict_draws")
}

#' @rdname ilm_predict_draws
#' @export
ilm_predict_draws.ilm_model <- function(fit, nsim, seed, uncertainty, natural, ...) {
  lmn_param_draws(fit, uncertainty, nsim, seed, illume::ilm_varcorr(fit), natural)
}

#' Refit a model to a fold's training rows
#'
#' `ilm_cv_refit()` is the step of [ilm_cv()] that refits the model to each
#' fold's training rows. You will not normally call it yourself: it is a
#' generic, so that a package fitting illume's models another way can refit
#' its own fits, and [ilm_cv()] then cross-validates them unchanged.
#'
#' The method for illume's fits is [illume::ilm_refit()], which repeats the
#' model's own call on the new rows. A method for another class of fit must
#' return a fit that [ilm_predict()] accepts, and must give the same fit
#' whether it runs in the main R session or in one of the worker processes
#' that `ilm_cv(ncores = )` starts. Each worker loads limn and the package
#' that registered the method, and runs one thread.
#'
#' @param fit A fitted model: an `"ilm_model"`, or a class that inherits from
#'   it.
#' @param data The fold's training rows.
#' @param ... Passed on to methods.
#' @return The model refitted to `data`.
#' @seealso [ilm_predict_draws()], the other extension point.
#' @examples
#' set.seed(1)
#' d <- data.frame(g = factor(rep(letters[1:8], each = 5)), x = rnorm(40))
#' d$y <- rpois(40, exp(0.5 + 0.3 * d$x))
#' fit <- illume::ilm_model(y ~ x + (1 | g), data = d, family = "poisson",
#'                          verbose = FALSE)
#' f1 <- ilm_cv_refit(fit, d[d$g != "a", ])
#' @export
ilm_cv_refit <- function(fit, data, ...) {
  UseMethod("ilm_cv_refit")
}

#' @rdname ilm_cv_refit
#' @export
ilm_cv_refit.ilm_model <- function(fit, data, ...) {
  illume::ilm_refit(fit, data = data)
}

## The shape ilm_predict() needs from ilm_predict_draws(): checked here so that
## a method returning something else stops with what is wrong, rather than
## deep in the assembly.
#' @keywords internal
#' @noRd
lmn_check_draws <- function(dr, nsim, uncertainty, natural) {
  bad <- function(what) lmn_signal("draws_shape", what = what, type = "error")
  if (!is.list(dr) || !is.matrix(dr$draws) || !is.numeric(dr$draws))
    bad("it has no numeric matrix `draws`")
  want <- if (uncertainty == "none") 1L else nsim
  if (ncol(dr$draws) != want)
    bad(paste0("`draws` has ", ncol(dr$draws), " columns, not ", want))
  if (!is.data.frame(dr$map) || !all(c("row", "block", "term") %in% names(dr$map)))
    bad("it has no `map` with columns row, block and term")
  if (nrow(dr$map) != nrow(dr$draws))
    bad(paste0("`map` describes ", nrow(dr$map), " rows, but `draws` has ", nrow(dr$draws)))
  if (natural && is.null(dr$natural))
    bad("natural-scale variances were asked for, and `natural` is missing")
  invisible(TRUE)
}

## The packages, other than limn, whose methods the extension generics dispatch
## to for this fit: each worker of ilm_cv() loads them, since a method
## registered by a package that is not loaded is not there to dispatch to.
#' @keywords internal
#' @noRd
lmn_method_packages <- function(fit, generics = c("ilm_predict_draws", "ilm_cv_refit")) {
  pk <- character(0)
  for (g in generics) for (cl in class(fit)) {
    m <- utils::getS3method(g, cl, optional = TRUE, envir = asNamespace("limn"))
    if (!is.null(m)) {
      e <- environment(m)
      if (isNamespace(e)) pk <- c(pk, getNamespaceName(e))
      break
    }
  }
  setdiff(unique(pk), "limn")
}
