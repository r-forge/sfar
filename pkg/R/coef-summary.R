# coefficients from summary.sfacross ----------

#' @title coef methods for class \code{summary} of frontier objects
#'
#' @name coef
#'
#' @aliases coef.summary.sfacross
#'
#' @description Extract the coefficients, their standard errors, z-values, and
#'   (asymptotic) P-values from models returned by the \code{summary} methods
#'   for objects of class \code{sfacross} or \code{lcmcross}. Optionally, the
#'   gradient is also returned and or with the coefficients confidence
#'   intervals.
#'
#' @param object An object of either class \code{"summary.sfacross"}
#' or class \code{"summary.lcmcross"}.
#' @param ... Currently ignored.
#'
#' @return The coef method for objects of class summary.sfacross or
#'   summary.lcmcross returns a matrix with four to seven columns (depending on
#'   the option \code{grad} and \code{ci}). These columns contain the estimated
#'   coefficients, their standard errors, z-values, and (asymptotic) P-values.
#'   Optionally, the gradient is also returned and or with the coefficients
#'   confidence intervals.
#' @export
#'
#' @examples
#' ## Data on fossil fuel fired steam electric power generation plants in the
#' ## United States
#'
#' # Translog SFA (cost function) truncated normal with scaling property
#'
#' tl_u_ts <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
#' log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
#' I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
#' udist = "tnormal", muhet = ~ regu, uhet = ~ regu, data = utility, S = -1,
#' scaling = TRUE, method = "mla")
#'
#' coef(summary(tl_u_ts))
#'
#' ## Using data on eighty-two countries production (DGP)
#'
#' # LCM Cobb Douglas (production function) half normal distribution
#'
#' cb_2c_h <- lcmcross(formula = ly ~ lk + ll + yr, udist = "hnormal",
#'                     data = worldprod, S = 1, method = "ucminf")
#'
#' coef(summary(cb_2c_h))
#' @keywords methods coefficients summary
coef.summary.sfacross <- function(object, ...) {
  object$mleRes
}


# coefficients from summary.lcmcross ----------
#' @rdname coef
#' @aliases coef.summary.lcmcross
#' @export
coef.summary.lcmcross <- function(object, ...) {
  object$mleRes
}
