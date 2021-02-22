# residuals from sfacross ----------

#' @title Residuals of frontier estimation
#'
#' @name residuals
#'
#' @aliases residuals.sfacross
#'
#' @description This method returns the residuals values
#' from stochastic frontier models estimated with \code{\link{sfacross}} or
#' \code{\link{lcmcross}}.
#'
#' @param object A (latent class) stochastic frontier model returned by
#' \code{\link{sfacross}} or \code{\link{lcmcross}}.
#' @param ... Currently ignored.
#'
#' @return \code{residuals} returns a vector of residuals values in the case of
#'   object of class \code{sfacross}. When class of object is \code{lcmcross} a
#'   data frame containing the residuals values for each latent class is
#'   returned. The residuals values are ordered in the same way as the
#'   corresponding observations in the data set used for the estimation.
#' @export
#'
#' @author K Hervé Dakpo, Yann Desjeux and Laure Latruffe
#'
#' @seealso
#'
#' \code{\link{sfacross}}, for the stochastic frontier analysis model fitting
#' function.
#'
#' \code{\link{lcmcross}}, for the latent class stochastic frontier analysis
#' model fitting function.
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
#' udist = 'tnormal', muhet = ~ regu, uhet = ~ regu, data = utility, S = -1,
#' scaling = TRUE, method = 'mla')
#'
#' residuals(tl_u_ts)
#'
#' ## Using data on eighty-two countries production (DGP)
#'
#' # LCM Cobb Douglas (production function) half normal distribution
#'
#' cb_2c_h <- lcmcross(formula = ly ~ lk + ll + yr, udist = 'hnormal',
#'                     data = worldprod, S = 1, method = 'ucminf')
#'
#' residuals(cb_2c_h)
#' @keywords methods residuals
residuals.sfacross <- function(object, ...) {
  object$dataTable$mleResiduals
}

# residuals from lcmcross ----------
#' @rdname residuals
#' @aliases residuals.lcmcross
#' @export
residuals.lcmcross <- function(object, ...) {
  if (object$nClasses == 2) {
    select(object$dataTable, "mleResiduals_c1",
      "mleResiduals_c2")
  } else {
    if (object$nClasses == 3) {
      select(object$dataTable, "mleResiduals_c1",
        "mleResiduals_c2", "mleResiduals_c3")
    } else {
      if (object$nClasses == 4) {
        select(object$dataTable, "mleResiduals_c1",
          "mleResiduals_c2", "mleResiduals_c3",
          "mleResiduals_c4")
      } else {
        if (object$nClasses == 5) {
          select(object$dataTable, "mleResiduals_c1",
          "mleResiduals_c2", "mleResiduals_c3",
          "mleResiduals_c4", "mleResiduals_c5")
        }
      }
    }
  }
}
