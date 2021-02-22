# fitted values for sfacross ----------

#' @title Fitted frontier value
#'
#' @name fitted
#'
#' @aliases fitted.sfacross
#'
#' @description This method returns the fitted frontier values
#' from stochastic frontier models estimated with \code{\link{sfacross}} or
#' \code{\link{lcmcross}}.
#' @param object A (latent class) stochastic frontier model returned by
#' \code{\link{sfacross}} or \code{\link{lcmcross}}.
#' @param ... Currently ignored.
#'
#' @return \code{fitted} returns a vector of fitted values in the case of object
#'   of class \code{sfacross}. When class of object is \code{lcmcross} a data
#'   frame containing the fitted values for each latent class is returned. The
#'   fitted values are ordered in the same way as the corresponding observations
#'   in the data set used for the estimation.
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
#'
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
#' fitted(tl_u_ts)
#'
#' ## Using data on eighty-two countries production (DGP)
#'
#' # LCM Cobb Douglas (production function) half normal distribution
#'
#' cb_2c_h <- lcmcross(formula = ly ~ lk + ll + yr, udist = 'hnormal',
#'                     data = worldprod, S = 1, method = 'ucminf')
#'
#' fitted(cb_2c_h)
#' @keywords methods fitted
fitted.sfacross <- function(object, ...) {
  object$dataTable$mleFitted
}

# fitted values for lcmcross ----------
#' @rdname fitted
#' @aliases fitted.lcmcross
#' @export
fitted.lcmcross <- function(object, ...) {
  if (object$nClasses == 2) {
    select(object$dataTable, "mleFitted_c1",
      "mleFitted_c2")
  } else {
    if (object$nClasses == 3) {
      select(object$dataTable, "mleFitted_c1",
        "mleFitted_c2", "mleFitted_c3")
    } else {
      if (object$nClasses == 4) {
        select(object$dataTable, "mleFitted_c1",
          "mleFitted_c2", "mleFitted_c3",
          "mleFitted_c4")
      } else {
        if (object$nClasses == 5) {
          select(object$dataTable, "mleFitted_c1",
          "mleFitted_c2", "mleFitted_c3",
          "mleFitted_c4", "mleFitted_c5")
        }
      }
    }
  }
}
