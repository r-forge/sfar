# information criteria for sfacross ----------

#' @title Information criteria
#'
#' @name ic
#'
#' @method ic sfacross
#'
#' @aliases ic.sfacross
#'
#' @description This method returns information criterion from (latent class)
#' stochastic frontier models estimated with \code{\link{sfacross}} or
#' \code{\link{lcmcross}}
#'
#' @param object A (latent class) stochastic frontier model returned by
#' \code{\link{sfacross}} or \code{\link{lcmcross}}.
#' @param IC charater string. Information criterion measure. Three criterion are
#'   available:
#'
#' - \code{"AIC"} for Akaike information criterion
#'
#' - \code{"BIC"} for Bayesian information criterion
#'
#' - \code{"HQIC"} for Hannan-Quinn information criterion
#'
#' @param ... Currently ignored.
#'
#' @return \code{ic} returns the value of the information criterion (AIC, BIC or
#'   HQIC) of the maximum likelihood coefficients.
#'
#' @details The different information criterion are computed as follows:
#'
#'  - AIC: \eqn{-2 \log{LL} + 2 * K}
#'
#'  - BIC: \eqn{-2 \log{LL} + \log{N} * K}
#'
#'  - HQIC: \eqn{-2 \log{LL} + 2 \log{\left[\log{N}\right]} * K}
#'
#'  where LL is the maximum likelihood value, K the number of parameters
#'  estimated and N the number of observations.
#' @export
#' @export ic
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
#' udist = "tnormal", muhet = ~ regu, uhet = ~ regu, data = utility, S = -1,
#' scaling = TRUE, method = "mla")
#'
#' ic(tl_u_ts, IC = "AIC")
#'
#' ## Using data on eighty-two countries production (DGP)
#'
#' # LCM Cobb Douglas (production function) half normal distribution
#'
#' cb_2c_h <- lcmcross(formula = ly ~ lk + ll + yr, udist = "hnormal",
#'                     data = worldprod, S = 1, method = "ucminf")
#'
#' ic(cb_2c_h, IC = "AIC")
#' @keywords methods AIC BIC HQIC
ic.sfacross <- function(object, IC = "AIC", ...) {
  if (!(IC %in% c("AIC", "BIC", "HQIC"))) {
    stop("Unknown information criteria: ", paste(IC), call. = FALSE)
  }
  if (IC == "AIC") {
    -2 * object$mleLoglik + 2 * object$nParm
  } else {
    if (IC == "BIC") {
      -2 * object$mleLoglik + log(object$Nobs) * object$nParm
    } else {
      if (IC == "HQIC") {
        -2 * object$mleLoglik + 2 * log(log(object$Nobs)) *
          object$nParm
      }
    }
  }
}

# information criteria for lcmcross ----------
#' @rdname ic
#' @method ic lcmcross
#' @aliases ic.lcmcross
#' @export
ic.lcmcross <- function(object, IC = "AIC", ...) {
  if (!(IC %in% c("AIC", "BIC", "HQIC"))) {
    stop("Unknown information criteria: ", paste(IC), call. = FALSE)
  }
  if (IC == "AIC") {
    -2 * object$mleLoglik + 2 * object$nParm
  } else {
    if (IC == "BIC") {
      -2 * object$mleLoglik + log(object$Nobs) * object$nParm
    } else {
      if (IC == "HQIC") {
        -2 * object$mleLoglik + 2 * log(log(object$Nobs)) *
          object$nParm
      }
    }
  }
}
