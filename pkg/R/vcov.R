# variance covariance matrix for sfacross ----------

#' @title Variance-covariance matrix for frontier objects
#'
#' @name vcov
#'
#' @aliases vcov.sfacross
#'
#' @description Extract the variance-covariance matrix of the maximum likelihood
#'   coefficients of frontier objects (\code{\link{sfacross}} or
#'   \code{\link{lcmcross}}).
#'
#' @param object A (latent class) stochastic frontier model returned by
#'   \code{\link{sfacross}} or \code{\link{lcmcross}}
#' @param extraPar Logical. Default = \code{FALSE}. If \code{TRUE}, and only in
#'   the case of object of class \code{"sfacross"}, variances and covariances of
#'   additional parameters are returned:
#'
#'   \code{sigmaSq} = \code{sigmauSq} + \code{sigmavSq},
#'
#'   \code{lambdaSq} = \code{sigmauSq}/\code{sigmavSq},
#'
#'   \code{sigmauSq} = \eqn{\exp{(Wu)}},
#'
#'   \code{sigmavSq} = \eqn{\exp{(Wv)}},
#'
#'   \code{sigma} = \code{sigmaSq}^0.5,
#'
#'   \code{lambda} = \code{lambdaSq}^0.5,
#'
#'   \code{sigmau} = \code{sigmauSq}^0.5,
#'
#'   \code{sigmav} = \code{sigmavSq}^0.5,
#'
#'   \code{gamma} = \code{sigmauSq}/(\code{sigmauSq} + \code{sigmavSq})
#' @param ... Currently ignored
#'
#' @return \code{vcov} returns the covariance matrix of the maximum likelihood
#'   coefficients.
#'
#' @export
#'
#' @details the variance-covariance matrix is obtained by the inversion of the
#'   negative hessian matrix. Depending on the distribution and the
#'   \code{"hessianType"} option the analytical/numeric hessian or the bhhh
#'   hessian or the robust hessian matrix is evaluated. When \code{"extraPar =
#'   TRUE"} the variance-covariance of the additional parameters (in the case of
#'   \code{\link{sfacross}}) is obtained using the delta method.
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
#' ## Data on Spanish dairy farms
#'
#' # Cobb Douglas (production function) half normal, exponential and Rayleigh
#' # distributions
#'
#' cb_s_h <- sfacross(
#'   formula = YIT ~ X1 + X2 + X3 + X4, udist = "hnormal", data
#'   = dairyspain, S = 1, method = "bfgs"
#' )
#'
#' vcov(cb_s_h)
#'
#' cb_s_e <- sfacross(
#'   formula = YIT ~ X1 + X2 + X3 + X4, udist = "exponential",
#'   data = dairyspain, S = 1, method = "mla"
#' )
#'
#' vcov(cb_s_e)
#'
#' cb_s_r <- sfacross(
#'   formula = YIT ~ X1 + X2 + X3 + X4, udist = "rayleigh",
#'   data = dairyspain, S = 1, method = "nlminb"
#' )
#'
#' vcov(cb_s_r)
#' @keywords methods vcov
vcov.sfacross <- function(object, extraPar = FALSE, ...) {
  if (length(extraPar) != 1 || !is.logical(extraPar[1])) {
    stop("argument 'extraPar' must be a single logical value",
      call. = FALSE
    )
  }
  resCov <- object$invHessian
  if (extraPar) {
    if (object$udist %in% c("tnormal", "lognormal")) {
      delta <- object$mleParam[(object$nXvar + object$nmuHvar +
        1):(object$nXvar + object$nmuHvar + object$nuHvar)]
      phi <- object$mleParam[(object$nXvar + object$nmuHvar +
        object$nuHvar + 1):(object$nXvar + object$nmuHvar +
        object$nuHvar + object$nvHvar)]
      uHvar <- model.matrix(object$formula,
        data = object$dataTable,
        rhs = 3
      )
      vHvar <- model.matrix(object$formula,
        data = object$dataTable,
        rhs = 4
      )
    } else {
      delta <- object$mleParam[(object$nXvar + 1):(object$nXvar +
        object$nuHvar)]
      phi <- object$mleParam[(object$nXvar + object$nuHvar +
        1):(object$nXvar + object$nuHvar + object$nvHvar)]
      uHvar <- model.matrix(object$formula,
        data = object$dataTable,
        rhs = 2
      )
      vHvar <- model.matrix(object$formula,
        data = object$dataTable,
        rhs = 3
      )
    }
    Wu <- mean(as.numeric(crossprod(matrix(delta), t(uHvar))))
    Wv <- mean(as.numeric(crossprod(matrix(phi), t(vHvar))))
    if (object$nuHvar > 1 || object$nvHvar > 1 || object$nmuHvar > 1) {
      stop("argument 'extraPar' is not available for heteroscedasctic models",
        call. = FALSE
      )
    }
    jac <- diag(nrow(resCov))
    jac <- rbind(jac, matrix(0, nrow = 9, ncol = ncol(resCov)))
    rownames(jac) <- c(
      rownames(resCov), "sigmaSq", "lambdaSq",
      "sigmauSq", "sigmavSq", "sigma", "lambda", "sigmau",
      "sigmav", "gamma"
    )
    colnames(jac) <- colnames(resCov)
    jac["sigmaSq", "Zu_(Intercept)"] <- exp(Wu)
    jac["sigmaSq", "Zv_(Intercept)"] <- exp(Wv)
    jac["lambdaSq", "Zu_(Intercept)"] <- exp(Wu) / exp(Wv)
    jac["lambdaSq", "Zv_(Intercept)"] <- -exp(Wu + Wv) / exp(2 *
      Wv)
    jac["sigmauSq", "Zu_(Intercept)"] <- exp(Wu)
    jac["sigmavSq", "Zv_(Intercept)"] <- exp(Wv)
    jac["sigma", "Zu_(Intercept)"] <- 1 / 2 * exp(Wu) * (exp(Wu) +
      exp(Wv))^(-1 / 2)
    jac["sigma", "Zv_(Intercept)"] <- 1 / 2 * exp(Wv) * (exp(Wu) +
      exp(Wv))^(-1 / 2)
    jac["lambda", "Zu_(Intercept)"] <- 1 / 2 * exp(Wu / 2) / exp(Wv / 2)
    jac["lambda", "Zv_(Intercept)"] <- -1 / 2 * exp(Wu / 2 +
      Wv / 2) / exp(Wv)
    jac["sigmau", "Zu_(Intercept)"] <- 1 / 2 * exp(Wu / 2)
    jac["sigmav", "Zv_(Intercept)"] <- 1 / 2 * exp(Wv / 2)
    jac["gamma", "Zu_(Intercept)"] <- (exp(Wu) * (exp(Wu) +
      exp(Wv)) - exp(2 * Wu)) / (exp(Wu) + exp(Wv))^2
    jac["gamma", "Zv_(Intercept)"] <- -exp(Wu + Wv) / (exp(Wu) +
      exp(Wv))^2
    resCov <- jac %*% resCov %*% t(jac)
  }
  return(resCov)
}


# variance covariance matrix for lcmcross ----------
#' @rdname vcov
#' @aliases vcov.lcmcross
#' @export
vcov.lcmcross <- function(object, extraPar = FALSE, ...) {
  if (length(extraPar) != 1 || !is.logical(extraPar[1])) {
    stop("argument 'extraPar' must be a single logical value",
      call. = FALSE
    )
  }
  if (extraPar) {
    warning("'extraPar' is not yet available for object of class 'lcmcross'",
      call. = FALSE
    )
  }
  resCov <- object$invHessian
  return(resCov)
}
