# summary for sfacross ----------

#' @title Summary methods for frontier objects
#'
#' @name summary
#'
#' @aliases summary.sfacross
#'
#' @description Create and print summary results of objects of class
#'   \code{'sfacross'} and \code{'lcmcross'} (returned by \code{\link{sfacross}}
#'   and \code{\link{lcmcross}}).
#'
#' @param object An object of either class \code{"sfacross"}, returned by the
#'   function \code{\link{sfacross}}, or class \code{"lcmcross"}, returned by
#'   the function \code{\link{lcmcross}}.
#' @param x An object of either class \code{"summary.sfacross"} or
#'   \code{"summary.lcmcross"} returned by the function \code{\link{summary}}.
#' @param grad Logical. If \code{TRUE} the gradient for the maximum likelihood
#'   (ML) estimates of the different parameters is returned. Default =
#'   \code{FALSE}.
#' @param ci Logical. If \code{TRUE} the 95% confidence interval for the
#'   different parameters (OLS and ML estimates) is returned. Default =
#'   \code{FALSE}.
#' @param digits Numeric. Number of digits displayed in values.
#' @param ... Currently ignored.
#'
#' @details The summary methods creates objects of either class
#'   \code{"summary.sfacross"} or class \code{"summary.lcmcross"} that extend
#'   the objects of class \code{sfacross} or \code{lcmcross} with various
#'   information about the estimated model like information criterion, see
#'   \bold{Value}.
#'
#'   \code{coef} and \code{print} methods are associated (i.e.
#'   \code{coef.summary.sfacross}, \code{print.summary.sfacross} and
#'   \code{coef.summary.lcmcross}, \code{print.summary.lcmcross}).
#'
#' @return
#' The \code{summary} method returns a list of class \code{"summary.sfacross"}
#' or \code{"summary.lcmcross"} that is identical to an object returned by
#' \code{\link{sfacross}} or \code{\link{lcmcross}}, with the following
#' additional elements:
#'
#' \item{AIC}{Akaike information criterion.}
#'
#' \item{BIC}{Bayesian information criterion.}
#'
#' \item{HQIC}{Hannan-Quinn information criterion.}
#'
#' \item{sigmavSq}{Only for \code{"sfacross"}. Variance of the two-sided error
#' term (\eqn{\sigma_v^2}).}
#'
#' \item{sigmauSq}{Only for \code{"sfacross"}. \eqn{\sigma_u^2}.}
#'
#' \item{theta}{Only for \code{"sfacross"} and \code{"udist = uniform"}.
#' \code{theta} value in the case of the unifom distribution is defined as:
#' \eqn{u_i \in [0, \theta]}.}
#'
#' \item{Varu}{Only for \code{"sfacross"}. Variance of the one-sided error term.}
#'
#' \item{Eu}{Only for \code{"sfacross"}. Expected unconditional inefficiency.}
#'
#' \item{Expu}{Only for \code{"sfacross"}. Expected unconditional efficiency.}
#'
#' \item{olsRes}{Only for \code{"sfacross"}. Matrix of OLS estimates, their
#' standard errors, t-values, P-values, and when \code{ci = TRUE} their
#' confidence intervals.}
#'
#' \item{mleRes}{matrix of ML estimates, their standard errors, z-values,
#' asymptotic P-values and when \code{grad = TRUE} and or \code{ci = TRUE} their
#' gradient and or their confidence intervals.}
#'
#' \item{chisq}{Only for \code{"sfacross"}. Chi-square statistics of the
#' difference between the stochastic frontier and the OLS.}
#'
#' \item{df}{Degree of freedom for the inefficiency model.}
#'
#' \item{digits}{Number of digits displayed.}
#' @export
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
#' summary(tl_u_ts)
#'
#' ## Using data on eighty-two countries production (DGP)
#'
#' # LCM Cobb Douglas (production function) half normal distribution
#'
#' cb_2c_h <- lcmcross(formula = ly ~ lk + ll + yr, udist = "hnormal", data =
#' worldprod, S = 1, method = "ucminf")
#'
#' summary(cb_2c_h)
#' @keywords methods summary
summary.sfacross <- function(object, grad = FALSE, ci = FALSE,
                             ...) {
  if (length(grad) != 1 || !is.logical(grad[1])) {
    stop("argument 'grad' must be a single logical value",
      call. = FALSE
    )
  }
  if (length(ci) != 1 || !is.logical(ci[1])) {
    stop("argument 'ci' must be a single logical value",
      call. = FALSE
    )
  }
  object$AIC <- ic.sfacross(object, IC = "AIC")
  object$BIC <- ic.sfacross(object, IC = "BIC")
  object$HQIC <- ic.sfacross(object, IC = "HQIC")
  if (object$udist == "tnormal") {
    if (object$scaling) {
      delta <- object$mleParam[(object$nXvar + 1):(object$nXvar +
        (object$nuHvar - 1))]
      tau <- object$mleParam[object$nXvar + (object$nuHvar -
        1) + 1]
      cu <- object$mleParam[object$nXvar + (object$nuHvar -
        1) + 2]
      phi <- object$mleParam[(object$nXvar + (object$nuHvar -
        1) + 2 + 1):(object$nXvar + (object$nuHvar -
        1) + 2 + object$nvHvar)]
      muHvar <- model.matrix(object$formula,
        data = object$dataTable,
        rhs = 2
      )
      uHvar <- model.matrix(object$formula,
        data = object$dataTable,
        rhs = 3
      )
      vHvar <- model.matrix(object$formula,
        data = object$dataTable,
        rhs = 4
      )
    } else {
      omega <- object$mleParam[(object$nXvar + 1):(object$nXvar +
        object$nmuHvar)]
      delta <- object$mleParam[(object$nXvar + object$nmuHvar +
        1):(object$nXvar + object$nmuHvar + object$nuHvar)]
      phi <- object$mleParam[(object$nXvar + object$nmuHvar +
        object$nuHvar + 1):(object$nXvar + object$nmuHvar +
        object$nuHvar + object$nvHvar)]
      muHvar <- model.matrix(object$formula,
        data = object$dataTable,
        rhs = 2
      )
      uHvar <- model.matrix(object$formula,
        data = object$dataTable,
        rhs = 3
      )
      vHvar <- model.matrix(object$formula,
        data = object$dataTable,
        rhs = 4
      )
    }
  } else {
    if (object$udist == "lognormal") {
      omega <- object$mleParam[(object$nXvar + 1):(object$nXvar +
        object$nmuHvar)]
      delta <- object$mleParam[(object$nXvar + object$nmuHvar +
        1):(object$nXvar + object$nmuHvar + object$nuHvar)]
      phi <- object$mleParam[(object$nXvar + object$nmuHvar +
        object$nuHvar + 1):(object$nXvar + object$nmuHvar +
        object$nuHvar + object$nvHvar)]
      muHvar <- model.matrix(object$formula,
        data = object$dataTable,
        rhs = 2
      )
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
  }
  mu <- if (object$udist == "tnormal") {
    if (object$scaling) {
      exp(as.numeric(crossprod(matrix(delta), t(uHvar[
        ,
        -1
      ])))) * tau
    } else {
      as.numeric(crossprod(matrix(omega), t(muHvar)))
    }
  } else {
    if (object$udist == "lognormal") {
      as.numeric(crossprod(matrix(omega), t(muHvar)))
    } else {
      NULL
    }
  }
  P <- if (object$udist == "gamma") {
    object$mleParam[object$nXvar + object$nuHvar + object$nvHvar +
      1]
  } else {
    NULL
  }
  k <- if (object$udist == "weibull") {
    object$mleParam[object$nXvar + object$nuHvar + object$nvHvar +
      1]
  } else {
    NULL
  }
  lambda <- if (object$udist == "tslaplace") {
    object$mleParam[object$nXvar + object$nuHvar + object$nvHvar +
      1]
  } else {
    NULL
  }
  Wu <- if (object$udist == "tnormal" & object$scaling == TRUE) {
    cu + 2 * as.numeric(crossprod(matrix(delta), t(uHvar[
      ,
      -1
    ])))
  } else {
    as.numeric(crossprod(matrix(delta), t(uHvar)))
  }
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  object$sigmavSq <- mean(exp(Wv))
  object$sigmauSq <- mean(exp(Wu))
  if (object$udist == "uniform") {
    object$theta <- sqrt(12 * object$sigmauSq)
  }
  object$Varu <- varuFun(
    object = object, mu = mu, P = P, k = k,
    lambda = lambda
  )
  object$Eu <- euFun(
    object = object, mu = mu, P = P, k = k,
    lambda = lambda
  )
  object$Expu <- eExpuFun(
    object = object, mu = mu, P = P,
    k = k, lambda = lambda
  )
  # OLS estimates and stder, p-values, CI
  dfOLS <- object$Nobs - object$nXvar
  if (ci) {
    olsRes <- matrix(nrow = object$nXvar + 1, ncol = 6)
    colnames(olsRes) <- c(
      "Coefficient", "Std. Error", "binf",
      "bsup", "t value", "Pr(>|t|)"
    )
    olsRes[, 1] <- c(object$olsParam, object$olsSigmasq)
    olsRes[, 2] <- c(object$olsStder, NA)
    olsRes[, 3] <- olsRes[, 1] - qt(0.975, df = dfOLS) *
      olsRes[, 2]
    olsRes[, 4] <- olsRes[, 1] + qt(0.975, df = dfOLS) *
      olsRes[, 2]
    olsRes[, 5] <- olsRes[, 1] / olsRes[, 2]
    olsRes[, 6] <- 2 * pt(-abs(olsRes[, 5]), df = dfOLS)
  } else {
    olsRes <- matrix(nrow = object$nXvar + 1, ncol = 4)
    colnames(olsRes) <- c(
      "Coefficient", "Std. Error", "t value",
      "Pr(>|t|)"
    )
    olsRes[, 1] <- c(object$olsParam, object$olsSigmasq)
    olsRes[, 2] <- c(object$olsStder, NA)
    olsRes[, 3] <- olsRes[, 1] / olsRes[, 2]
    olsRes[, 4] <- 2 * pt(-abs(olsRes[, 3]), df = dfOLS)
  }
  row.names(olsRes) <- c(names(object$olsParam), "sigmaSq")
  object$olsRes <- olsRes
  # MLE estimates and stder, p-values, CI, Gradient
  if (grad && ci) {
    mleRes <- matrix(nrow = object$nParm, ncol = 7)
    colnames(mleRes) <- c(
      "Coefficient", "Std. Error", "binf",
      "bsup", "gradient", "z-value", "Pr(>|z|)"
    )
    mleRes[, 1] <- object$mleParam
    mleRes[, 2] <- sqrt(diag(object$invHessian))
    mleRes[, 3] <- mleRes[, 1] - qnorm(0.975) * mleRes[
      ,
      2
    ]
    mleRes[, 4] <- mleRes[, 1] + qnorm(0.975) * mleRes[
      ,
      2
    ]
    mleRes[, 5] <- object$gradient
    mleRes[, 6] <- mleRes[, 1] / mleRes[, 2]
    mleRes[, 7] <- 2 * pnorm(-abs(mleRes[, 6]))
  } else {
    if (grad == TRUE && ci == FALSE) {
      mleRes <- matrix(nrow = object$nParm, ncol = 5)
      colnames(mleRes) <- c(
        "Coefficient", "Std. Error",
        "gradient", "z-value", "Pr(>|z|)"
      )
      mleRes[, 1] <- object$mleParam
      mleRes[, 2] <- sqrt(diag(object$invHessian))
      mleRes[, 3] <- object$gradient
      mleRes[, 4] <- mleRes[, 1] / mleRes[, 2]
      mleRes[, 5] <- 2 * pnorm(-abs(mleRes[, 4]))
    } else {
      if (grad == FALSE && ci == TRUE) {
        mleRes <- matrix(nrow = object$nParm, ncol = 6)
        colnames(mleRes) <- c(
          "Coefficient", "Std. Error",
          "binf", "bsup", "z-value", "Pr(>|z|)"
        )
        mleRes[, 1] <- object$mleParam
        mleRes[, 2] <- sqrt(diag(object$invHessian))
        mleRes[, 3] <- mleRes[, 1] - qnorm(0.975) * mleRes[
          ,
          2
        ]
        mleRes[, 4] <- mleRes[, 1] + qnorm(0.975) * mleRes[
          ,
          2
        ]
        mleRes[, 5] <- mleRes[, 1] / mleRes[, 2]
        mleRes[, 6] <- 2 * pnorm(-abs(mleRes[, 5]))
      } else {
        mleRes <- matrix(nrow = object$nParm, ncol = 4)
        colnames(mleRes) <- c(
          "Coefficient", "Std. Error",
          "z value", "Pr(>|z|)"
        )
        mleRes[, 1] <- object$mleParam
        mleRes[, 2] <- sqrt(diag(object$invHessian))
        mleRes[, 3] <- mleRes[, 1] / mleRes[, 2]
        mleRes[, 4] <- 2 * pnorm(-abs(mleRes[, 3]))
      }
    }
  }
  row.names(mleRes) <- names(object$startVal)
  object$mleRes <- mleRes
  object$chisq <- 2 * (object$mleLoglik - object$olsLoglik)
  object$df <- object$nParm - object$nXvar - object$nvHvar
  class(object) <- "summary.sfacross"
  return(object)
}

# print summary for sfacross ----------
#' @rdname summary
#' @aliases print.summary.sfacross
#' @export

print.summary.sfacross <- function(x, digits = max(3, getOption("digits") - 2), ...) {
  mleRes <- x$mleRes
  if (dim(mleRes)[2] == 4) {
    mleRes[, 1] <- as.numeric(formatC(x$mleRes[, 1],
      digits = digits,
      format = "f"
    ))
    mleRes[, 2] <- as.numeric(formatC(x$mleRes[, 2],
      digits = digits,
      format = "f"
    ))
    mleRes[, 3] <- as.numeric(formatC(x$mleRes[, 3],
      digits = digits,
      format = "f"
    ))
    mleRes[, 4] <- as.numeric(formatC(x$mleRes[, 4],
      digits = digits,
      format = "e"
    ))
  } else {
    if (dim(mleRes)[2] == 5) {
      mleRes[, 1] <- as.numeric(formatC(x$mleRes[, 1],
        digits = digits, format = "f"
      ))
      mleRes[, 2] <- as.numeric(formatC(x$mleRes[, 2],
        digits = digits, format = "f"
      ))
      mleRes[, 3] <- as.numeric(formatC(x$mleRes[, 3],
        digits = digits, format = "e"
      ))
      mleRes[, 4] <- as.numeric(formatC(x$mleRes[, 4],
        digits = digits, format = "f"
      ))
      mleRes[, 5] <- as.numeric(formatC(x$mleRes[, 5],
        digits = digits, format = "e"
      ))
    } else {
      if (dim(mleRes)[2] == 6) {
        mleRes[, 1] <- as.numeric(formatC(x$mleRes[
          ,
          1
        ], digits = digits, format = "f"))
        mleRes[, 2] <- as.numeric(formatC(x$mleRes[
          ,
          2
        ], digits = digits, format = "f"))
        mleRes[, 3] <- as.numeric(formatC(x$mleRes[
          ,
          3
        ], digits = digits, format = "f"))
        mleRes[, 4] <- as.numeric(formatC(x$mleRes[
          ,
          4
        ], digits = digits, format = "f"))
        mleRes[, 5] <- as.numeric(formatC(x$mleRes[
          ,
          5
        ], digits = digits, format = "f"))
        mleRes[, 6] <- as.numeric(formatC(x$mleRes[
          ,
          6
        ], digits = digits, format = "e"))
      } else {
        if (dim(mleRes)[2] == 7) {
          mleRes[, 1] <- as.numeric(formatC(x$mleRes[
            ,
            1
          ], digits = digits, format = "f"))
          mleRes[, 2] <- as.numeric(formatC(x$mleRes[
            ,
            2
          ], digits = digits, format = "f"))
          mleRes[, 3] <- as.numeric(formatC(x$mleRes[
            ,
            3
          ], digits = digits, format = "f"))
          mleRes[, 4] <- as.numeric(formatC(x$mleRes[
            ,
            4
          ], digits = digits, format = "f"))
          mleRes[, 5] <- as.numeric(formatC(x$mleRes[
            ,
            5
          ], digits = digits, format = "e"))
          mleRes[, 6] <- as.numeric(formatC(x$mleRes[
            ,
            6
          ], digits = digits, format = "f"))
          mleRes[, 7] <- as.numeric(formatC(x$mleRes[
            ,
            7
          ], digits = digits, format = "e"))
        }
      }
    }
  }
  row.names(mleRes) <- formatC(row.names(mleRes),
    width = max(nchar(row.names(mleRes))),
    flag = "-"
  )
  mleRes1 <- mleRes[1:x$nXvar, , drop = FALSE]
  if (x$udist == "tnormal") {
    if (x$scaling) {
      mleRes2 <- mleRes[(x$nXvar + 1):(x$nXvar + (x$nuHvar -
        1)), , drop = FALSE]
      mleRes3 <- mleRes[x$nXvar + (x$nuHvar - 1) + 1, ,
        drop = FALSE
      ]
      mleRes4 <- mleRes[x$nXvar + (x$nuHvar - 1) + 2, ,
        drop = FALSE
      ]
      mleRes5 <- mleRes[(x$nXvar + (x$nuHvar - 1) + 2 +
        1):(x$nXvar + (x$nuHvar - 1) + 2 + x$nvHvar), ,
      drop = FALSE
      ]
    } else {
      mleRes2 <- mleRes[(x$nXvar + 1):(x$nXvar + x$nmuHvar), ,
        drop = FALSE
      ]
      mleRes3 <- mleRes[(x$nXvar + x$nmuHvar + 1):(x$nXvar +
        x$nmuHvar + x$nuHvar), , drop = FALSE]
      mleRes4 <- mleRes[(x$nXvar + x$nmuHvar + x$nuHvar +
        1):(x$nXvar + x$nmuHvar + x$nuHvar + x$nvHvar), ,
      drop = FALSE
      ]
    }
  } else {
    if (x$udist == "lognormal") {
      mleRes2 <- mleRes[(x$nXvar + 1):(x$nXvar + x$nmuHvar), ,
        drop = FALSE
      ]
      mleRes3 <- mleRes[(x$nXvar + x$nmuHvar + 1):(x$nXvar +
        x$nmuHvar + x$nuHvar), , drop = FALSE]
      mleRes4 <- mleRes[(x$nXvar + x$nmuHvar + x$nuHvar +
        1):(x$nXvar + x$nmuHvar + x$nuHvar + x$nvHvar), ,
      drop = FALSE
      ]
    } else {
      mleRes2 <- mleRes[(x$nXvar + 1):(x$nXvar + x$nuHvar), ,
        drop = FALSE
      ]
      mleRes3 <- mleRes[(x$nXvar + x$nuHvar + 1):(x$nXvar +
        x$nuHvar + x$nvHvar), , drop = FALSE]
      if (x$udist %in% c("gamma", "weibull", "tslaplace")) {
        mleRes4 <- mleRes[x$nXvar + x$nuHvar + x$nvHvar +
          1, , drop = FALSE]
      }
    }
  }
  lengthSum <- nchar(sfadist(x$udist)) + 10
  dimCoefTable <- as.character(dim(x$mleRes)[2])
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  cat(sfadist(x$udist), "\n")
  cat(
    "Dependent Variable:", paste0(rep(" ", lengthSum - nchar("Dependent Variable:") -
      nchar(paste0(attr(x$formula, "lhs")))), collapse = ""),
    paste0(attr(x$formula, "lhs")), "\n"
  )
  cat("Log likelihood solver:", paste0(rep(" ", lengthSum -
    nchar("Log likelihood solver:") - nchar(x$optType)),
  collapse = ""
  ), x$optType, "\n")
  cat("Log likelihood iter:", paste0(rep(" ", lengthSum - nchar("Log likelihood iter:") -
    nchar(x$nIter)), collapse = ""), x$nIter, "\n")
  cat("Log likelihood value:", paste0(rep(" ", lengthSum -
    nchar("Log likelihood value:") - nchar(formatC(x$mleLoglik,
      digits = digits, format = "f"
    ))), collapse = ""), formatC(x$mleLoglik,
    digits = digits, format = "f"
  ), "\n")
  cat("Log likelihood gradient norm:", paste0(rep(" ", lengthSum -
    nchar("Log likelihood gradient norm:") - nchar(formatC(x$gradientNorm,
      digits = digits, format = "e"
    ))), collapse = ""), formatC(x$gradientNorm,
    digits = digits, format = "e"
  ), "\n")
  cat(
    "Estimation based on:", paste0(rep(" ", lengthSum - nchar("Estimation based on:") -
      nchar(x$Nobs) - nchar(x$nParm) - nchar("N = ") - nchar("and K = ") -
      3), collapse = ""), "N = ", x$Nobs, "and K = ", x$nParm,
    "\n"
  )
  cat("Inf. Cr:", paste0(rep(" ", lengthSum - nchar("Inf. Cr:") -
    nchar("AIC  = ") - nchar(formatC(x$AIC, digits = 1, format = "f")) -
    nchar("AIC/N  = ") - nchar(formatC(x$AIC / x$Nobs,
      digits = 3,
      format = "f"
    )) - 3), collapse = ""), "AIC  = ", formatC(x$AIC,
    digits = 1, format = "f"
  ), "AIC/N  = ", formatC(x$AIC / x$Nobs,
    digits = 3, format = "f"
  ), "\n")
  cat(
    paste0(rep(" ", lengthSum - nchar("BIC  = ") - nchar(formatC(x$BIC,
      digits = 1, format = "f"
    )) - nchar("BIC/N  = ") - nchar(formatC(x$BIC / x$Nobs,
      digits = 3, format = "f"
    )) - 2), collapse = ""), "BIC  = ",
    formatC(x$BIC, digits = 1, format = "f"), "BIC/N  = ",
    formatC(x$BIC / x$Nobs, digits = 3, format = "f"), "\n"
  )
  cat(
    paste0(rep(" ", lengthSum - nchar("HQIC = ") - nchar(formatC(x$HQIC,
      digits = 1, format = "f"
    )) - nchar("HQIC/N = ") - nchar(formatC(x$HQIC / x$Nobs,
      digits = 3, format = "f"
    )) - 2), collapse = ""), "HQIC = ",
    formatC(x$HQIC, digits = 1, format = "f"), "HQIC/N = ",
    formatC(x$HQIC / x$Nobs, digits = 3, format = "f"), "\n"
  )
  cat(paste0(rep("-", lengthSum + 2), collapse = ""), "\n")
  cat("Variances: Sigma-squared(v)   = ", paste0(rep(" ", lengthSum -
    nchar("Variances: Sigma-squared(v)   = ") - nchar(formatC(x$sigmavSq,
      digits = digits, format = "f"
    ))), collapse = ""), formatC(x$sigmavSq,
    digits = digits, format = "f"
  ), "\n")
  cat("           Sigma(v)           = ", paste0(rep(" ", lengthSum -
    nchar("Variances: Sigma-squared(v)   = ") - nchar(formatC(sqrt(x$sigmavSq),
      digits = digits, format = "f"
    ))), collapse = ""), formatC(sqrt(x$sigmavSq),
    digits = digits, format = "f"
  ), "\n")
  cat("           Sigma-squared(u)   = ", paste0(rep(" ", lengthSum -
    nchar("Variances: Sigma-squared(u)   = ") - nchar(formatC(x$sigmauSq,
      digits = digits, format = "f"
    ))), collapse = ""), formatC(x$sigmauSq,
    digits = digits, format = "f"
  ), "\n")
  cat("           Sigma(u)           = ", paste0(rep(" ", lengthSum -
    nchar("Variances: Sigma-squared(u)   = ") - nchar(formatC(sqrt(x$sigmauSq),
      digits = digits, format = "f"
    ))), collapse = ""), formatC(sqrt(x$sigmauSq),
    digits = digits, format = "f"
  ), "\n")
  cat(
    "Sigma = Sqrt[(s^2(u)+s^2(v))] = ", paste0(rep(" ", lengthSum -
      nchar("Sigma = Sqrt[(s^2(u)+s^2(v))] = ") - nchar(formatC(sqrt(x$sigmavSq +
        x$sigmauSq), digits = digits, format = "f"))), collapse = ""),
    formatC(sqrt(x$sigmavSq + x$sigmauSq),
      digits = digits,
      format = "f"
    ), "\n"
  )
  cat(
    "Gamma = sigma(u)^2/sigma^2    = ", paste0(rep(" ", lengthSum -
      nchar("Gamma = sigma(u)^2/sigma^2    = ") - nchar(formatC(x$sigmauSq / (x$sigmavSq +
        x$sigmauSq), digits = digits, format = "f"))), collapse = ""),
    formatC(x$sigmauSq / (x$sigmavSq + x$sigmauSq),
      digits = digits,
      format = "f"
    ), "\n"
  )
  cat("Lambda = sigma(u)/sigma(v)    = ", paste0(rep(" ", lengthSum -
    nchar("Lambda = sigma(u)/sigma(v)    = ") - nchar(formatC(sqrt(x$sigmauSq / x$sigmavSq),
      digits = digits, format = "f"
    ))), collapse = ""), formatC(sqrt(x$sigmauSq / x$sigmavSq),
    digits = digits, format = "f"
  ), "\n")
  cat(
    "Var[u]/{Var[u]+Var[v]}        = ", paste0(rep(" ", lengthSum -
      nchar("Var[u]/{Var[u]+Var[v]}        = ") - nchar(formatC(x$Varu / (x$Varu +
        x$sigmavSq), digits = digits, format = "f"))), collapse = ""),
    formatC(x$Varu / (x$Varu + x$sigmavSq),
      digits = digits,
      format = "f"
    ), "\n"
  )
  cat(
    "Var[e]                        = ", paste0(rep(" ", lengthSum -
      nchar("Var[e]                        = ") - nchar(formatC(x$Varu +
        x$sigmavSq, digits = digits, format = "f"))), collapse = ""),
    formatC(x$Varu + x$sigmavSq, digits = digits, format = "f"),
    "\n"
  )
  if (x$udist == "uniform") {
    cat("theta                         = ", paste0(rep(
      " ",
      lengthSum - nchar("theta                         = ") -
        nchar(formatC(x$theta, digits = digits, format = "f"))
    ),
    collapse = ""
    ), formatC(x$theta,
      digits = digits,
      format = "f"
    ), "\n")
  }
  if (x$nuHvar > 1 || x$nvHvar > 1) {
    cat("Variances averaged over observations \n")
  }
  cat(paste0(rep("-", lengthSum + 2), collapse = ""), "\n")
  cat("Average inefficiency E[u]     = ", paste0(rep(" ", lengthSum -
    nchar("Average inefficiency E[u]     = ") - nchar(formatC(x$Eu,
      digits = digits, format = "f"
    ))), collapse = ""), formatC(x$Eu,
    digits = digits, format = "f"
  ), "\n")
  cat("Average efficiency E[exp(-u)] = ", paste0(rep(" ", lengthSum -
    nchar("Average efficiency E[exp(-u)] = ") - nchar(formatC(x$Expu,
      digits = digits, format = "f"
    ))), collapse = ""), formatC(x$Expu,
    digits = digits, format = "f"
  ), "\n")
  cat(paste0(rep("-", lengthSum + 2), collapse = ""), "\n")
  cat(x$typeSfa, "\n")
  cat("-----[ Tests vs. No Inefficiency ]-----\n")
  cat("Likelihood Ratio Test of Inefficiency\n")
  cat(
    "Deg. freedom for inefficiency model", paste0(rep(
      " ",
      lengthSum - nchar("Deg. freedom for inefficiency model") -
        nchar(formatC(x$df, digits = digits, format = "d"))
    ),
    collapse = ""
    ), formatC(x$df, digits = digits, format = "d"),
    "\n"
  )
  cat("Log Likelihood for OLS Log(H0) = ", paste0(rep(
    " ",
    lengthSum - nchar("Log Likelihood for OLS Log(H0) = ") -
      nchar(formatC(x$olsLoglik, digits = digits, format = "f"))
  ),
  collapse = ""
  ), formatC(x$olsLoglik,
    digits = digits,
    format = "f"
  ), "\n")
  cat("LR statistic: ", "\n")
  cat(
    "Chisq = 2*[LogL(H0)-LogL(H1)]  = ", paste0(rep(
      " ",
      lengthSum - nchar("Chisq = 2*[LogL(H0)-LogL(H1)]  = ") -
        nchar(formatC(x$chisq, digits = digits, format = "f"))
    ),
    collapse = ""
    ), formatC(x$chisq, digits = digits, format = "f"),
    "\n"
  )
  cat("Kodde-Palm C*:       95%:", formatC(qchibarsq(0.95,
    df = x$df
  ), digits = digits, format = "f"), paste0(rep(
    " ",
    lengthSum - nchar("Kodde-Palm C*:       95%:") - nchar(formatC(qchibarsq(0.95,
      df = x$df
    ), digits = digits, format = "f")) - nchar(formatC(qchibarsq(0.99,
      df = x$df
    ), digits = digits, format = "f")) - nchar("99%") -
      3
  ), collapse = ""), "99%:", formatC(qchibarsq(0.99,
    df = x$df
  ), digits = digits, format = "f"), "\n")
  cat("Coelli (1995) skewness test on OLS residuals\n")
  cat("M3T                            = ", paste0(rep(
    " ",
    lengthSum - nchar("M3T                            = ") -
      nchar(formatC(x$CoelliM3Test[1],
        digits = digits,
        format = "f"
      ))
  ), collapse = ""), formatC(x$CoelliM3Test[1],
    digits = digits, format = "f"
  ), "\n")
  cat("final maximum likelihood estimates \n")
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  cat(centerText("Deterministic Component of SFA", width = lengthSum +
    2 + switch(dimCoefTable, `4` = 18, `5` = 31, `6` = 43,
      `7` = 57
    )), "\n")
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  printCoefmat(mleRes1, P.values = TRUE, digits = digits, signif.legend = FALSE)
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  if (x$udist == "tnormal") {
    if (x$scaling) {
      cat(centerText("Scaling property parameters", width = lengthSum +
        2 + switch(dimCoefTable, `4` = 18, `5` = 31,
          `6` = 43, `7` = 57
        )), "\n")
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      printCoefmat(mleRes2,
        P.values = TRUE, digits = digits,
        signif.legend = FALSE
      )
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      cat(centerText("Location parameter [offset mu] in u (one sided error)",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )
      ), "\n")
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      printCoefmat(mleRes3,
        P.values = TRUE, digits = digits,
        signif.legend = FALSE
      )
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      cat(centerText("Parameter in variance of u (one sided error)",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )
      ), "\n")
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      printCoefmat(mleRes4,
        P.values = TRUE, digits = digits,
        signif.legend = FALSE
      )
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      cat(centerText("Parameters in variance of v (symmetric error)",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )
      ), "\n")
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      printCoefmat(mleRes5,
        P.values = TRUE, digits = digits,
        signif.legend = TRUE
      )
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
    } else {
      cat(centerText("Location parameter [offset mu] in u (one sided error)",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )
      ), "\n")
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      printCoefmat(mleRes2,
        P.values = TRUE, digits = digits,
        signif.legend = FALSE
      )
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      cat(centerText("Parameter in variance of u (one sided error)",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )
      ), "\n")
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      printCoefmat(mleRes3,
        P.values = TRUE, digits = digits,
        signif.legend = FALSE
      )
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      cat(centerText("Parameters in variance of v (symmetric error)",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )
      ), "\n")
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      printCoefmat(mleRes4,
        P.values = TRUE, digits = digits,
        signif.legend = TRUE
      )
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
    }
  } else {
    if (x$udist == "lognormal") {
      cat(centerText("Location parameter [offset mu] in u (one sided error)",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )
      ), "\n")
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      printCoefmat(mleRes2,
        P.values = TRUE, digits = digits,
        signif.legend = FALSE
      )
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      cat(centerText("Parameter in variance of u (one sided error)",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )
      ), "\n")
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      printCoefmat(mleRes3,
        P.values = TRUE, digits = digits,
        signif.legend = FALSE
      )
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      cat(centerText("Parameters in variance of v (symmetric error)",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )
      ), "\n")
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      printCoefmat(mleRes4,
        P.values = TRUE, digits = digits,
        signif.legend = TRUE
      )
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
    } else {
      if (x$udist == "gamma") {
        cat(
          centerText("Parameter in variance of u (one sided error)",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )
          ),
          "\n"
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        printCoefmat(mleRes2,
          P.values = TRUE, digits = digits,
          signif.legend = FALSE
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        cat(
          centerText("Parameters in variance of v (symmetric error)",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )
          ),
          "\n"
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        printCoefmat(mleRes3,
          P.values = TRUE, digits = digits,
          signif.legend = TRUE
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        cat(
          centerText("Location parameter P in u (one sided error)",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )
          ),
          "\n"
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        printCoefmat(mleRes4,
          P.values = TRUE, digits = digits,
          signif.legend = TRUE
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
      } else {
        if (x$udist == "weibull") {
          cat(
            centerText("Parameter in variance of u (one sided error)",
              width = lengthSum + 2 + switch(dimCoefTable,
                `4` = 18, `5` = 31, `6` = 43, `7` = 57
              )
            ),
            "\n"
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          printCoefmat(mleRes2,
            P.values = TRUE, digits = digits,
            signif.legend = FALSE
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          cat(
            centerText("Parameters in variance of v (symmetric error)",
              width = lengthSum + 2 + switch(dimCoefTable,
                `4` = 18, `5` = 31, `6` = 43, `7` = 57
              )
            ),
            "\n"
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          printCoefmat(mleRes3,
            P.values = TRUE, digits = digits,
            signif.legend = TRUE
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          cat(
            centerText("Shape parameter k in u (one sided error)",
              width = lengthSum + 2 + switch(dimCoefTable,
                `4` = 18, `5` = 31, `6` = 43, `7` = 57
              )
            ),
            "\n"
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          printCoefmat(mleRes4,
            P.values = TRUE, digits = digits,
            signif.legend = TRUE
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
        } else {
          if (x$udist == "tslaplace") {
            cat(
              centerText("Parameter in variance of u (one sided error)",
                width = lengthSum + 2 + switch(dimCoefTable,
                  `4` = 18, `5` = 31, `6` = 43, `7` = 57
                )
              ),
              "\n"
            )
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )),
            collapse = ""
            ), "\n")
            printCoefmat(mleRes2,
              P.values = TRUE, digits = digits,
              signif.legend = FALSE
            )
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )),
            collapse = ""
            ), "\n")
            cat(
              centerText("Parameters in variance of v (symmetric error)",
                width = lengthSum + 2 + switch(dimCoefTable,
                  `4` = 18, `5` = 31, `6` = 43, `7` = 57
                )
              ),
              "\n"
            )
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )),
            collapse = ""
            ), "\n")
            printCoefmat(mleRes3,
              P.values = TRUE, digits = digits,
              signif.legend = TRUE
            )
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )),
            collapse = ""
            ), "\n")
            cat(
              centerText("Skewness parameter 'lambda' in u (one sided error)",
                width = lengthSum + 2 + switch(dimCoefTable,
                  `4` = 18, `5` = 31, `6` = 43, `7` = 57
                )
              ),
              "\n"
            )
            printCoefmat(mleRes4,
              P.values = TRUE, digits = digits,
              signif.legend = TRUE
            )
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )),
            collapse = ""
            ), "\n")
          } else {
            cat(
              centerText("Parameter in variance of u (one sided error)",
                width = lengthSum + 2 + switch(dimCoefTable,
                  `4` = 18, `5` = 31, `6` = 43, `7` = 57
                )
              ),
              "\n"
            )
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )),
            collapse = ""
            ), "\n")
            printCoefmat(mleRes2,
              P.values = TRUE, digits = digits,
              signif.legend = FALSE
            )
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )),
            collapse = ""
            ), "\n")
            cat(
              centerText("Parameters in variance of v (symmetric error)",
                width = lengthSum + 2 + switch(dimCoefTable,
                  `4` = 18, `5` = 31, `6` = 43, `7` = 57
                )
              ),
              "\n"
            )
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )),
            collapse = ""
            ), "\n")
            printCoefmat(mleRes3,
              P.values = TRUE, digits = digits,
              signif.legend = TRUE
            )
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )),
            collapse = ""
            ), "\n")
          }
        }
      }
    }
  }
  cat(x$mleDate, "\n")
  cat("Log likelihood status:", x$optStatus, "\n")
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  invisible(x)
}

# summary for lcmcross ----------
#' @rdname summary
#' @aliases summary.lcmcross
#' @export
summary.lcmcross <- function(object, grad = FALSE, ci = FALSE,
                             ...) {
  if (length(grad) != 1 || !is.logical(grad[1])) {
    stop("argument 'grad' must be a single logical value",
      call. = FALSE
    )
  }
  if (length(ci) != 1 || !is.logical(ci[1])) {
    stop("argument 'ci' must be a single logical value",
      call. = FALSE
    )
  }
  object$AIC <- ic.lcmcross(object, IC = "AIC")
  object$BIC <- ic.lcmcross(object, IC = "BIC")
  object$HQIC <- ic.lcmcross(object, IC = "HQIC")
  # MLE estimates and stder, p-values, CI, Gradient
  if (grad && ci) {
    mleRes <- matrix(nrow = object$nParm, ncol = 7)
    colnames(mleRes) <- c(
      "Coefficient", "Std. Error", "binf",
      "bsup", "gradient", "z-value", "Pr(>|z|)"
    )
    mleRes[, 1] <- object$mleParam
    mleRes[, 2] <- sqrt(diag(object$invHessian))
    mleRes[, 3] <- mleRes[, 1] - qnorm(0.975) * mleRes[
      ,
      2
    ]
    mleRes[, 4] <- mleRes[, 1] + qnorm(0.975) * mleRes[
      ,
      2
    ]
    mleRes[, 5] <- object$gradient
    mleRes[, 6] <- mleRes[, 1] / mleRes[, 2]
    mleRes[, 7] <- 2 * pnorm(-abs(mleRes[, 6]))
  } else {
    if (grad == TRUE && ci == FALSE) {
      mleRes <- matrix(nrow = object$nParm, ncol = 5)
      colnames(mleRes) <- c(
        "Coefficient", "Std. Error",
        "gradient", "z-value", "Pr(>|z|)"
      )
      mleRes[, 1] <- object$mleParam
      mleRes[, 2] <- sqrt(diag(object$invHessian))
      mleRes[, 3] <- object$gradient
      mleRes[, 4] <- mleRes[, 1] / mleRes[, 2]
      mleRes[, 5] <- 2 * pnorm(-abs(mleRes[, 4]))
    } else {
      if (grad == FALSE && ci == TRUE) {
        mleRes <- matrix(nrow = object$nParm, ncol = 6)
        colnames(mleRes) <- c(
          "Coefficient", "Std. Error",
          "binf", "bsup", "z-value", "Pr(>|z|)"
        )
        mleRes[, 1] <- object$mleParam
        mleRes[, 2] <- sqrt(diag(object$invHessian))
        mleRes[, 3] <- mleRes[, 1] - qnorm(0.975) * mleRes[
          ,
          2
        ]
        mleRes[, 4] <- mleRes[, 1] + qnorm(0.975) * mleRes[
          ,
          2
        ]
        mleRes[, 5] <- mleRes[, 1] / mleRes[, 2]
        mleRes[, 6] <- 2 * pnorm(-abs(mleRes[, 5]))
      } else {
        mleRes <- matrix(nrow = object$nParm, ncol = 4)
        colnames(mleRes) <- c(
          "Coefficient", "Std. Error",
          "z value", "Pr(>|z|)"
        )
        mleRes[, 1] <- object$mleParam
        mleRes[, 2] <- sqrt(diag(object$invHessian))
        mleRes[, 3] <- mleRes[, 1] / mleRes[, 2]
        mleRes[, 4] <- 2 * pnorm(-abs(mleRes[, 3]))
      }
    }
  }
  row.names(mleRes) <- names(object$startVal)
  object$mleRes <- mleRes
  # object$chisq <- 2 * (object$mleLoglik - object$olsLoglik)
  object$df <- object$nParm - object$nClasses * object$nXvar -
    object$nClasses * object$nvHvar - object$nZvar *
    (object$nClasses - 1)
  class(object) <- "summary.lcmcross"
  return(object)
}


# print summary for lcmcross ----------
#' @rdname summary
#' @aliases print.summary.lcmcross
#' @export
print.summary.lcmcross <- function(x, digits = max(3, getOption("digits") - 2), ...) {
  mleRes <- x$mleRes
  if (dim(mleRes)[2] == 4) {
    mleRes[, 1] <- as.numeric(formatC(x$mleRes[, 1],
      digits = digits,
      format = "f"
    ))
    mleRes[, 2] <- as.numeric(formatC(x$mleRes[, 2],
      digits = digits,
      format = "f"
    ))
    mleRes[, 3] <- as.numeric(formatC(x$mleRes[, 3],
      digits = digits,
      format = "f"
    ))
    mleRes[, 4] <- as.numeric(formatC(x$mleRes[, 4],
      digits = digits,
      format = "e"
    ))
  } else {
    if (dim(mleRes)[2] == 5) {
      mleRes[, 1] <- as.numeric(formatC(x$mleRes[, 1],
        digits = digits, format = "f"
      ))
      mleRes[, 2] <- as.numeric(formatC(x$mleRes[, 2],
        digits = digits, format = "f"
      ))
      mleRes[, 3] <- as.numeric(formatC(x$mleRes[, 3],
        digits = digits, format = "e"
      ))
      mleRes[, 4] <- as.numeric(formatC(x$mleRes[, 4],
        digits = digits, format = "f"
      ))
      mleRes[, 5] <- as.numeric(formatC(x$mleRes[, 5],
        digits = digits, format = "e"
      ))
    } else {
      if (dim(mleRes)[2] == 6) {
        mleRes[, 1] <- as.numeric(formatC(x$mleRes[
          ,
          1
        ], digits = digits, format = "f"))
        mleRes[, 2] <- as.numeric(formatC(x$mleRes[
          ,
          2
        ], digits = digits, format = "f"))
        mleRes[, 3] <- as.numeric(formatC(x$mleRes[
          ,
          3
        ], digits = digits, format = "f"))
        mleRes[, 4] <- as.numeric(formatC(x$mleRes[
          ,
          4
        ], digits = digits, format = "f"))
        mleRes[, 5] <- as.numeric(formatC(x$mleRes[
          ,
          5
        ], digits = digits, format = "f"))
        mleRes[, 6] <- as.numeric(formatC(x$mleRes[
          ,
          6
        ], digits = digits, format = "e"))
      } else {
        if (dim(mleRes)[2] == 7) {
          mleRes[, 1] <- as.numeric(formatC(x$mleRes[
            ,
            1
          ], digits = digits, format = "f"))
          mleRes[, 2] <- as.numeric(formatC(x$mleRes[
            ,
            2
          ], digits = digits, format = "f"))
          mleRes[, 3] <- as.numeric(formatC(x$mleRes[
            ,
            3
          ], digits = digits, format = "f"))
          mleRes[, 4] <- as.numeric(formatC(x$mleRes[
            ,
            4
          ], digits = digits, format = "f"))
          mleRes[, 5] <- as.numeric(formatC(x$mleRes[
            ,
            5
          ], digits = digits, format = "e"))
          mleRes[, 6] <- as.numeric(formatC(x$mleRes[
            ,
            6
          ], digits = digits, format = "f"))
          mleRes[, 7] <- as.numeric(formatC(x$mleRes[
            ,
            7
          ], digits = digits, format = "e"))
        }
      }
    }
  }
  row.names(mleRes) <- formatC(row.names(mleRes),
    width = max(nchar(row.names(mleRes))),
    flag = "-"
  )
  mleRes1 <- mleRes[1:x$nXvar, ]
  mleRes2 <- mleRes[(x$nXvar + 1):(x$nXvar + x$nuHvar), , drop = FALSE]
  mleRes3 <- mleRes[(x$nXvar + x$nuHvar + 1):(x$nXvar + x$nuHvar +
    x$nvHvar), , drop = FALSE]
  sfaModel <- "Normal-Half Normal Latent Class Stochastic Frontier Model"
  lengthSum <- nchar(sfaModel) # + 10
  dimCoefTable <- as.character(dim(x$mleRes)[2])
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  cat(sfaModel, "\n")
  cat(
    "Dependent Variable:", paste0(rep(" ", lengthSum - nchar("Dependent Variable:") -
      nchar(paste0(attr(x$formula, "lhs")))), collapse = ""),
    paste0(attr(x$formula, "lhs")), "\n"
  )
  cat("Log likelihood solver:", paste0(rep(" ", lengthSum -
    nchar("Log likelihood solver:") - nchar(x$optType)),
  collapse = ""
  ), x$optType, "\n")
  cat("Log likelihood iter:", paste0(rep(" ", lengthSum - nchar("Log likelihood iter:") -
    nchar(x$nIter)), collapse = ""), x$nIter, "\n")
  cat("Log likelihood value:", paste0(rep(" ", lengthSum -
    nchar("Log likelihood value:") - nchar(formatC(x$mleLoglik,
      digits = digits, format = "f"
    ))), collapse = ""), formatC(x$mleLoglik,
    digits = digits, format = "f"
  ), "\n")
  cat("Log likelihood gradient norm:", paste0(rep(" ", lengthSum -
    nchar("Log likelihood gradient norm:") - nchar(formatC(x$gradientNorm,
      digits = digits, format = "e"
    ))), collapse = ""), formatC(x$gradientNorm,
    digits = digits, format = "e"
  ), "\n")
  cat(
    "Estimation based on:", paste0(rep(" ", lengthSum - nchar("Estimation based on:") -
      nchar(x$Nobs) - nchar(x$nParm) - nchar("N = ") - nchar("and K = ") -
      3), collapse = ""), "N = ", x$Nobs, "and K = ", x$nParm,
    "\n"
  )
  cat("Inf. Cr:", paste0(rep(" ", lengthSum - nchar("Inf. Cr:") -
    nchar("AIC  = ") - nchar(formatC(x$AIC, digits = 1, format = "f")) -
    nchar("AIC/N  = ") - nchar(formatC(x$AIC / x$Nobs,
      digits = 3,
      format = "f"
    )) - 3), collapse = ""), "AIC  = ", formatC(x$AIC,
    digits = 1, format = "f"
  ), "AIC/N  = ", formatC(x$AIC / x$Nobs,
    digits = 3, format = "f"
  ), "\n")
  cat(
    paste0(rep(" ", lengthSum - nchar("BIC  = ") - nchar(formatC(x$BIC,
      digits = 1, format = "f"
    )) - nchar("BIC/N  = ") - nchar(formatC(x$BIC / x$Nobs,
      digits = 3, format = "f"
    )) - 2), collapse = ""), "BIC  = ",
    formatC(x$BIC, digits = 1, format = "f"), "BIC/N  = ",
    formatC(x$BIC / x$Nobs, digits = 3, format = "f"), "\n"
  )
  cat(
    paste0(rep(" ", lengthSum - nchar("HQIC = ") - nchar(formatC(x$HQIC,
      digits = 1, format = "f"
    )) - nchar("HQIC/N = ") - nchar(formatC(x$HQIC / x$Nobs,
      digits = 3, format = "f"
    )) - 2), collapse = ""), "HQIC = ",
    formatC(x$HQIC, digits = 1, format = "f"), "HQIC/N = ",
    formatC(x$HQIC / x$Nobs, digits = 3, format = "f"), "\n"
  )
  cat(paste0(rep("-", lengthSum + 2), collapse = ""), "\n")
  cat(x$typeSfa, "\n")
  cat("Latent class model with", x$nClasses, "latent classes \n")
  cat("final maximum likelihood estimates \n")
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  cat(centerText("Deterministic Component of SFA for latent class 1",
    width = lengthSum + 2 + switch(dimCoefTable, `4` = 18,
      `5` = 31, `6` = 43, `7` = 57
    )
  ), "\n")
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  printCoefmat(mleRes[1:x$nXvar, , drop = FALSE],
    P.values = TRUE,
    digits = digits, signif.legend = FALSE
  )
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  cat(centerText("Parameter in variance of u (one sided error) for latent class 1",
    width = lengthSum + 2 + switch(dimCoefTable, `4` = 18,
      `5` = 31, `6` = 43, `7` = 57
    )
  ), "\n")
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  printCoefmat(mleRes[(x$nXvar + 1):(x$nXvar + x$nuHvar), ,
    drop = FALSE
  ], P.values = TRUE, digits = digits, signif.legend = FALSE)
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  cat(centerText("Parameters in variance of v (symmetric error) for latent class 1",
    width = lengthSum + 2 + switch(dimCoefTable, `4` = 18,
      `5` = 31, `6` = 43, `7` = 57
    )
  ), "\n")
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  printCoefmat(mleRes[(x$nXvar + x$nuHvar + 1):(x$nXvar + x$nuHvar +
    x$nvHvar), , drop = FALSE],
  P.values = TRUE, digits = digits,
  signif.legend = FALSE
  )
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  cat(centerText("Deterministic Component of SFA for latent class 2",
    width = lengthSum + 2 + switch(dimCoefTable, `4` = 18,
      `5` = 31, `6` = 43, `7` = 57
    )
  ), "\n")
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  printCoefmat(mleRes[(x$nXvar + x$nuHvar + x$nvHvar + 1):(2 *
    x$nXvar + x$nuHvar + x$nvHvar), , drop = FALSE],
  P.values = TRUE,
  digits = digits, signif.legend = FALSE
  )
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  cat(centerText("Parameter in variance of u (one sided error) for latent class 2",
    width = lengthSum + 2 + switch(dimCoefTable, `4` = 18,
      `5` = 31, `6` = 43, `7` = 57
    )
  ), "\n")
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  printCoefmat(mleRes[(2 * x$nXvar + x$nuHvar + x$nvHvar +
    1):(2 * x$nXvar + 2 * x$nuHvar + x$nvHvar), , drop = FALSE],
  P.values = TRUE, digits = digits, signif.legend = FALSE
  )
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  cat(centerText("Parameters in variance of v (symmetric error) for latent class 2",
    width = lengthSum + 2 + switch(dimCoefTable, `4` = 18,
      `5` = 31, `6` = 43, `7` = 57
    )
  ), "\n")
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  printCoefmat(mleRes[(2 * x$nXvar + 2 * x$nuHvar + x$nvHvar +
    1):(2 * x$nXvar + 2 * x$nuHvar + 2 * x$nvHvar), , drop = FALSE],
  P.values = TRUE, digits = digits, signif.legend = FALSE
  )
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  if (x$nClasses == 2) {
    cat(centerText("Estimated prior probabilities for class membership",
      width = lengthSum + 2 + switch(dimCoefTable, `4` = 18,
        `5` = 31, `6` = 43, `7` = 57
      )
    ), "\n")
    cat(
      paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57
      )), collapse = ""),
      "\n"
    )
    printCoefmat(mleRes[(2 * x$nXvar + 2 * x$nuHvar + 2 *
      x$nvHvar + 1):(2 * x$nXvar + 2 * x$nuHvar + 2 * x$nvHvar +
      x$nZvar), , drop = FALSE],
    P.values = TRUE, digits = digits,
    signif.legend = TRUE
    )
    cat(
      paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57
      )), collapse = ""),
      "\n"
    )
  } else {
    if (x$nClasses == 3) {
      cat(centerText("Deterministic Component of SFA for latent class 3",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )
      ), "\n")
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      printCoefmat(mleRes[(2 * x$nXvar + 2 * x$nuHvar +
        2 * x$nvHvar + 1):(3 * x$nXvar + 2 * x$nuHvar +
        2 * x$nvHvar), , drop = FALSE],
      P.values = TRUE,
      digits = digits, signif.legend = FALSE
      )
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      cat(centerText("Parameter in variance of u (one sided error) for latent class 3",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )
      ), "\n")
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      printCoefmat(mleRes[(3 * x$nXvar + 2 * x$nuHvar +
        2 * x$nvHvar + 1):(3 * x$nXvar + 3 * x$nuHvar +
        2 * x$nvHvar), , drop = FALSE],
      P.values = TRUE,
      digits = digits, signif.legend = FALSE
      )
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      cat(centerText("Parameters in variance of v (symmetric error) for latent class 3",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )
      ), "\n")
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      printCoefmat(mleRes[(3 * x$nXvar + 3 * x$nuHvar +
        2 * x$nvHvar + 1):(3 * x$nXvar + 3 * x$nuHvar +
        3 * x$nvHvar), , drop = FALSE],
      P.values = TRUE,
      digits = digits, signif.legend = FALSE
      )
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      cat(centerText("Estimated prior probabilities for class membership",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )
      ), "\n")
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      printCoefmat(mleRes[(3 * x$nXvar + 3 * x$nuHvar +
        3 * x$nvHvar + 1):(3 * x$nXvar + 3 * x$nuHvar +
        3 * x$nvHvar + 2 * x$nZvar), , drop = FALSE],
      P.values = TRUE, digits = digits, signif.legend = TRUE
      )
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
    } else {
      if (x$nClasses == 4) {
        cat(
          centerText("Deterministic Component of SFA for latent class 3",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )
          ),
          "\n"
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        printCoefmat(mleRes[(2 * x$nXvar + 2 * x$nuHvar +
          2 * x$nvHvar + 1):(3 * x$nXvar + 2 * x$nuHvar +
          2 * x$nvHvar), , drop = FALSE],
        P.values = TRUE,
        digits = digits, signif.legend = FALSE
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        cat(
          centerText("Parameter in variance of u (one sided error) for latent class 3",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )
          ),
          "\n"
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        printCoefmat(mleRes[(3 * x$nXvar + 2 * x$nuHvar +
          2 * x$nvHvar + 1):(3 * x$nXvar + 3 * x$nuHvar +
          2 * x$nvHvar), , drop = FALSE],
        P.values = TRUE,
        digits = digits, signif.legend = FALSE
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        cat(
          centerText("Parameters in variance of v (symmetric error) for latent class 3",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )
          ),
          "\n"
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        printCoefmat(mleRes[(3 * x$nXvar + 3 * x$nuHvar +
          2 * x$nvHvar + 1):(3 * x$nXvar + 3 * x$nuHvar +
          3 * x$nvHvar), , drop = FALSE],
        P.values = TRUE,
        digits = digits, signif.legend = FALSE
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        cat(
          centerText("Deterministic Component of SFA for latent class 4",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )
          ),
          "\n"
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        printCoefmat(mleRes[(3 * x$nXvar + 3 * x$nuHvar +
          3 * x$nvHvar + 1):(4 * x$nXvar + 3 * x$nuHvar +
          3 * x$nvHvar), , drop = FALSE],
        P.values = TRUE,
        digits = digits, signif.legend = FALSE
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        cat(
          centerText("Parameter in variance of u (one sided error) for latent class 4",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )
          ),
          "\n"
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        printCoefmat(mleRes[(4 * x$nXvar + 3 * x$nuHvar +
          3 * x$nvHvar + 1):(4 * x$nXvar + 4 * x$nuHvar +
          3 * x$nvHvar), , drop = FALSE],
        P.values = TRUE,
        digits = digits, signif.legend = FALSE
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        cat(
          centerText("Parameters in variance of v (symmetric error) for latent class 4",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )
          ),
          "\n"
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        printCoefmat(mleRes[(4 * x$nXvar + 4 * x$nuHvar +
          3 * x$nvHvar + 1):(4 * x$nXvar + 4 * x$nuHvar +
          4 * x$nvHvar), , drop = FALSE],
        P.values = TRUE,
        digits = digits, signif.legend = FALSE
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        cat(
          centerText("Estimated prior probabilities for class membership",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )
          ),
          "\n"
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        printCoefmat(mleRes[(4 * x$nXvar + 4 * x$nuHvar +
          4 * x$nvHvar + 1):(4 * x$nXvar + 4 * x$nuHvar +
          4 * x$nvHvar + 3 * x$nZvar), , drop = FALSE],
        P.values = TRUE, digits = digits, signif.legend = TRUE
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
      } else {
        if (x$nClasses == 5) {
          cat(
            centerText("Deterministic Component of SFA for latent class 3",
              width = lengthSum + 2 + switch(dimCoefTable,
                `4` = 18, `5` = 31, `6` = 43, `7` = 57
              )
            ),
            "\n"
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          printCoefmat(mleRes[(2 * x$nXvar + 2 * x$nuHvar +
            2 * x$nvHvar + 1):(3 * x$nXvar + 2 * x$nuHvar +
            2 * x$nvHvar), , drop = FALSE],
          P.values = TRUE,
          digits = digits, signif.legend = FALSE
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          cat(
            centerText("Parameter in variance of u (one sided error) for latent class 3",
              width = lengthSum + 2 + switch(dimCoefTable,
                `4` = 18, `5` = 31, `6` = 43, `7` = 57
              )
            ),
            "\n"
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          printCoefmat(mleRes[(3 * x$nXvar + 2 * x$nuHvar +
            2 * x$nvHvar + 1):(3 * x$nXvar + 3 * x$nuHvar +
            2 * x$nvHvar), , drop = FALSE],
          P.values = TRUE,
          digits = digits, signif.legend = FALSE
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          cat(
            centerText("Parameters in variance of v (symmetric error) for latent class 3",
              width = lengthSum + 2 + switch(dimCoefTable,
                `4` = 18, `5` = 31, `6` = 43, `7` = 57
              )
            ),
            "\n"
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          printCoefmat(mleRes[(3 * x$nXvar + 3 * x$nuHvar +
            2 * x$nvHvar + 1):(3 * x$nXvar + 3 * x$nuHvar +
            3 * x$nvHvar), , drop = FALSE],
          P.values = TRUE,
          digits = digits, signif.legend = FALSE
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          cat(
            centerText("Deterministic Component of SFA for latent class 4",
              width = lengthSum + 2 + switch(dimCoefTable,
                `4` = 18, `5` = 31, `6` = 43, `7` = 57
              )
            ),
            "\n"
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          printCoefmat(mleRes[(3 * x$nXvar + 3 * x$nuHvar +
            3 * x$nvHvar + 1):(4 * x$nXvar + 3 * x$nuHvar +
            3 * x$nvHvar), , drop = FALSE],
          P.values = TRUE,
          digits = digits, signif.legend = FALSE
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          cat(
            centerText("Parameter in variance of u (one sided error) for latent class 4",
              width = lengthSum + 2 + switch(dimCoefTable,
                `4` = 18, `5` = 31, `6` = 43, `7` = 57
              )
            ),
            "\n"
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          printCoefmat(mleRes[(4 * x$nXvar + 3 * x$nuHvar +
            3 * x$nvHvar + 1):(4 * x$nXvar + 4 * x$nuHvar +
            3 * x$nvHvar), , drop = FALSE],
          P.values = TRUE,
          digits = digits, signif.legend = FALSE
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          cat(
            centerText("Parameters in variance of v (symmetric error) for latent class 4",
              width = lengthSum + 2 + switch(dimCoefTable,
                `4` = 18, `5` = 31, `6` = 43, `7` = 57
              )
            ),
            "\n"
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          printCoefmat(mleRes[(4 * x$nXvar + 4 * x$nuHvar +
            3 * x$nvHvar + 1):(4 * x$nXvar + 4 * x$nuHvar +
            4 * x$nvHvar), , drop = FALSE],
          P.values = TRUE,
          digits = digits, signif.legend = FALSE
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          cat(
            centerText("Deterministic Component of SFA for latent class 5",
              width = lengthSum + 2 + switch(dimCoefTable,
                `4` = 18, `5` = 31, `6` = 43, `7` = 57
              )
            ),
            "\n"
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          printCoefmat(mleRes[(4 * x$nXvar + 4 * x$nuHvar +
            4 * x$nvHvar + 1):(5 * x$nXvar + 4 * x$nuHvar +
            4 * x$nvHvar), , drop = FALSE],
          P.values = TRUE,
          digits = digits, signif.legend = FALSE
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          cat(
            centerText("Parameter in variance of u (one sided error) for latent class 5",
              width = lengthSum + 2 + switch(dimCoefTable,
                `4` = 18, `5` = 31, `6` = 43, `7` = 57
              )
            ),
            "\n"
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          printCoefmat(mleRes[(5 * x$nXvar + 4 * x$nuHvar +
            4 * x$nvHvar + 1):(5 * x$nXvar + 5 * x$nuHvar +
            4 * x$nvHvar), , drop = FALSE],
          P.values = TRUE,
          digits = digits, signif.legend = FALSE
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          cat(
            centerText("Parameters in variance of v (symmetric error) for latent class 5",
              width = lengthSum + 2 + switch(dimCoefTable,
                `4` = 18, `5` = 31, `6` = 43, `7` = 57
              )
            ),
            "\n"
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          printCoefmat(mleRes[(5 * x$nXvar + 5 * x$nuHvar +
            4 * x$nvHvar + 1):(5 * x$nXvar + 5 * x$nuHvar +
            5 * x$nvHvar), , drop = FALSE],
          P.values = TRUE,
          digits = digits, signif.legend = FALSE
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          cat(
            centerText("Estimated prior probabilities for class membership",
              width = lengthSum + 2 + switch(dimCoefTable,
                `4` = 18, `5` = 31, `6` = 43, `7` = 57
              )
            ),
            "\n"
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          printCoefmat(mleRes[(5 * x$nXvar + 5 * x$nuHvar +
            5 * x$nvHvar + 1):(5 * x$nXvar + 5 * x$nuHvar +
            5 * x$nvHvar + 4 * x$nZvar), , drop = FALSE],
          P.values = TRUE, digits = digits, signif.legend = TRUE
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
        }
      }
    }
  }
  cat(x$mleDate, "\n")
  cat("Log likelihood status:", x$optStatus, "\n")
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  invisible(x)
}
