# marginal effects computation sfacross ----------

#' @title Marginal effects of the determinants of inefficiency
#'
#' @name marginal
#'
#' @aliases marginal.sfacross
#'
#' @description This function returns marginal effects of the determinants of
#'   efficiency from (Latent class) stochastic frontier models estimated with
#'   \code{\link{sfacross}} or \code{\link{lcmcross}}
#'
#' @param object A (latent class) stochastic frontier model returned by
#'   \code{\link{sfacross}} or \code{\link{lcmcross}}.
#' @param indataTable Logical. If \code{TRUE}, marginal effects are returned in
#'   addition to the information in the \code{dataTable} object. When \code{FALSE}
#'   (Default) only the marginal effects are returned.
#'   \code{FALSE}.
#' @param ... Currently ignored.
#'
#' @details The function \code{marginal} operates in the presence of exogenous
#'   variables that explain inefficiency (\eqn{uhet = ~ Z_u} or \eqn{muhet = ~
#'   Z_{mu}}).
#'
#'   Two components are computed for each variables: the marginal effects on the
#'   expected inefficiency (\eqn{\frac{\partial E[u]}{\partial Z_{(m)u}}}) and the
#'   marginal effects on the variance of inefficiency (\eqn{\frac{\partial
#'   V[u]}{\partial Z_{(m)u}}}).
#'
#' The model also accommodates the Wang (2002) parametrization of
#' \eqn{\mu} and \eqn{\sigma_u^2} by the same vector of exogenous variables.
#' This double parameterization accounts for non-monotonic relationships
#' between the inefficiency and its determinants.
#'
#' @return A dataframe containing the marginal effects of the \code{Z} variables
#'   on the expected inefficiency (each variable has the prefix \code{"Eu_"})
#'   and on the variance of the inefficiency (each variable has the prefix
#'   \code{"Vu_"}) is returned.
#'
#'   In the case of the latent class model, each variable terminates with the
#'   suffix \code{"_c#"} where \code{"#"} is the class number.
#'
#' @author K Hervé Dakpo, Yann Desjeux and Laure Latruffe
#'
#' @references
#'
#' Wang, H.J. 2002. Heteroscedasticity and non-monotonic efficiency effects of a
#' stochastic frontier model. emph{Journal of Productivity Analysis},
#' \bold{18}:241--253.
#'
#' @export
#' @export marginal
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
#' udist = "tnormal", muhet = ~ regu, uhet = ~ regu, data = utility, S = -1,
#' scaling = TRUE, method = "mla")
#'
#' marginal(tl_u_ts)
#'
#' ## Using data on eighty-two countries production (DGP)
#'
#' # LCM Cobb Douglas (production function) half normal distribution
#'
#' cb_2c_h <- lcmcross(formula = ly ~ lk + ll + yr, udist = "hnormal",
#' data = worldprod, uhet = ~ initStat, S = 1, method = "ucminf")
#'
#' marginal(cb_2c_h, indataTable = TRUE)
#' @keywords methods marginal
marginal.sfacross <- function(object, indataTable = FALSE, ...) {
  if (length(indataTable) != 1 || !is.logical(indataTable[1])) {
    stop("argument 'indataTable' must be a single logical value",
      call. = FALSE
    )
  }
  if (object$udist == "hnormal") {
    if (object$nuHvar == 1) {
      warning("marginal effects cannot be computed for homoscedastic models",
        call. = FALSE
      )
    } else {
      EffMarg <- bind_cols(
        as_tibble(cmarghalfnorm_Eu(object = object)),
        as_tibble(cmarghalfnorm_Vu(object = object))
      )
    }
  } else {
    if (object$udist == "exponential") {
      if (object$nuHvar == 1) {
        warning("marginal effects cannot be computed for homoscedastic models",
          call. = FALSE
        )
      } else {
        EffMarg <- bind_cols(
          as_tibble(cmargexponorm_Eu(object = object)),
          as_tibble(cmargexponorm_Vu(object = object))
        )
      }
    } else {
      if (object$udist == "gamma") {
        if (object$nuHvar == 1) {
          warning("marginal effects cannot be computed for homoscedastic models",
            call. = FALSE
          )
        } else {
          EffMarg <- bind_cols(
            as_tibble(cmarggammanorm_Eu(object = object)),
            as_tibble(cmarggammanorm_Vu(object = object))
          )
        }
      } else {
        if (object$udist == "rayleigh") {
          if (object$nuHvar == 1) {
            warning("marginal effects cannot be computed for homoscedastic models",
              call. = FALSE
            )
          } else {
            EffMarg <- bind_cols(
              as_tibble(cmargraynorm_Eu(object = object)),
              as_tibble(cmargraynorm_Vu(object = object))
            )
          }
        } else {
          if (object$udist == "uniform") {
            if (object$nuHvar == 1) {
              warning("marginal effects cannot be computed for homoscedastic models",
                call. = FALSE
              )
            } else {
              EffMarg <- bind_cols(
                as_tibble(cmarguninorm_Eu(object = object)),
                as_tibble(cmarguninorm_Vu(object = object))
              )
            }
          } else {
            if (object$udist == "tnormal") {
              if (object$scaling) {
                if (object$nuHvar == 1) {
                  warning("marginal effects cannot be computed for homogeneous or homoscedastic models",
                    call. = FALSE
                  )
                } else {
                  EffMarg <- bind_cols(
                    as_tibble(cmargtruncnormscal_Eu(object = object)),
                    as_tibble(cmargtruncnormscal_Vu(object = object))
                  )
                }
              } else {
                if (object$nmuHvar == 1 & object$nuHvar ==
                  1) {
                  warning("marginal effects cannot be computed for homogeneous or homoscedastic models",
                    call. = FALSE
                  )
                } else {
                  EffMarg <- bind_cols(
                    as_tibble(cmargtruncnorm_Eu(object = object)),
                    as_tibble(cmargtruncnorm_Vu(object = object))
                  )
                }
              }
            } else {
              if (object$udist == "lognormal") {
                if (object$nmuHvar == 1 & object$nuHvar ==
                  1) {
                  warning("marginal effects cannot be computed for homogeneous or homoscedastic models",
                    call. = FALSE
                  )
                } else {
                  EffMarg <- bind_cols(
                    as_tibble(cmarglognorm_Eu(object = object)),
                    as_tibble(cmarglognorm_Vu(object = object))
                  )
                }
              }
            }
          }
        }
      }
    }
  }
  if (indataTable) {
    return(bind_cols(object$dataTable, EffMarg))
  } else {
    return(EffMarg)
  }
}

# marginal effects computation lcmcross ----------
#' @rdname marginal
#' @aliases marginal.lcmcross
#' @export
marginal.lcmcross <- function(object, indataTable = FALSE, ...) {
  if (length(indataTable) != 1 || !is.logical(indataTable[1])) {
    stop("argument 'indataTable' must be a single logical value",
      call. = FALSE
    )
  }
  if (object$nuHvar == 1) {
    warning("marginal effects cannot be computed for homoscedastic models",
      call. = FALSE
    )
  } else {
    if (object$nClasses == 2) {
      EffMarg <- bind_cols(
        as_tibble(cmargLCM2Chalfnorm_Eu(object = object)),
        as_tibble(cmargLCM2Chalfnorm_Vu(object = object))
      )
    } else {
      if (object$nClasses == 3) {
        EffMarg <- bind_cols(
          as_tibble(cmargLCM3Chalfnorm_Eu(object = object)),
          as_tibble(cmargLCM3Chalfnorm_Vu(object = object))
        )
      } else {
        if (object$nClasses == 4) {
          EffMarg <- bind_cols(
            as_tibble(cmargLCM4Chalfnorm_Eu(object = object)),
            as_tibble(cmargLCM4Chalfnorm_Vu(object = object))
          )
        } else {
          if (object$nClasses == 5) {
            EffMarg <- bind_cols(
              as_tibble(cmargLCM5Chalfnorm_Eu(object = object)),
              as_tibble(cmargLCM5Chalfnorm_Vu(object = object))
            )
          }
        }
      }
    }
  }
  if (indataTable) {
    return(bind_cols(object$dataTable, EffMarg))
  } else {
    return(EffMarg)
  }
}
