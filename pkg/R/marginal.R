# marginal effects computation sfacross ----------

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
