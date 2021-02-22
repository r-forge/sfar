# log likelihood extraction for sfacross ----------

logLik.sfacross <- function(object, individual = FALSE, ...) {
  if (length(individual) != 1 || !is.logical(individual[1]))
    stop("argument 'individual' must be a single logical value",
         call. = FALSE)
  if (individual) {
    LL <- object$dataTable$logL_OBS
  } else {
    LL <- object$mleLoglik
  }
  attr(LL, "Nobs") <- object$Nobs
  attr(LL, "df") <- object$nParm
  return(LL)
}

# log likelihood extraction for lcmcross ----------

logLik.lcmcross <- function(object, individual = FALSE, ...) {
  if (length(individual) != 1 || !is.logical(individual[1]))
    stop("argument 'individual' must be a single logical value",
         call. = FALSE)
  if (individual) {
    LL <- object$dataTable$logL_OBS
  } else {
    LL <- object$mleLoglik
  }
  attr(LL, "Nobs") <- object$Nobs
  attr(LL, "df") <- object$nParm
  return(LL)
}
