# residuals from sfacross ----------

residuals.sfacross <- function(object, ...) {
  object$dataTable$mleResiduals
}

# residuals from lcmcross ----------

residuals.lcmcross <- function(object, ...) {
  if (object$nClasses == 2) {
    data.frame(select(object$dataTable, "mleResiduals_c1",
      "mleResiduals_c2"))
  } else {
    if (object$nClasses == 3) {
      data.frame(select(object$dataTable, "mleResiduals_c1",
        "mleResiduals_c2", "mleResiduals_c3"))
    } else {
      if (object$nClasses == 4) {
        data.frame(select(object$dataTable, "mleResiduals_c1",
          "mleResiduals_c2", "mleResiduals_c3",
          "mleResiduals_c4"))
      } else {
        if (object$nClasses == 5) {
          data.frame(select(object$dataTable, "mleResiduals_c1",
          "mleResiduals_c2", "mleResiduals_c3",
          "mleResiduals_c4", "mleResiduals_c5"))
        }
      }
    }
  }
}
