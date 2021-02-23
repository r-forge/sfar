# fitted values for sfacross ----------

fitted.sfacross <- function(object, ...) {
  object$dataTable$mleFitted
}

# fitted values for lcmcross ----------

fitted.lcmcross <- function(object, ...) {
  if (object$nClasses == 2) {
    data.frame(select(object$dataTable, "mleFitted_c1",
      "mleFitted_c2"))
  } else {
    if (object$nClasses == 3) {
      data.frame(select(object$dataTable, "mleFitted_c1",
        "mleFitted_c2", "mleFitted_c3"))
    } else {
      if (object$nClasses == 4) {
        data.frame(select(object$dataTable, "mleFitted_c1",
          "mleFitted_c2", "mleFitted_c3",
          "mleFitted_c4"))
      } else {
        if (object$nClasses == 5) {
          data.frame(select(object$dataTable, "mleFitted_c1",
          "mleFitted_c2", "mleFitted_c3",
          "mleFitted_c4", "mleFitted_c5"))
        }
      }
    }
  }
}
