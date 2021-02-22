# skewness test for sfacross ----------

skewnessTest.sfacross <- function(object, test = "agostino", ...) {
  if (test == "agostino") {
    object$AgostinoTest
  } else {
    if (test == "coelli") {
      object$CoelliM3Test
    } else {
      stop("argument 'test' must be either 'agostino', or 'coelli'",
           call. = FALSE)
    }
  }
}

