# skewness test for sfacross ----------

skewnessTest <- function(object, test = "agostino") {
if (! inherits(object,"sfacross")) {
        stop("Argument 'object' must be a of class 'sfacross'")
}
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

