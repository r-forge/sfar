# skewness test for sfacross ----------

#' @title Skewness test for object of class \code{"sfacross"}
#'
#' @name skewnessTest
#'
#' @aliases skewnessTest.sfacross
#'
#' @param object An object of either class \code{"sfacross"}, returned by the
#'   function \code{\link{sfacross}}.
#' @param test character string. If \code{"agostino"} D'agostino skewness test is
#' returned. If \code{"coelli"} Coelli (1995) skewness test is returned
#' @param ... Currently ignored.
#'
#' @return \code{skewnessTest} returns either the D'agostino or the Coelli
#' skewness test
#'
#' @export
#' @export skewnessTest
#'
#' @note \code{skewnessTest} is only available for object of class \code{"sfacross"}.
#'
#' @references
#'
#' Coelli, T. 1995. Estimators and Hypothesis Tests for a Stochastic Frontier
#' Function - a Monte-Carlo Analysis. \emph{Journal of Productivity Analysis},
#' \bold{6}:247-268.
#'
#' D'Agostino, R., and E.S. Pearson. 1973. Tests for departure from normality.
#' Empirical results for the distributions of \eqn{b_2} and \eqn{\sqrt{b_1}}.
#' \emph{Biometrika}, \bold{60}:613-622.
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
#' skewnessTest(tl_u_ts)
#' @keywords methods skewness
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

