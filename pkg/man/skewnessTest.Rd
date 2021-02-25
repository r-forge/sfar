\name{skewnessTest}
\alias{skewnessTest}

\title{Skewness test for object of class \code{'sfacross'}}

\description{Skewness test for object of class \code{'sfacross'}.}

\usage{
skewnessTest(object, test = "agostino")
}

\arguments{
\item{object}{An object of class \code{'sfacross'}, returned by \code{\link{sfacross}}.}

\item{test}{A character string specifying the test to implement. If \code{"agostino"} 
(default), D'agostino skewness test is implemented (D'Agostino et Pearson, 1973). 
If \code{"coelli"}, Coelli skewness test is implemented (Coelli, 1995).}
}

\value{
\code{skewnessTest} returns the results of either the D'Agostino's or the 
Coelli's skewness test.
}

\note{
\code{\link{skewnessTest}} is currently only available for object of class \code{'sfacross'}.
}

\examples{
## Data on fossil fuel fired steam electric power generation plants in the
## United States
# Translog SFA (cost function) truncated normal with scaling property
tl_u_ts <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
    log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
    I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
    udist = "tnormal", muhet = ~ regu, uhet = ~ regu, data = utility, S = -1,
    scaling = TRUE, method = "mla")
  skewnessTest(tl_u_ts)
}

\references{
Coelli, T. 1995. Estimators and hypothesis tests for a stochastic frontier
function - a Monte-Carlo analysis. \emph{Journal of Productivity Analysis},
\bold{6}:247--268.

D'Agostino, R., and E.S. Pearson. 1973. Tests for departure from normality.
Empirical results for the distributions of \eqn{b_2} and \eqn{\sqrt{b_1}}.
\emph{Biometrika}, \bold{60}:613--622.
}

\keyword{methods}
