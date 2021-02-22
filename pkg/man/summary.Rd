\name{summary}

\alias{summary}
\alias{summary.sfacross}
\alias{print.summary.sfacross}
\alias{summary.lcmcross}
\alias{print.summary.lcmcross}

\title{Summary methods for frontier objects}

\usage{
\method{summary}{sfacross}(object, grad = FALSE, ci = FALSE, ...)

\method{print}{summary.sfacross}(x, digits = max(3, getOption("digits") - 2), ...)

\method{summary}{lcmcross}(object, grad = FALSE, ci = FALSE, ...)

\method{print}{summary.lcmcross}(x, digits = max(3, getOption("digits") - 2), ...)
}

\arguments{
\item{object}{An object of either class \code{"sfacross"}, returned by the
function \code{\link{sfacross}}, or class \code{"lcmcross"}, returned by
the function \code{\link{lcmcross}}.}

\item{grad}{Logical. If \code{TRUE} the gradient for the maximum likelihood
(ML) estimates of the different parameters is returned. Default =
\code{FALSE}.}

\item{ci}{Logical. If \code{TRUE} the 95\% confidence interval for the
different parameters (OLS and ML estimates) is returned. Default =
\code{FALSE}.}

\item{...}{Currently ignored.}

\item{x}{An object of either class \code{"summary.sfacross"} or
\code{"summary.lcmcross"} returned by the function \code{\link{summary}}.}

\item{digits}{Numeric. Number of digits displayed in values.}
}

\value{
The \code{summary} method returns a list of class \code{"summary.sfacross"}
or \code{"summary.lcmcross"} that is identical to an object returned by
\code{\link{sfacross}} or \code{\link{lcmcross}}, with the following
additional elements:

\item{AIC}{Akaike information criterion.}

\item{BIC}{Bayesian information criterion.}

\item{HQIC}{Hannan-Quinn information criterion.}

\item{sigmavSq}{Only for \code{"sfacross"}. Variance of the two-sided error
term (\eqn{\sigma_v^2}).}

\item{sigmauSq}{Only for \code{"sfacross"}. \eqn{\sigma_u^2}.}

\item{theta}{Only for \code{"sfacross"} and \code{"udist = uniform"}.
\code{theta} value in the case of the unifom distribution is defined as:
\eqn{u_i \in [0, \theta]}.}

\item{Varu}{Only for \code{"sfacross"}. Variance of the one-sided error term.}

\item{Eu}{Only for \code{"sfacross"}. Expected unconditional inefficiency.}

\item{Expu}{Only for \code{"sfacross"}. Expected unconditional efficiency.}

\item{olsRes}{Only for \code{"sfacross"}. Matrix of OLS estimates, their
standard errors, t-values, P-values, and when \code{ci = TRUE} their
confidence intervals.}

\item{mleRes}{matrix of ML estimates, their standard errors, z-values,
asymptotic P-values and when \code{grad = TRUE} and or \code{ci = TRUE} their
gradient and or their confidence intervals.}

\item{chisq}{Only for \code{"sfacross"}. Chi-square statistics of the
difference between the stochastic frontier and the OLS.}

\item{df}{Degree of freedom for the inefficiency model.}

\item{digits}{Number of digits displayed.}
}

\description{Create and print summary results of objects of class
\code{'sfacross'} and \code{'lcmcross'} (returned by \code{\link{sfacross}}
and \code{\link{lcmcross}}).
}

\details{
The summary methods creates objects of either class
\code{"summary.sfacross"} or class \code{"summary.lcmcross"} that extend
the objects of class \code{sfacross} or \code{lcmcross} with various
information about the estimated model like information criterion, see
\bold{Value}.

\code{coef} and \code{print} methods are associated (i.e.
\code{coef.summary.sfacross}, \code{print.summary.sfacross} and
\code{coef.summary.lcmcross}, \code{print.summary.lcmcross}).
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
  summary(tl_u_ts)

## Using data on eighty-two countries production (DGP)
# LCM Cobb Douglas (production function) half normal distribution
cb_2c_h <- lcmcross(formula = ly ~ lk + ll + yr, udist = "hnormal", data =
worldprod, S = 1, method = "ucminf")
  summary(cb_2c_h)
}

\seealso{
\code{\link{sfacross}}, for the stochastic frontier analysis model fitting
function.

\code{\link{lcmcross}}, for the latent class stochastic frontier analysis
model fitting function.
}

\keyword{methods}
\keyword{summary}
