\name{summary}

\alias{summary}
\alias{summary.sfacross}
\alias{print.summary.sfacross}
\alias{summary.lcmcross}
\alias{print.summary.lcmcross}

\title{Summary methods for frontier objects}

\description{Create and print summary results of objects of class 
\code{'sfacross'} and \code{'lcmcross'} returned by \code{\link{sfacross}}
and \code{\link{lcmcross}}.
}

\usage{
\method{summary}{sfacross}(object, grad = FALSE, ci = FALSE, ...)

\method{print}{summary.sfacross}(x, digits = max(3, getOption("digits") - 2), ...)

\method{summary}{lcmcross}(object, grad = FALSE, ci = FALSE, ...)

\method{print}{summary.lcmcross}(x, digits = max(3, getOption("digits") - 2), ...)
}

\arguments{
\item{object}{An object of either class \code{'sfacross'}, returned by the
function \code{\link{sfacross}}, or class \code{'lcmcross'}, returned by
the function \code{\link{lcmcross}}.}

\item{grad}{Logical. Default = \code{FALSE}. If \code{TRUE}, the gradient 
for the maximum likelihood (ML) estimates of the different parameters is returned.}

\item{ci}{Logical. Default = \code{FALSE}. If \code{TRUE}, the 95\% confidence 
interval for the different parameters (OLS and ML estimates) is returned.}

\item{...}{Currently ignored.}

\item{x}{An object of either class \code{'summary.sfacross'} or
\code{'summary.lcmcross'}.}

\item{digits}{Numeric. Number of digits displayed in values.}
}

\value{
The \code{\link{summary}} method returns a list of class \code{'summary.sfacross'}
or \code{'summary.lcmcross'} that contains the same elements as an object returned by
\code{\link{sfacross}} or \code{\link{lcmcross}} but the following additional elements:

\item{AIC}{Akaike information criterion.}

\item{BIC}{Bayesian information criterion.}

\item{HQIC}{Hannan-Quinn information criterion.}

\item{sigmavSq}{For \code{object} of class \code{'sfacross'}. Variance of the two-sided error
term (\eqn{\sigma_v^2}).}

\item{sigmauSq}{For \code{object} of class \code{'sfacross'}. Variance of the one-sided error
term (\eqn{\sigma_u^2}).}

\item{theta}{For \code{object} of class \code{'sfacross'} with \code{"udist = uniform"}.
\code{theta} value in the case the unifom distribution is defined as: \eqn{u_i \in [0, \theta]}.}

\item{Varu}{For \code{object} of class \code{'sfacross'}. Variance of the one-sided error term.}

\item{Eu}{For \code{object} of class \code{'sfacross'}. Expected unconditional inefficiency.}

\item{Expu}{For \code{object} of class \code{'sfacross'}. Expected unconditional efficiency.}

\item{olsRes}{For \code{object} of class \code{'sfacross'}. Matrix of OLS estimates, their
standard errors, t-values, P-values, and when \code{ci = TRUE} their
confidence intervals.}

\item{mleRes}{Matrix of ML estimates, their standard errors, z-values,
asymptotic P-values, and when \code{grad = TRUE} their
gradient, \code{ci = TRUE} their confidence intervals.}

\item{chisq}{For \code{object} of class \code{'sfacross'}. Chi-square statistics of the
difference between the stochastic frontier and the OLS.}

\item{df}{Degree of freedom for the inefficiency model.}

\item{digits}{Number of digits displayed.}
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
  summary(cb_2c_h, grad = TRUE, ci = TRUE)
}

\seealso{
\code{\link{sfacross}}, for the stochastic frontier analysis model fitting
function.

\code{\link{lcmcross}}, for the latent class stochastic frontier analysis
model fitting function.
}

\keyword{methods}
