\name{ic}

\alias{ic}
\alias{ic.sfacross}
\alias{ic.lcmcross}
\title{Extract information criteria}

\description{
\code{\link{ic}} returns information criterion from classic or latent class stochastic frontier 
models estimated with \code{\link{sfacross}} or \code{\link{lcmcross}.}
}

\usage{
\method{ic}{sfacross}(object, IC = "AIC", ...)

\method{ic}{lcmcross}(object, IC = "AIC", ...)
}

\arguments{
\item{object}{A classic or latent class stochastic frontier model returned by
\code{\link{sfacross}} or \code{\link{lcmcross}}.}

\item{IC}{Character string. Information criterion measure. Three criteria are available:
\itemize{
\item \code{"AIC"} for Akaike information criterion (default)
\item \code{"BIC"} for Bayesian information criterion
\item \code{"HQIC"} for Hannan-Quinn information criterion
}}

\item{...}{Currently ignored.}
}

\details{
The different information criteria are computed as follows:
\itemize{
\item AIC: \eqn{-2 \log{LL} + 2 * K}
\item BIC: \eqn{-2 \log{LL} + \log{N} * K}
\item HQIC: \eqn{-2 \log{LL} + 2 \log{\left[\log{N}\right]} * K}
}
where \eqn{LL} is the maximum likelihood value, \eqn{K} the number of parameters
estimated and \eqn{N} the number of observations.
}

\value{
\code{\link{ic}} returns the value of the information criterion (AIC, BIC or
HQIC) of the maximum likelihood coefficients.
}

\examples{
## Using data on eighty-two countries production (DGP)
# LCM Cobb Douglas (production function) half normal distribution
cb_2c_h <- lcmcross(formula = ly ~ lk + ll + yr, udist = "hnormal", data = worldprod)
  ic(cb_2c_h)
  ic(cb_2c_h, IC = "BIC")
  ic(cb_2c_h, IC = "HQIC")
}

\seealso{
\code{\link{sfacross}}, for the stochastic frontier analysis model fitting
function.

\code{\link{lcmcross}}, for the latent class stochastic frontier analysis
model fitting function.
}

\author{K Hervé Dakpo, Yann Desjeux and Laure Latruffe}

\keyword{methods}
