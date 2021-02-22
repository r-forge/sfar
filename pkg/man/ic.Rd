\name{ic}

\alias{ic}
\alias{ic.sfacross}
\alias{ic.lcmcross}
\title{Information criteria}

\usage{
\method{ic}{sfacross}(object, IC = "AIC", ...)

\method{ic}{lcmcross}(object, IC = "AIC", ...)
}

\arguments{
\item{object}{A (latent class) stochastic frontier model returned by
\code{\link{sfacross}} or \code{\link{lcmcross}}.}

\item{IC}{charater string. Information criterion measure. Three criterion are
available:
\itemize{
\item \code{"AIC"} for Akaike information criterion
\item \code{"BIC"} for Bayesian information criterion
\item \code{"HQIC"} for Hannan-Quinn information criterion
}}

\item{...}{Currently ignored.}
}

\value{
\code{ic} returns the value of the information criterion (AIC, BIC or
HQIC) of the maximum likelihood coefficients.
}
\description{
This method returns information criterion from (latent class)
stochastic frontier models estimated with \code{\link{sfacross}} or
\code{\link{lcmcross}}
}
\details{
The different information criterion are computed as follows:
\itemize{
\item AIC: \eqn{-2 \log{LL} + 2 * K}
\item BIC: \eqn{-2 \log{LL} + \log{N} * K}
\item HQIC: \eqn{-2 \log{LL} + 2 \log{\left[\log{N}\right]} * K}
}

where LL is the maximum likelihood value, K the number of parameters
estimated and N the number of observations.
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
  ic(tl_u_ts, IC = "AIC")

## Using data on eighty-two countries production (DGP)
# LCM Cobb Douglas (production function) half normal distribution
cb_2c_h <- lcmcross(formula = ly ~ lk + ll + yr, udist = "hnormal",
                    data = worldprod, S = 1, method = "ucminf")
  ic(cb_2c_h, IC = "AIC")
}

\seealso{
\code{\link{sfacross}}, for the stochastic frontier analysis model fitting
function.

\code{\link{lcmcross}}, for the latent class stochastic frontier analysis
model fitting function.
}

\author{K Hervé Dakpo, Yann Desjeux and Laure Latruffe}

\keyword{AIC}
\keyword{BIC}
\keyword{HQIC}
\keyword{methods}
