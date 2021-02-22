\name{logLik}

\alias{logLik}
\alias{logLik.sfacross}
\alias{logLik.lcmcross}

\title{Extract Log-Likelihood Value}

\usage{
\method{logLik}{sfacross}(object, individual = FALSE, ...)

\method{logLik}{lcmcross}(object, individual = FALSE, ...)
}

\arguments{
\item{object}{A (latent class) stochastic frontier model returned by
\code{\link{sfacross}} or \code{\link{lcmcross}}.}

\item{individual}{Logical. If \code{TRUE} a vector of individual likelihood is
returned. If \code{FALSE} (Default) the sum of all observations likelihood is
returned.}

\item{...}{Currently ignored.}
}

\value{
\code{logLik} returns an object of class logLik, which is either a
numeric scalar (the log-likelihood value when \code{"individual = FALSE"})
or a numeric vector (each observation likelihood when \code{"individual =
  TRUE"}) with 2 attributes: Nobs (total number of observations in all
equations) and df (number of free parameters, i.e. length of the coefficient
vector)
}

\description{
This method extracts the log-likelihood value(s) from stochastic
frontier models estimated with \code{\link{sfacross}} or \code{\link{lcmcross}}.
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
  logLik(tl_u_ts)

## Using data on eighty-two countries production (DGP)
# LCM Cobb Douglas (production function) half normal distribution
cb_2c_h <- lcmcross(formula = ly ~ lk + ll + yr, udist = "hnormal",
                    data = worldprod, S = 1, method = "ucminf")
  logLik(cb_2c_h)
}

\keyword{likelihood}
\keyword{methods}
