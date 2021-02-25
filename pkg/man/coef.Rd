\name{coef}

\alias{coef}
\alias{coef.summary.sfacross}
\alias{coef.summary.lcmcross}
\alias{coef.sfacross}
\alias{coef.lcmcross}

\title{coef methods for class \code{summary} of frontier objects}

\description{
Extract the coefficients, their standard errors, z-values, and
(asymptotic) P-values from models returned by the \code{summary} methods
for objects of class \code{sfacross} or \code{lcmcross}.

\code{coef} extracts coefficients for the likelihood estimation
}
\usage{
\method{coef}{summary.sfacross}(object, ...)

\method{coef}{summary.lcmcross}(object, ...)

\method{coef}{sfacross}(object, extraPar = FALSE, ...)

\method{coef}{lcmcross}(object, extraPar = FALSE, ...)
}

\arguments{
\item{object}{A (latent class) stochastic frontier model returned by
\code{\link{sfacross}} or \code{\link{lcmcross}}.}

\item{...}{Currently ignored.}

\item{extraPar}{Logical (default = \code{FALSE}). If \code{TRUE} additional parameters are also
returned:

\code{sigmaSq} = \code{sigmauSq} + \code{sigmavSq}

\code{lambdaSq} = \code{sigmauSq}/\code{sigmavSq}

\code{sigmauSq} = \eqn{\exp{(Wu)}}

\code{sigmavSq} = \eqn{\exp{(Wv)}}

\code{sigma} = \code{sigmaSq}^0.5

\code{lambda} = \code{lambdaSq}^0.5

\code{sigmau} = \code{sigmauSq}^0.5

\code{sigmav} = \code{sigmavSq}^0.5

\code{gamma} = \code{sigmauSq}/(\code{sigmauSq} + \code{sigmavSq})

In the case of object of class \code{lcmcross}, each additional
parameter terminates with \code{"#"} that represents the class
number.}
}
\value{
The \code{coef} method for objects of class summary.sfacross or
summary.lcmcross returns a matrix with four. These columns contain the estimated
coefficients, their standard errors, z-values, and (asymptotic) P-values.
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
  coef(tl_u_ts, extraPar = TRUE)
  coef(summary(tl_u_ts))
}

\seealso{
\code{\link{sfacross}}, for the stochastic frontier analysis model fitting
function.

\code{\link{lcmcross}}, for the latent class stochastic frontier analysis
model fitting function.
}

\author{K Hervé Dakpo, Yann Desjeux and Laure Latruffe}

\keyword{methods}
