\name{vcov}

\alias{vcov}
\alias{vcov.sfacross}
\alias{vcov.lcmcross}

\title{Variance-covariance matrix for frontier objects}

\description{
Extract the variance-covariance matrix of the maximum likelihood (ML) 
coefficients of classic or latent class stochastic frontier models estimated by 
\code{\link{sfacross}} or \code{\link{lcmcross}}.
}

\usage{
\method{vcov}{sfacross}(object, extraPar = FALSE, ...)

\method{vcov}{lcmcross}(object, ...)
}

\arguments{
\item{object}{A classic or latent class stochastic frontier model 
returned by \code{\link{sfacross}} or \code{\link{lcmcross}}}

\item{extraPar}{Logical. Only available for objects of class \code{'sfacross'}. 
Default = \code{FALSE}. If \code{TRUE}, variances and covariances of additional 
parameters are returned:

\code{sigmaSq} = \code{sigmauSq} + \code{sigmavSq}

\code{lambdaSq} = \code{sigmauSq}/\code{sigmavSq}

\code{sigmauSq} = \eqn{\exp{(Wu)}}

\code{sigmavSq} = \eqn{\exp{(Wv)}}

\code{sigma} = \code{sigmaSq}^0.5

\code{lambda} = \code{lambdaSq}^0.5

\code{sigmau} = \code{sigmauSq}^0.5

\code{sigmav} = \code{sigmavSq}^0.5

\code{gamma} = \code{sigmauSq}/(\code{sigmauSq} + \code{sigmavSq})}

\item{...}{Currently ignored}
}

\value{
\code{\link{vcov}} returns the covariance matrix of the maximum likelihood
coefficients.
}

\details{
The variance-covariance matrix is obtained by the inversion of the
negative hessian matrix. Depending on the distribution and the
\code{"hessianType"} option, the analytical/numeric hessian or the bhhh
hessian or the robust hessian matrix is evaluated.

The argument \code{extraPar}, is currently available for objects of class \code{'sfacross'}. 
When \code{"extraPar = TRUE"}, the variance-covariance of the additional 
parameters is obtained using the delta method.
}

\examples{
## Data on Spanish dairy farms
# Cobb Douglas (production function) half normal distribution
cb_s_h <- sfacross(formula = YIT ~ X1 + X2 + X3 + X4, udist = "hnormal", 
    data = dairyspain, S = 1, method = "bfgs")
  vcov(cb_s_h)
  vcov(cb_s_h, extraPar = TRUE)


## Data on eighty-two countries production (DGP)
# LCM Cobb Douglas (production function) half normal distribution
cb_2c_h <- lcmcross(formula = ly ~ lk + ll + yr, udist = "hnormal",
    data = worldprod, uhet = ~ initStat, S = 1, method = "ucminf")
  vcov(cb_2c_h)
}

\seealso{
\code{\link{sfacross}}, for the stochastic frontier analysis model fitting function.

\code{\link{lcmcross}}, for the latent class stochastic frontier analysis model fitting function.
}

\author{K Hervé Dakpo, Yann Desjeux and Laure Latruffe}

\keyword{methods}
