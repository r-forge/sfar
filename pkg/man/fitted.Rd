\name{fitted}

\alias{fitted}
\alias{fitted.sfacross}
\alias{fitted.lcmcross}
\title{Fitted frontier value}

\description{
This method returns the fitted frontier values from classic or latent class 
stochastic frontier models estimated respectively with \code{\link{sfacross}} or
\code{\link{lcmcross}}.
}

\usage{
\method{fitted}{sfacross}(object, ...)

\method{fitted}{lcmcross}(object, ...)
}

\arguments{
\item{object}{A classic or a latent class stochastic frontier model returned by
\code{\link{sfacross}} or \code{\link{lcmcross}}, respectively.}

\item{...}{Currently ignored.}
}

\value{
\code{\link{fitted}} returns a vector of fitted values in the case of object
of class \code{sfacross}.

\code{\link{fitted}} returns a data frame containing the fitted values for each 
 class in the case of object of class \code{lcmcross}.

The fitted values are ordered in the same way as the corresponding observations
in the dataset used for the estimation.
}

\examples{
## Using data on eighty-two countries production (DGP)
# LCM Cobb Douglas (production function) half normal distribution
cb_2c_h <- lcmcross(formula = ly ~ lk + ll + yr, udist = 'hnormal', data = worldprod)
  fit.val <- fitted(cb_2c_h)
  head(fit.val)
}

\seealso{
\code{\link{sfacross}}, for the stochastic frontier analysis model fitting
function.

\code{\link{lcmcross}}, for the latent class stochastic frontier analysis
model fitting function.
}

\author{K Hervé Dakpo, Yann Desjeux and Laure Latruffe}

\keyword{methods}
