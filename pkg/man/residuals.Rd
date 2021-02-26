\name{residuals}

\alias{residuals}
\alias{residuals.sfacross}
\alias{residuals.lcmcross}

\title{Residuals of frontier estimation}

\description{
It returns the residuals' values from classic or latent class 
stochastic frontier models estimated with \code{\link{sfacross}} or
\code{\link{lcmcross}}.
}

\usage{
\method{residuals}{sfacross}(object, ...)

\method{residuals}{lcmcross}(object, ...)
}

\arguments{
\item{object}{A classic or a latent class stochastic frontier model 
returned by \code{\link{sfacross}} or \code{\link{lcmcross}}.}

\item{...}{Currently ignored.}
}

\value{
When the \code{object} is of class \code{'sfacross'}, \code{link{residuals}} 
returns a vector of residuals values.

When the \code{object} is of class \code{'lcmcross'}, \code{link{residuals}} 
returns a data frame containing the residuals values for each latent class.
}

\note{The residuals values are ordered in the same way as the
corresponding observations in the dataset used for the estimation.
}

\examples{
## Data on fossil fuel fired steam electric power generation plants in the
## United States
# Translog SFA (cost function) truncated normal with scaling property
tl_u_ts <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
    log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
    I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
    udist = 'tnormal', muhet = ~ regu, uhet = ~ regu, data = utility, S = -1,
    scaling = TRUE, method = 'mla')
  residuals(tl_u_ts)

## Using data on eighty-two countries production (DGP)
# LCM Cobb Douglas (production function) half normal distribution
cb_2c_h <- lcmcross(formula = ly ~ lk + ll + yr, udist = 'hnormal', data = worldprod, 
    S = 1, method = 'ucminf')
  resid <- residuals(cb_2c_h)
  head(resid)
}

\seealso{
\code{\link{sfacross}}, for the stochastic frontier analysis model fitting
function.

\code{\link{lcmcross}}, for the latent class stochastic frontier analysis
model fitting function.
}

\author{K Hervé Dakpo, Yann Desjeux and Laure Latruffe}

\keyword{methods}
