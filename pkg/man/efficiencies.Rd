\name{efficiencies}

\alias{efficiencies}
\alias{efficiencies.sfacross}
\alias{efficiencies.lcmcross}

\title{Returning conditional (in-)efficiency estimates}

\usage{
\method{efficiencies}{sfacross}(object, indataTable = FALSE, level = 0.95, ...)

\method{efficiencies}{lcmcross}(object, indataTable = FALSE, level = 0.95, ...)
}

\arguments{
\item{object}{A (latent class) stochastic frontier model returned by
\code{\link{sfacross}} or \code{\link{lcmcross}}.}

\item{indataTable}{Logical. If \code{TRUE}, (in)efficiency scores are
returned in addition to the information in the \code{dataTable} object.
When \code{FALSE}, only the (in)efficiency scores are returned.}

\item{level}{A number between between 0 and 0.9999 used for the computation
of (in)efficiency confidence intervals. Only when \code{udist} is
\code{"hnormal"}, \code{"exponential"}, \code{"tnormal"} or
\code{"uniform"}.}

\item{...}{Currently ignored.}
}

\value{
A dataframe that contains (in-)efficiency estimates.

\bold{- For object of class \code{"sfacross"} the following elements are
returned:}

\item{u}{Conditional inefficiency. In the case argument \code{udist} of
\link{sfacross} is set to \code{"uniform"} two conditional inefficiency
estimates are returned: \code{u1} for the classic conditional inefficiency
following Jondrow \emph{et al.} (1982) and \code{u2} which is obtained when
\eqn{\theta/\sigma_v \longrightarrow \infty} (see Nguyen,
2010).}

\item{uLB}{Lower bound for conditional inefficiency. Only when the argument
\code{udist} of \link{sfacross} is set to \code{"hnormal"},
\code{"exponential"}, \code{"tnormal"} or \code{"uniform"}.}

\item{uUB}{Upper bound for conditional inefficiency. Only when the argument
\code{udist} of \link{sfacross} is set to \code{"hnormal"},
\code{"exponential"}, \code{udist = "tnormal"} or \code{"uniform"}.}

\item{teJLMS}{\eqn{\exp{(-u)}}. When the argument \code{udist} of
\link{sfacross} is set to \code{"uniform"} \code{teJLMS1} = \eqn{\exp{(-u1)}}
and \code{teJLMS2} = \eqn{\exp{(-u2)}}. Only when \code{logDepVar = TRUE}.}

\item{m}{Conditional model. Only when the argument \code{udist} of
\link{sfacross} is set to \code{"hnormal"}, \code{"exponential"},
\code{"tnormal"}, or \code{"rayleigh"}.}

\item{teMO}{\eqn{\exp{(-m)}}. Only when, in the function \link{sfacross},
\code{logDepVar = TRUE} and \code{udist = "hnormal"}, \code{"exponential"},
\code{"tnormal"}, \code{"uniform"}, or \code{"rayleigh"}.}

\item{teBC}{Battese and Coelli (1988) conditional efficiency. Only when, in
the function \link{sfacross}, \code{logDepVar = TRUE} and \code{udist =
"hnormal"}, \code{"exponential"}, \code{"tnormal"}, \code{"genexponential"},
\code{"rayleigh"}, or \code{"tslaplace"}. In the case \code{udist =
"uniform"}, two conditional efficiency estimates are returned: \code{teBC1}
which is the classic conditional efficiency following and Battese and Coelli
(1988) and \code{teBC2} when \eqn{\theta/\sigma_v \longrightarrow \infty}
(see Nguyen, 2010).}

\item{teBCLB}{Lower bound for Battese and Coelli (1988) conditional
efficiency. Only when, in the function \link{sfacross}, \code{logDepVar =
TRUE} and when \code{udist = "hnormal"}, \code{"exponential"},
\code{"tnormal"}, or \code{"uniform"}.}

\item{teBCUB}{Upper bound for Battese and Coelli (1988) conditional
efficiency. Only when, in the function \link{sfacross}, \code{logDepVar =
TRUE} and when \code{udist = "hnormal"}, \code{"exponential"},
\code{"tnormal"}, or \code{"uniform"}.}

\bold{- For object of class \code{"lcmcross"} the following elements are
returned:}

\item{u_c}{Conditional inefficiency of the most probable class given the
posterior probability.}

\item{teJLMS_c}{\eqn{\exp{(-u_c)}}. Only when, in the function
\link{lcmcross}, \code{logDepVar = TRUE}.}

\item{PosteriorProb_c}{Highesht posterior probability.}

\item{Group_c}{Class of each observation.} Class assignment is based on the
largest posterior probability. This probability is obtained using Bayes' rule
and for class j whe have

\deqn{w\left(j|i\right)=\frac{P\left(i|j\right)\pi\left(i,
j\right)}{\sum_{m=1}^JP\left(i|m\right)\pi\left(i, m\right)}}

\item{u_c#}{conditional inefficiency associated to class #.}

\item{PosteriorProb_c#}{Posterior probability of class #.}

\item{PriorProb_c#}{Prior probability of class #.}
}
\description{
This method returns (in-)efficiency estimates from (latent
class) stochastic frontier models estimated with \code{\link{sfacross}} or
\code{\link{lcmcross}}
}
\details{
Conditional inefficiency are obtained following Jondrow \emph{et
al.} (1982) and the conditional efficiency is computed following Battese
and Coelli (1988). In some cases the conditional mode is also returned
(Jondrow \emph{et al.} 1982). In other cases the confidence interval is
computed following Horrace and Schmidt (1996) but also Hjalmarsson \emph{et
al.} (1996) and Berra and Sharma (1999).

In the case of the half normal distribution for the one-sided error term:
\itemize{
\item The conditional inefficiency is
}

\deqn{E\left[u_i|\epsilon_i\right]=\mu_{i*} +
\sigma_*\frac{\phi\left(\frac{\mu_{i*}}{\sigma_*}\right)}{\Phi\left(\frac{\mu_{i*}}{\sigma_*}\right)}}

where

\deqn{\mu_{i*}=\frac{-S\epsilon_i\sigma_u^2}{\sigma_u^2 + \sigma_v^2}}

and

\deqn{\sigma_*^2 = \frac{\sigma_u^2 \sigma_v^2}{\sigma_u^2 + \sigma_v^2}}
\itemize{
\item The Battese and Coelli (1988) conditional efficiency score is obtained by:
}

\deqn{E\left[e^{-u_i}|\epsilon_i\right] = e^{\left(-\mu_{i*}+\frac{1}{2}\sigma_*^2\right)}
 \frac{\Phi\left(\frac{\mu_{i*}}{\sigma_*}-\sigma_*\right)}{\Phi\left(\frac{\mu_{i*}}{\sigma_*}\right)}}
\itemize{
\item The conditional mode is computed using:
}

\deqn{M\left[u_i|\epsilon_i\right]= \mu_{i*} \quad For \quad \mu_{i*} > 0}

and

\deqn{M\left[u_i|\epsilon_i\right]= 0 \quad For \quad \mu_{i*} \leq 0}
\itemize{
\item The confidence intervals are obtained by:
}

\deqn{\mu_{i*} + Z_L\sigma_* \leq E\left[u_i|\epsilon_i\right] \leq \mu_{i*} + Z_U\sigma_*}

with \eqn{L_i = \mu_{i*} + Z_L\sigma_*} and \eqn{U_i = \mu_{i*} + Z_U\sigma_*}

and

\deqn{Z_L = \Phi^{-1}\left\{1 - \left(1-\frac{\alpha}{2}\right)\left[1-\Phi\left(-\frac{\mu_{i*}}{\sigma_*}\right)\right]\right\}}

and

\deqn{Z_U = \Phi^{-1}\left\{1-\frac{\alpha}{2}\left[1-\Phi\left(-\frac{\mu_{i*}}{\sigma_*}\right)\right]\right\}}

For the efficiency score we have

\deqn{\exp{\left(-U_i\right)} \leq E\left[\exp{\left(-u_i\right)}|\epsilon_i\right] \leq \exp{\left(-L_i\right)}}
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
  efficiencies(tl_u_ts)

## Using data on eighty-two countries production (DGP)
# LCM Cobb Douglas (production function) half normal distribution
cb_2c_h <- lcmcross(formula = ly ~ lk + ll + yr, udist = "hnormal",
                    data = worldprod, S = 1, method = "ucminf")
  efficiencies(cb_2c_h,  indataTable = TRUE)
}

\references{
Battese, G.E., and T.J. Coelli. 1988. Prediction of Firm-Level Technical
Efficiencies with a Generalized Frontier Production Function and Panel Data.
\emph{Journal of econometrics}, \bold{38}:387--399.

Bera, A.K., and S.C. Sharma. 1999. Estimating production uncertainty in
stochastic frontier production function models. \emph{Journal of Productivity
Analysis}, \bold{12}:187-210.

Hjalmarsson, L., S.C. Kumbhakar, and A. Heshmati. 1996. DEA, DFA and SFA: A
comparison. \emph{Journal of Productivity Analysis}, \bold{7}:303-327.

Horrace, W.C., and P. Schmidt. 1996. Confidence statements for efficiency
estimates from stochastic frontier models. \emph{Journal of Productivity
Analysis}, \bold{7}:257-282.

Jondrow, J., C.A.K. Lovell, I.S. Materov, and P. Schmidt. 1982. On the
estimation of technical inefficiency in the stochastic frontier production
function model. \emph{Journal of econometrics}, \bold{19}:233--238.

Nguyen, N.B. 2010. Estimation of technical efficiency in stochastic frontier
analysis. Bowling Green State University.
}

\seealso{
\code{\link{sfacross}}, for the stochastic frontier analysis model fitting
function.

\code{\link{lcmcross}}, for the latent class stochastic frontier analysis
model fitting function.
}

\author{K Hervé Dakpo, Yann Desjeux and Laure Latruffe}

\keyword{efficiencies}
\keyword{methods}
