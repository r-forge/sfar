\name{sfacross}

\alias{sfacross}
\alias{print.sfacross}

\title{Stochastic frontier estimation using cross/pooled sectional data}

\description{\code{sfacross} is a symbolic formula-based function for the estimation of
stochastic frontier models in the case of cross-sectional or pooled cross
section data, using maximum (simulated) likelihood - M(S)L.

Ten distributions are considered for the one-sided error term and nine
optimization algorithms are also available.

The function also accounts for heteroscedasticity in both one-sided and
two-sided error terms as in Reifschneider and Stevenson (1991), Caudill and
Ford (1993), Caudill \emph{et al.} (1995) and Hadri (1999), but also
heterogeneity in the mean of the pre-truncated distribution as in Kumbhakar
\emph{et al.} (1991), Huang and Liu (1994), Battese and Coelli (1995).

Finally, the truncated normal - normal distribution with scaling property as
in Wang and Schmidt (2002) is also implemented.
}

\usage{
sfacross(formula, muhet, uhet, vhet, data, subset, S = 1L, udist = "hnormal",
  scaling = FALSE, start = NULL, logDepVar = TRUE, method = "bfgs", hessianType = 1L,
  simType = "halton", Nsim = 100, prime = 2L, burn = 10, antithetics = FALSE,
  seed = 12345, itermax = 2000, printInfo = FALSE, tol = 1e-12, gradtol = 1e-06,
  stepmax = 0.1, qac = "marquardt")

\method{print}{sfacross}(x, digits = max(3, getOption("digits") - 2), ...)
}

\arguments{
\item{formula}{A symbolic description of the model to be estimated based on
the generic function \code{formula} (see section \sQuote{Details}).}

\item{muhet}{A one-part formula to consider heterogeneity in the mean of the
pre-truncated distribution (see section \sQuote{Details}).}

\item{uhet}{A one-part formula to consider heteroscedasticity in the
one-sided error variance (see section \sQuote{Details}).}

\item{vhet}{A one-part formula to consider heteroscedasticity in the
two-sided error variance (see section \sQuote{Details}).}

\item{data}{An optional data frame containing the data or the the variables
in the model.}

\item{subset}{An optional vector specifying a subset of observations to be
used in the optimization process.}

\item{S}{Integer. If \code{1} (Default) a production (profit) frontier is
estimated \eqn{\epsilon_i = v_i-u_i}. If \code{-1} a cost frontier is
estimated \eqn{\epsilon_i = v_i+u_i}.}

\item{udist}{Character string. Distribution specification for the one-sided error
term. Default = \code{"hnormal"}. 10 different distributions are available:
\itemize{
\item \code{'hnormal'}, for the half normal distribution (Aigner \emph{et al.}
1977, Meeusen and Vandenbroeck 1977)
\item \code{'exponential'}, for the exponential distribution
\item \code{tnormal} for the truncated normal distribution (Stevenson 1980)
\item \code{'rayleigh'}, for the rayleigh distribution (Hajargasht 2015)
\item \code{'uniform'}, for the uniform distribution (Li 1996, Nguyen 2010)
\item \code{'gamma'}, for the gamma distribution (Greene 2003)
\item \code{'lognormal'}, for the log normal distribution (Migon and Medici
2001, Wang and Ye 2020)
\item \code{'weibull'}, for the weibull distribution (Tsionas 2007)
\item \code{'genexponential'}, for the generalized exponential distribution
(Papadopoulos 2020)
\item \code{'tslaplace'}, for the truncated skewed laplace distribution (Wang
2012).
}}

\item{scaling}{Logical. Only when \code{udist = 'tnormal'} and \code{scaling
  = TRUE} that the scaling property model (Wang and Schmidt 2002) is
estimated. Default = \code{FALSE}. (see section \sQuote{Details}).}

\item{start}{Numeric vector. Optional starting values for the maximum
likelihood (ML) estimation.}

\item{logDepVar}{Logical. Inform whether the dependent variable is logged or
not. Default = \code{TRUE}.}

\item{method}{Character string. Optimization algorithm used for the
estimation. Default = \code{'bfgs'}. 9 algorithms are available:
\itemize{
\item \code{'bfgs'}, for Broyden-Fletcher-Goldfarb-Shanno
(\code{\link[maxLik]{maxBFGS}})
\item \code{'bhhh'}, for Berndt-Hall-Hall-Hausman
(\code{\link[maxLik]{maxBHHH}})
\item \code{'nr'}, for Newton-Raphson (\code{\link[maxLik]{maxNR}})
\item \code{'nm'}, for Nelder-Mead (\code{\link[maxLik]{maxNM}})
\item \code{'ucminf'}, implements a quasi-Newton type with BFGS updating of the
inverse Hessian and soft line search with a trust region type monitoring of
the input to the line search algorithm (\code{\link[ucminf]{ucminf}})
\item \code{'mla'}, for general-purpose optimization based on
Marquardt-Levenberg algorithm (\code{\link[marqLevAlg]{mla}})
\item \code{'sr1'}, for Symmetric Rank 1
(\code{\link[trustOptim]{trust.optim}})
\item \code{'sparse'}, for trust regions and sparse Hessian
(\code{\link[trustOptim]{trust.optim}})
\item \code{'nlminb'}, for optimization using PORT routines
(\code{\link[stats]{nlminb}})
}}

\item{hessianType}{Integer. If \code{1} (Default) analytic hessian is
returned for all the distributions except \code{'gamma'},
\code{'lognormal'} and \code{'weibull'} for which the numeric hessian is
returned. If \code{2} bhhh hessian is estimated (\eqn{g'g}) and if \code{3}
robust hessian is computed (\eqn{H^{-1}GH^{-1}}).}

\item{simType}{Character string. If \code{simType = 'halton'} (Default) Halton draws are use
for maximum simulated likelihood (MSL), if \code{simType = 'ghalton'}
Generalized-Halton draws are use for MSL, if \code{simType = 'sobol'} Sobol
draws are use for MSL, if \code{simType = 'uniform'} uniform draws are use
for MSL. (see section \sQuote{Details}).}

\item{Nsim}{Numeric. Number for draws for MSL.}

\item{prime}{Integer. Prime number considered for Halton and
Generalized-Halton draws. Default = \code{2}.}

\item{burn}{Numeric. Number of the first observations discarded in the case
of the Halton draws. Default = \code{10}.}

\item{antithetics}{Logical. If \code{TRUE} (Default), antithetics counterpart
of the uniform draws is computed. (see section \sQuote{Details}).}

\item{seed}{Numeric. Seed considered for the random draws.}

\item{itermax}{Numeric. Maximum number of iterations allowed for
optimization.}

\item{printInfo}{Logical. Print information during optimization. Default =
\code{FALSE}.}

\item{tol}{Numeric. Convergence tolerance.}

\item{gradtol}{Numeric. Convergence tolerance for gradient.}

\item{stepmax}{Numeric. Step max for \code{ucminf} algorithm.}

\item{qac}{Character. Quadratic Approximation Correction for \code{'bhhh'}
and \code{'nr'} algorithms. If \code{'stephalving'} step is decreased but
the direction is kept, if \code{'marquardt'} the step length is decreased
while also moving closer to the pure gradient direction. Default =
\code{'marquardt'}. See \code{\link[maxLik]{maxBHHH}} and
\code{\link[maxLik]{maxNR}}.}

\item{x}{An object of class \code{sfacross} returned by the function
\code{sfacross}}

\item{digits}{Numeric. Number of digits displayed in values}

\item{...}{Currently ignored.}
}

\value{
\code{sfacross} return a list of class \code{'sfacross'} containing
the following elements: \item{call}{The matched call.}

\item{formula}{The estimated model.}

\item{S}{Integer. Argument \code{'S'} (see above).}

\item{typeSfa}{Character string. \sQuote{Stochastic Production/Profit
Frontier, e = v - u} when \code{S = 1} and \sQuote{Stochastic Cost
Frontier, e = v + u} when \code{S = -1}.}

\item{Nobs}{Numeric. Number of observations used for optimization.}

\item{nXvar}{Numeric. Number of main explanatory variables.}

\item{nmuHvar}{Numeric. If \code{udist = 'tnormal'} or \code{udist =
  'lognormal'}, number of variables explaining heterogeneity in the truncated
mean.}

\item{scaling}{Logical. Argument \code{'scaling'} (see above).}

\item{logDepVar}{Logical. Argument \code{'logDepVar'} (see above).}
\item{nuHvar}{Numeric. Number of variables explaining heteroscedasticity in
the one-sided error term.}

\item{nvHvar}{Numeric. Number of variables explaining heteroscedasticity in
the two-sided error term.} \item{nParm}{Numeric. Total number of parameters
estimated.}

\item{udist}{Character string. Argument \code{'udist'} (see above).}

\item{startVal}{Numeric vector. Starting value for M(S)L estimation.}

\item{dataTable}{Tibble. Data frame containing information on the dependent
and explanatory variables use for optimization along with residuals and
fitted values of the OLS and M(S)L estimations, and the individual
observation log likelihood.}

\item{olsParam}{Numeric vector. OLS estimates.}

\item{olsStder}{Numeric vector. Standard errors of OLS estimates.}

\item{olsSigmasq}{Numeric. Estimated variance of OLS random error.}

\item{olsLoglik}{Numeric. Log likelihood value of OLS estimation.}
\item{olsSkew}{Numeric. Skewness of the residuals of the OLS estimation.}

\item{olsM3Okay}{Logical. Indicating if the residuals of the OLS estimation
have the expected skewness.}

\item{CoelliM3Test}{Numeric vector. Coelli's test for OLS residuals
skewness. See Coelli (1995).}

\item{AgostinoTest}{List. D'Agostino's test for OLS residuals skewness. See
D'agostino and Pearson (1973).}

\item{optType}{Character string. Optimization algorithm used.}

\item{nIter}{Numeric. Number of iterations of the ML estimation.}

\item{optStatus}{Character string. Optimization algorithm termination message.}

\item{startLoglik}{Numeric. Log-likelihood at the starting values.}

\item{mleLoglik}{Numeric. Log-likelihood value of the M(S)L estimation.}

\item{mleParam}{Numeric vector. Parameters obtained from M(S)L estimation.}

\item{gradient}{Numeric vector. Each variable gradient of the M(S)L
estimation.}

\item{gradL_OBS}{Matrix. Each variable individual observation gradient of
the M(S)L estimation.}

\item{gradientNorm}{Numeric. Gradient norm of the M(S)L estimation.}

\item{invHessian}{Matrix. Covariance matrix of the parameters obtained from
the M(S)L estimation.}

\item{hessianType}{Integer. Argument \code{'hessianType'} (see above).}

\item{mleDate}{Character string. Date and time of the estimated model.}

\item{simDist}{Character string. Argument \code{'simDist'} (see above).}

\item{Nsim}{Integer. Argument \code{'Nsim'} (see above).}

\item{FiMat}{Matrix. Matrix of random draws used for MSL.}
}

\details{The stochastic frontier model is defined as:
\deqn{y_i = \alpha + \mathbf{x}'_i\beta + v_i - Su_i}

\deqn{\epsilon_i = v_i -Su_i}.

where \eqn{y} is the output (cost, revenue, profit) and \eqn{x} is the vector
of inputs and other control variables. \code{S = 1} in the case of production
(profit) frontier function and \code{S = -1} in the case of cost frontier
function.

Model is estimated using maximum likelihood (ML) for most distributions
except the gamma, weibull and log-normal distributions for which maximum
simulated likelihood (MSL) is used. For this latter several draws can be
implemented namely Halton, Generalized Halton, Sobol and uniform. In the case
of uniform draws, antithetics can also be computed: first \code{Nsim/2} draws
are obtained then in second the \code{Nsim/2} other draws are obtained as
counterpart of one (\code{1-draw}).

To accommodate heteroscedasticity in the variance parameters of the error
terms, a single part (right) formula can also be specified. To impose the
positivity to these parameters, the variances are modelled as:
\eqn{\sigma^2_u = \exp{(\delta'Z_u)}} or \eqn{\sigma^2_v = \exp{(\phi'Z_v)}}.
In the case of heterogeneity in the truncated mean \eqn{\mu}, it is modelled
as \eqn{\mu=\omega'Z_{mu}}. The scaling property can be applied for the
truncated normal distribution: \eqn{u \sim h(Z_u, \delta)u} where \eqn{u}
follows a truncated normal distribution \eqn{N^+(\tau, \exp{(cu)})}.

In the case of the truncated normal distribution, the convolution of \eqn{u_i}
and \eqn{v_i} is:

\deqn{f(\epsilon_i)=\frac{1}{\sqrt{\sigma_u^2 + \sigma_v^2}}
  \phi\left(\frac{S\epsilon_i + \mu}{\sqrt{\sigma_u^2 + \sigma_v^2}}\right)
  \Phi\left(\frac{\mu_{i*}}{\sigma_*}\right)/\Phi\left(\frac{\mu}{\sigma_u}\right)}

where

\deqn{\mu_{i*}=\frac{\mu\sigma_v^2 - S\epsilon_i\sigma_u^2}{\sigma_u^2 +
  \sigma_v^2}}

and

\deqn{\sigma_*^2 = \frac{\sigma_u^2 \sigma_v^2}{\sigma_u^2 + \sigma_v^2}}

In the case of the half normal distribution the convolution is obtained by
setting \eqn{\mu=0}
}

\note{For the halton draws, we have adapted the code from the \pkg{mlogit} package.}

\examples{
## Using data on fossil fuel fired steam electric power generation plants in
## the United States
# Translog (cost function) half normal with heteroscedasticity
tl_u_h <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
udist = 'hnormal', uhet = ~ regu, data = utility, S = -1, method = 'sparse')
  summary(tl_u_h)

# Translog (cost function) truncated normal with heteroscedasticity
tl_u_t <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
udist = 'tnormal', muhet = ~ regu, data = utility, S = -1, method = 'ucminf')
  summary(tl_u_t)

# Translog (cost function) truncated normal with scaling property
tl_u_ts <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
udist = 'tnormal', muhet = ~ regu, uhet = ~ regu, data = utility, S = -1,
scaling = TRUE, method = 'mla')
  summary(tl_u_ts)

## Data on Philippine rice producers
# Cobb Douglas (production function) generalized exponential, weibull and
# lognormal distributions
cb_p_ge <- sfacross(formula = log(PROD) ~ log(AREA) + log(LABOR) + log(NPK) +
log(OTHER), udist = 'genexponential', data = ricephil, S = 1, method = 'bfgs')
  summary(cb_p_ge)

cb_p_w <- sfacross(formula = log(PROD) ~ log(AREA) + log(LABOR) + log(NPK) +
log(OTHER), udist = 'weibull', data = ricephil, S = 1, method = 'bfgs',
hessianType = 2, simType = 'halton', Nsim = 150)
  summary(cb_p_w)

cb_p_ln <- sfacross(formula = log(PROD) ~ log(AREA) + log(LABOR) + log(NPK) +
log(OTHER), udist = 'lognormal', data = ricephil, S = 1, method = 'bfgs',
hessianType = 2, simType = 'halton', Nsim = 150)
  summary(cb_p_ln)
}

\references{
Aigner, D., Lovell, C. A. K., & Schmidt, P. (1977). Formulation and
estimation of stochastic frontier production function models. \emph{Journal
of econometrics}, \bold{6}(1), 21--37.

Battese, G. E., & Coelli, T. J. (1995). A model for technical inefficiency
effects in a stochastic frontier production function for panel data.
\emph{Empirical Economics}, \bold{20}(2), 325--332.

Caudill, S. B., & Ford, J. M. (1993). Biases in Frontier Estimation Due to
Heteroscedasticity. \emph{Economics Letters}, \bold{41}(1), 17--20.

Caudill, S. B., Ford, J. M., & Gropper, D. M. (1995). Frontier Estimation and
Firm-Specific Inefficiency Measures in the Presence of Heteroscedasticity.
\emph{Journal of Business & Economic Statistics}, \bold{13}(1), 105--111.

Coelli, T. 1995. Estimators and Hypothesis Tests for a Stochastic Frontier
Function - a Monte-Carlo Analysis. \emph{Journal of Productivity Analysis},
\bold{6}:247--268.

D'Agostino, R., and E.S. Pearson. 1973. Tests for departure from normality.
Empirical results for the distributions of \eqn{b_2} and \eqn{\sqrt{b_1}}.
\emph{Biometrika}, \bold{60}:613--622.

Greene, W. H. (2003). Simulated likelihood estimation of the normal-gamma
stochastic frontier function. \emph{Journal of Productivity Analysis},
\bold{19}(2-3), 179--190.

Hadri, K. (1999). Estimation of a doubly heteroscedastic stochastic frontier
cost function. \emph{Journal of Business & Economic Statistics},
\bold{17}(3), 359--363.

Hajargasht, G. (2015). Stochastic frontiers with a Rayleigh distribution.
\emph{Journal of Productivity Analysis}, \bold{44}(2), 199--208.

Huang, C. J., & Liu, J.-T. (1994). Estimation of a non-neutral stochastic
frontier production function. \emph{Journal of Productivity Analysis},
\bold{5}(2), 171--180.

Kumbhakar, S. C., Ghosh, S., & McGuckin, J. T. (1991). A Generalized
Production Frontier Approach for Estimating Determinants of Inefficiency in
U.S. Dairy Farms. \emph{Journal of Business & Economic Statistics},
\bold{9}(3), 279--286.

Li, Q. (1996). Estimating a stochastic production frontier when the adjusted
error is symmetric. \emph{Economics Letters}, \bold{52}(3), 221--228.

Meeusen, W., & Vandenbroeck, J. (1977). Efficiency Estimation from
Cobb-Douglas Production Functions with Composed Error. \emph{International
Economic Review}, \bold{18}(2), 435--445.

Migon, H. S., & Medici, E. V. (2001). Bayesian hierarchical models for
stochastic production frontier.

Nguyen, N. B. (2010). Estimation of technical efficiency in stochastic
frontier analysis. Bowling Green State University.

Papadopoulos, A. (2020). Stochastic frontier models using the generalized
exponential distribution.

Reifschneider, D., & Stevenson, R. (1991). Systematic Departures from the
Frontier: A Framework for the Analysis of Firm Inefficiency.
\emph{International Economic Review}, \bold{32}(3), 715--723.

Stevenson, R. E. (1980). Likelihood Functions for Generalized Stochastic
Frontier Estimation. \emph{Journal of econometrics}, \bold{13}(1), 57--66.

Tsionas, E. G. (2007). Efficiency measurement with the weibull stochastic
frontier. \emph{Oxford Bulletin of Economics and Statistics}, \bold{69}(5),
693--706.

Wang, K., & Ye, X. (2020). Development of alternative stochastic frontier
models for estimating time-space prism vertices. \emph{Transportation}.

Wang, H.J., & Schmidt, P. (2002). One-step and two-step estimation of the
effects of exogenous variables on technical efficiency levels. \emph{Journal
of Productivity Analysis}, \bold{18}:129--144.

Wang, J. (2012). A Normal Truncated Skewed-Laplace Model in Stochastic
Frontier Analysis. Western Kentucky University.
}

\seealso{
\code{\link{summary.sfacross}} for creating and printing summary
results

\code{\link{coef.sfacross}} for extracting coefficients of the estimation

\code{\link{efficiencies.sfacross}} for calculating efficiency estimates

\code{\link{fitted.sfacross}} for obtaining the fitted frontier value

\code{\link{ic.sfacross}} for computing information criterion

\code{\link{logLik.sfacross}} for extracting log-likelihood value of the
estimation

\code{\link{marginal.sfacross}} for calculating marginal effects of
Z-variables on inefficiency

\code{\link{residuals.sfacross}} for extracting residuals of the estimation

\code{\link{vcov.sfacross}} for extracting the variance-covariance matrix
of the coefficients
}

\author{K Hervé Dakpo, Yann Desjeux and Laure Latruffe}

\keyword{cross-section}
\keyword{likelihood}
\keyword{models}
\keyword{optimize}
