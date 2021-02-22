\name{lcmcross}

\alias{lcmcross}
\alias{print.lcmcross}

\title{Latent class stochastic frontier using cross/pooled sectional data}

\usage{
lcmcross(formula, uhet, vhet, thet, data, subset, S = 1L, udist = "hnormal", 
  start = NULL, lcmClasses = 2, logDepVar = TRUE, method = "bfgs", hessianType = 1,
  itermax = 2000L, printInfo = FALSE, tol = 1e-12, gradtol = 1e-06, stepmax = 0.1, 
  qac = "marquardt", initStart = FALSE, initAlg = "nlminb", initIter = 100,
  initFactorLB = 0.5, initFactorUB = 1.5)

\method{print}{lcmcross}(x, digits = max(3, getOption("digits") - 2), ...)
}

\arguments{
\item{formula}{A symbolic description of the model to be estimated based on
the generic function \code{formula} (see section \sQuote{Details}).}

\item{uhet}{A one-part formula to consider heteroscedasticity in the
one-sided error variance (see section \sQuote{Details}).}

\item{vhet}{A one-part formula to consider heteroscedasticity in the
two-sided error variance (see section \sQuote{Details}).}

\item{thet}{A one-part formula to consider technological heterogeneity in the
construction of the classes.}

\item{data}{An optional data frame containing the data or the the variables
in the model.}

\item{subset}{An optional vector specifying a subset of observations to be
used in the optimization process.}

\item{S}{Integer. If \code{1} (Default) a production (profit) frontier is
estimated \eqn{\epsilon_i = v_i-u_i}. If \code{-1} a cost frontier is
estimated \eqn{\epsilon_i = v_i+u_i}.}

\item{udist}{Character string. Distribution specification for the one-sided
error term. Only \code{"hnormal"} for the half normal distribution (Aigner
\emph{et al.} 1977, Meeusen and Vandenbroeck 1977) is is currently
implemented.}

\item{start}{Numeric vector. Optional starting values for the maximum
likelihood (ML) estimation.}

\item{lcmClasses}{Numeric. Number of classes to be estimated. Maximum of five
classes can be estimated. Default = \code{2}.}

\item{logDepVar}{Logical. Inform whether the dependent variable is logged or
not. Default = \code{TRUE}.}

\item{method}{Character. Optimization algorithm used for the estimation.
Default = \code{"bfgs"}. 9 algorithms are available:
\itemize{
\item \code{"bfgs"}, for Broyden-Fletcher-Goldfarb-Shanno
(\code{\link[maxLik]{maxBFGS}})
\item \code{"bhhh"}, for Berndt-Hall-Hall-Hausman
(\code{\link[maxLik]{maxBHHH}})
\item \code{"nr"}, for Newton-Raphson (\code{\link[maxLik]{maxNR}})
\item \code{"nm"}, for Nelder-Mead (\code{\link[maxLik]{maxNM}})
\item \code{"ucminf"}, implements a quasi-Newton type with BFGS updating of the
inverse Hessian and soft line search with a trust region type monitoring of
the input to the line search algorithm (\code{\link[ucminf]{ucminf}})
\item \code{"mla"}, for general-purpose optimization based on
Marquardt-Levenberg algorithm (\code{\link[marqLevAlg]{mla}})
\item \code{"sr1"}, for Symmetric Rank 1
(\code{\link[trustOptim]{trust.optim}})
\item \code{"sparse"}, for trust regions and sparse Hessian
(\code{\link[trustOptim]{trust.optim}})
\item \code{"nlminb"}, for optimization using PORT routines
(\code{\link[stats]{nlminb}})
}}

\item{hessianType}{Integer. If \code{1} (Default) analytic hessian is
returned for all the distributions except \code{"gamma"},
\code{"lognormal"} and \code{"weibull"} for which the numeric hessian is
returned. If \code{2} bhhh hessian is estimated (\eqn{g'g}) and if \code{3}
robust hessian is computed (\eqn{H^{-1}GH^{-1}}).}

\item{itermax}{Numeric. Maximum number of iterations allowed for
optimization.}

\item{printInfo}{Logical. Print information during optimization. Default =
\code{FALSE}.}

\item{tol}{Numeric. Convergence tolerance.}

\item{gradtol}{Numeric. Convergence tolerance for gradient.}

\item{stepmax}{Numeric. Step max for \code{ucminf} algorithm.}

\item{qac}{Character. Quadratic Approximation Correction for \code{"bhhh"}
and \code{"nr"} algorithms. If \code{"stephalving"} step is decreased but
the direction is kept, if \code{"marquardt"} the step length is decreased
while also moving closer to the pure gradient direction. Default =
\code{"marquardt"}. See \code{\link[maxLik]{maxBHHH}} and
\code{\link[maxLik]{maxNR}}.}

\item{initStart}{Logical. if \code{TRUE} the model is jump-started using an
alternative algorithm (\code{"nlminb"}) within certain bounds. Default =
\code{FALSE}.}

\item{initAlg}{Character. Algorithm used to jump-start the latent class
model. Only \code{"nlminb"} is currently available.}

\item{initIter}{Numeric. Maximum number of iterations for the algorihtm when
\code{initStart} is \code{TRUE}.}

\item{initFactorLB}{Numeric. Value indicating by which factor the
starting value should be multiplied by to define the lower (or upper)
bounds  for the jump-start algorithm. Default = \code{0.5}.}

\item{initFactorUB}{Numeric. Value indicating by which factor the
starting value should be multiplied by to define the upper (or lower)
bounds  for the jump-start algorithm. Default = \code{1.5}.}

\item{x}{An object of class \code{lcmcross} returned by the function
\code{lcmcross}}

\item{digits}{Numeric. Number of digits displayed in values}

\item{...}{Currently ignored.}
}
\value{
\code{lcmcross} return a list of class \code{'lcmcross'} containing
the following elements:

\item{call}{The matched call.}

\item{formula}{Multi parts formula describing the estimated model.}

\item{S}{Integer. Argument \code{'S'} (see above).}

\item{typeSfa}{Character string. \sQuote{Latent Class Production/Profit
Frontier, e = v - u} when \code{S = 1} and \sQuote{Latent Class
Production/Profit Frontier, e = v + u} when \code{S = -1}.}

\item{Nobs}{Numeric. Number of observations used for optimization.}

\item{nXvar}{Numeric. Number of main explanatory variables.}

\item{nZvar}{Numeric. Number of variables in the logit specification of the
finite mixture model.}

\item{logDepVar}{Logical. Argument \code{'logDepVar'} (see above).}

\item{nuHvar}{Numeric. Number of variables explaining heteroscedasticity in
the one-sided error term.}

\item{nvHvar}{Numeric. Number of variables explaining heteroscedasticity in
the two-sided error term.}

\item{nParm}{Numeric. Total number of parameters estimated.}

\item{udist}{Character. Argument \code{'udist'} (see above).}

\item{startVal}{Numeric vector. Starting value for ML estimation.}

\item{dataTable}{Tibble. Data frame containing information on the dependent
and explanatory variables use for optimization along with residuals and
fitted values of the OLS and ML estimations, and the individual observation
log likelihood.}

\item{InitHalf}{Initial ML estimation with half normal distribution for the
one-sided error term. Model to construct the starting values for the latent
class estimation. Object of class \code{'maxLik'} and \code{'maxim'}
returned.}
\item{optType}{Character. Optimization algorithm used.}

\item{nIter}{Numeric. Number of iterations of the ML estimation.}

\item{optStatus}{Character. Optimization algorithm termination message.}

\item{startLoglik}{Numeric. Log-likelihood at the starting values.}

\item{nClasses}{Numeric. Number of classes estimated.}

\item{mleLoglik}{Numeric. Log-likelihood value of the ML estimation.}

\item{mleParam}{Numeric vector. Parameters obtained from ML estimation.}

\item{gradient}{Numeric vector. Each variable gradient of the ML
estimation.}

\item{gradL_OBS}{Matrix. Each variable individual observation gradient of
the ML estimation.}

\item{gradientNorm}{Numeric. Gradient norm of the ML estimation.}

\item{invHessian}{Matrix. Covariance matrix of the parameters obtained from
the ML estimation.}

\item{hessianType}{Integer. Argument \code{'hessianType'} (see above).}

\item{mleDate}{Character. Date and time of the estimated model.}
}
\description{
Function \code{lcmcross} is a symbolic formula based function for the
estimation of the latent class stochastic frontier model (LCM) in the case of
cross-sectional or pooled cross section data. Model is estimated using
maximum likelihood (ML). See Orea and Kumbhakar (2004), Parmeter and
Kumbhakar (2014, p282)

Only the half-normal distribution is considered for the one-sided error term.
Nine optimization algorithms are also available.

The function also accounts for heteroscedasticity in both one-sided and
two-sided error terms, as in Reifschneider and Stevenson (1991), Caudill and
Ford (1993), Caudill \emph{et al.} (1995) and Hadri (1999).

Finally, the model can estimate up to five classes.
}
\details{
LCM is an estimation of a finite mixture of production functions:

\deqn{y_i = \alpha_j + x'_i\beta_j + v_{i|j} - Su_{i|j}}

\deqn{\epsilon_{i|j} = v_{i|j} -Su_{i|j}}

where \eqn{y} is the output (cost, revenue, profit) and \eqn{x} is the vector
of inputs and other control variables. \code{S = 1} in the case of production
(profit) frontier function and \code{S = -1} in the case of cost frontier
function.

The contribution of observation i to the likelihood conditional on class j
is defined as

\deqn{P(i|j) = \frac{2}{\sqrt{\sigma_{u|j}^2 + \sigma_{v|j}^2}}
  \phi\left(\frac{S\epsilon_{i|j}}{\sqrt{\sigma_{u|j}^2 +
  \sigma_{v|j}^2}}\right) \Phi\left(\frac{\mu_{i*|j}}{\sigma_{*|j}}\right)}

where

\deqn{\mu_{i*|j}=\frac{- S\epsilon_{i|j}\sigma_{u|j}^2}{\sigma_{u|j}^2 +
  \sigma_{v|j}^2}}

and

\deqn{\sigma_*^2 = \frac{\sigma_{u|j}^2 \sigma_{v|j}^2}{\sigma_{u|j}^2 +
  \sigma_{v|j}^2}}

The prior probability of using a particular technology can depend on some
covariates (namely the variables separating the observations into classes)
using a logit specification:

\deqn{\pi(i,j) = \frac{\exp{(\theta_j'Z_h)}}{\sum_{m=1}^{J}\exp{(\theta_m'Z_h)}}}

with \eqn{\exp(\theta_J'Z_h)=1}.

In the case of panel data, \code{lcmcross} estimates a pooled cross section
where the probability of beloging to a class a priori is not permanent
(fixed) over time.

The unconditional likelihood of observation i is simply the average over
the J classes:

\eqn{P(i) = \sum_{m=1}^{J}\pi(i,m)P(i|m)}

The number of classes can be retained based on information criterion (see
for instance \code{\link{ic.lcmcross}}).

To accommodate heteroscedasticity in the variance parameters of the error
terms, a single part (right) formula can also be specified. To impose the
positivity on these parameters, the variances are modelled as:
\eqn{\sigma^2_{u|j} = \exp{(\delta_j'Z_u)}} or \eqn{\sigma^2_{v|j} =
  \exp{(\phi_j'Z_v)}}.
}

\examples{
## Using data on eighty-two countries production (DGP)
# Cobb Douglas (production function) half normal distribution
# Only the intercept is used as the separating variables
cb_2c_h <- lcmcross(formula = ly ~ lk + ll + yr, udist = "hnormal", 
           data = worldprod, S = 1, method = "ucminf")
  summary(cb_2c_h)

# summary of the initial MLE model
  summary(cb_2c_h$InitHalf)

# same result by jump-starting the estimation
cb_2c_h <- lcmcross(formula = ly ~ lk + ll + yr, udist = "hnormal", 
           data = worldprod, S = 1, method = "bfgs", initStart = TRUE)
  summary(cb_2c_h)

# Intercept and initStat used as separating variables
cb_2c_h <- lcmcross(formula = ly ~ lk + ll + yr, udist = "hnormal", 
           thet = ~initStat, data = worldprod, S = 1, method = "ucminf")
  summary(cb_2c_h)

# Only the intercept is used as the separating variables + inefficiency
# determinants
cb_2c_h <- lcmcross(formula = ly ~ lk + ll + yr, udist = "hnormal", 
           uhet = ~initStat, data = worldprod, S = 1, method = "ucminf")
  summary(cb_2c_h)
}

\references{
Aigner, D., Lovell, C. A. K., and P. Schmidt. (1977). Formulation and
estimation of stochastic frontier production function models. \emph{Journal
of econometrics}, \bold{6}(1), 21--37.

Caudill, S. B., and J. M. Ford. (1993). Biases in Frontier Estimation Due to
Heteroscedasticity. \emph{Economics Letters}, \bold{41}(1), 17--20.

Caudill, S. B., Ford, J. M., and D. M. Gropper. (1995). Frontier Estimation
and Firm-Specific Inefficiency Measures in the Presence of
Heteroscedasticity. \emph{Journal of Business & Economic Statistics},
\bold{13}(1), 105--111.

Hadri, K. (1999). Estimation of a doubly heteroscedastic stochastic frontier
cost function. \emph{Journal of Business & Economic Statistics},
\bold{17}(3), 359--363.

Meeusen, W., and J. Vandenbroeck. (1977). Efficiency Estimation from
Cobb-Douglas Production Functions with Composed Error. \emph{International
Economic Review}, \bold{18}(2), 435--445.

Orea, L., and S.C. Kumbhakar. (2004). Efficiency measurement using a latent
class stochastic frontier model. \emph{Empirical Economics}, \bold{29},
169--183.

Parmeter, C.F., and S.C. Kumbhakar. (2014). Efficiency analysis: A primer on
recent advances. \emph{Foundations and Trends in Econometrics}, \bold{7},
191--385.

Reifschneider, D., and R. Stevenson. (1991). Systematic Departures from the
Frontier: A Framework for the Analysis of Firm Inefficiency.
\emph{International Economic Review}, \bold{32}(3), 715--723.
}
\seealso{
\code{\link{summary.lcmcross}} for creating and printing summary results

\code{\link{coef.lcmcross}} for extracting coefficients of the estimation

\code{\link{efficiencies.lcmcross}} for calculating efficiency estimates

\code{\link{fitted.lcmcross}} for obtaining the fitted frontier value

\code{\link{ic.lcmcross}} for computing information criterion

\code{\link{logLik.lcmcross}} for extracting log-likelihood value of the estimation

\code{\link{marginal.lcmcross}} for calculating marginal effects of
Z-variables on inefficiency

\code{\link{residuals.lcmcross}} for extracting residuals of the estimation

\code{\link{vcov.lcmcross}} for extracting the variance-covariance matrix
of the coefficients
}

\author{K Hervé Dakpo, Yann Desjeux and Laure Latruffe}

\keyword{cross-section}
\keyword{latent}
\keyword{likelihood}
\keyword{models}
\keyword{optimize}
