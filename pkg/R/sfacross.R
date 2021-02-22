# SFA estimation for cross sectional data ----------
#' @name sfacross
#'
#' @aliases sfacross
#'
#' @title Stochastic frontier estimation using cross/pooled sectional data
#'
#' @description
#' \code{sfacross} is a symbolic formula-based function for the estimation of
#' stochastic frontier models in the case of cross-sectional or pooled cross
#' section data, using maximum (simulated) likelihood - M(S)L.
#'
#' Ten distributions are considered for the one-sided error term and nine
#' optimization algorithms are also available.
#'
#' The function also accounts for heteroscedasticity in both one-sided and
#' two-sided error terms as in Reifschneider and Stevenson (1991), Caudill and
#' Ford (1993), Caudill \emph{et al.} (1995) and Hadri (1999), but also
#' heterogeneity in the mean of the pre-truncated distribution as in Kumbhakar
#' \emph{et al.} (1991), Huang and Liu (1994), Battese and Coelli (1995).
#'
#' Finally, the truncated normal - normal distribution with scaling property as
#' in Wang and Schmidt (2002) is also implemented.
#'
#' @param formula A symbolic description of the model to be estimated based on
#'   the generic function \code{formula} (see section \sQuote{Details}).
#' @param muhet A one-part formula to consider heterogeneity in the mean of the
#'   pre-truncated distribution (see section \sQuote{Details}).
#' @param uhet A one-part formula to consider heteroscedasticity in the
#'   one-sided error variance (see section \sQuote{Details}).
#' @param vhet A one-part formula to consider heteroscedasticity in the
#'   two-sided error variance (see section \sQuote{Details}).
#' @param data An optional data frame containing the data or the the variables
#'   in the model.
#' @param subset An optional vector specifying a subset of observations to be
#'   used in the optimization process.
#' @param S Integer. If \code{1} (Default) a production (profit) frontier is
#'   estimated \eqn{\epsilon_i = v_i-u_i}. If \code{-1} a cost frontier is
#'   estimated \eqn{\epsilon_i = v_i+u_i}.
#' @param udist Character string. Distribution specification for the one-sided error
#'   term. Default = \code{"hnormal"}. 10 different distributions are available:
#'
#'   - \code{'hnormal'}, for the half normal distribution (Aigner \emph{et al.}
#'   1977, Meeusen and Vandenbroeck 1977)
#'
#'   - \code{'exponential'}, for the exponential distribution
#'
#'   - \code{tnormal} for the truncated normal distribution (Stevenson 1980)
#'
#'   - \code{'rayleigh'}, for the rayleigh distribution (Hajargasht 2015)
#'
#'   - \code{'uniform'}, for the uniform distribution (Li 1996, Nguyen 2010)
#'
#'   - \code{'gamma'}, for the gamma distribution (Greene 2003)
#'
#'   - \code{'lognormal'}, for the log normal distribution (Migon and Medici
#'   2001, Wang and Ye 2020)
#'
#'   - \code{'weibull'}, for the weibull distribution (Tsionas 2007)
#'
#'   - \code{'genexponential'}, for the generalized exponential distribution
#'   (Papadopoulos 2020)
#'
#'   - \code{'tslaplace'}, for the truncated skewed laplace distribution (Wang
#'   2012).
#'
#' @param scaling Logical. Only when \code{udist = 'tnormal'} and \code{scaling
#'   = TRUE} that the scaling property model (Wang and Schmidt 2002) is
#'   estimated. Default = \code{FALSE}. (see section \sQuote{Details}).
#' @param start Numeric vector. Optional starting values for the maximum
#'   likelihood (ML) estimation.
#' @param logDepVar Logical. Inform whether the dependent variable is logged or
#'   not. Default = \code{TRUE}.
#' @param method Character string. Optimization algorithm used for the
#'   estimation. Default = \code{'bfgs'}. 9 algorithms are available:
#'
#'   - \code{'bfgs'}, for Broyden-Fletcher-Goldfarb-Shanno
#'   (\code{\link[maxLik]{maxBFGS}})
#'
#'   - \code{'bhhh'}, for Berndt-Hall-Hall-Hausman
#'   (\code{\link[maxLik]{maxBHHH}})
#'
#'   - \code{'nr'}, for Newton-Raphson (\code{\link[maxLik]{maxNR}})
#'
#'   - \code{'nm'}, for Nelder-Mead (\code{\link[maxLik]{maxNM}})
#'
#'   - \code{'ucminf'}, implements a quasi-Newton type with BFGS updating of the
#'   inverse Hessian and soft line search with a trust region type monitoring of
#'   the input to the line search algorithm (\code{\link[ucminf]{ucminf}})
#'
#'   - \code{'mla'}, for general-purpose optimization based on
#'   Marquardt-Levenberg algorithm (\code{\link[marqLevAlg]{mla}})
#'
#'   - \code{'sr1'}, for Symmetric Rank 1
#'   (\code{\link[trustOptim]{trust.optim}})
#'
#'   - \code{'sparse'}, for trust regions and sparse Hessian
#'   (\code{\link[trustOptim]{trust.optim}})
#'
#'   - \code{'nlminb'}, for optimization using PORT routines
#'   (\code{\link[stats]{nlminb}})
#' @param hessianType Integer. If \code{1} (Default) analytic hessian is
#'   returned for all the distributions except \code{'gamma'},
#'   \code{'lognormal'} and \code{'weibull'} for which the numeric hessian is
#'   returned. If \code{2} bhhh hessian is estimated (\eqn{g'g}) and if \code{3}
#'   robust hessian is computed (\eqn{H^{-1}GH^{-1}}).
#' @param simType Character string. If \code{simType = 'halton'} (Default) Halton draws are use
#'   for maximum simulated likelihood (MSL), if \code{simType = 'ghalton'}
#'   Generalized-Halton draws are use for MSL, if \code{simType = 'sobol'} Sobol
#'   draws are use for MSL, if \code{simType = 'uniform'} uniform draws are use
#'   for MSL. (see section \sQuote{Details}).
#' @param Nsim Numeric. Number for draws for MSL.
#' @param prime Integer. Prime number considered for Halton and
#'   Generalized-Halton draws. Default = \code{2}.
#' @param burn Numeric. Number of the first observations discarded in the case
#'   of the Halton draws. Default = \code{10}.
#' @param antithetics Logical. If \code{TRUE} (Default), antithetics counterpart
#'   of the uniform draws is computed. (see section \sQuote{Details}).
#' @param seed Numeric. Seed considered for the random draws.
#' @param itermax Numeric. Maximum number of iterations allowed for
#'   optimization.
#' @param printInfo Logical. Print information during optimization. Default =
#'   \code{FALSE}.
#' @param tol Numeric. Convergence tolerance.
#' @param gradtol Numeric. Convergence tolerance for gradient.
#' @param stepmax Numeric. Step max for \code{ucminf} algorithm.
#' @param qac Character. Quadratic Approximation Correction for \code{'bhhh'}
#'   and \code{'nr'} algorithms. If \code{'stephalving'} step is decreased but
#'   the direction is kept, if \code{'marquardt'} the step length is decreased
#'   while also moving closer to the pure gradient direction. Default =
#'   \code{'marquardt'}. See \code{\link[maxLik]{maxBHHH}} and
#'   \code{\link[maxLik]{maxNR}}.
#' @param x An object of class \code{sfacross} returned by the function
#'   \code{sfacross}
#' @param digits Numeric. Number of digits displayed in values
#' @param ... Currently ignored.
#'
#' @details
#' The stochastic frontier model is defined as:
#'   \deqn{y_i = \alpha + \mathbf{x}'_i\beta + v_i - Su_i}
#'
#'   \deqn{\epsilon_i = v_i -Su_i}.
#'
#' where \eqn{y} is the output (cost, revenue, profit) and \eqn{x} is the vector
#' of inputs and other control variables. \code{S = 1} in the case of production
#' (profit) frontier function and \code{S = -1} in the case of cost frontier
#' function.
#'
#' Model is estimated using maximum likelihood (ML) for most distributions
#' except the gamma, weibull and log-normal distributions for which maximum
#' simulated likelihood (MSL) is used. For this latter several draws can be
#' implemented namely Halton, Generalized Halton, Sobol and uniform. In the case
#' of uniform draws, antithetics can also be computed: first \code{Nsim/2} draws
#' are obtained then in second the \code{Nsim/2} other draws are obtained as
#' counterpart of one (\code{1-draw}).
#'
#' To accommodate heteroscedasticity in the variance parameters of the error
#' terms, a single part (right) formula can also be specified. To impose the
#' positivity to these parameters, the variances are modelled as:
#' \eqn{\sigma^2_u = \exp{(\delta'Z_u)}} or \eqn{\sigma^2_v = \exp{(\phi'Z_v)}}.
#' In the case of heterogeneity in the truncated mean \eqn{\mu}, it is modelled
#' as \eqn{\mu=\omega'Z_{mu}}. The scaling property can be applied for the
#' truncated normal distribution: \eqn{u \sim h(Z_u, \delta)u} where \eqn{u}
#' follows a truncated normal distribution \eqn{N^+(\tau, \exp{(cu)})}.
#'
#' In the case of the truncated normal distribution, the convolution of \eqn{u_i}
#'   and \eqn{v_i} is:
#'
#'   \deqn{f(\epsilon_i)=\frac{1}{\sqrt{\sigma_u^2 + \sigma_v^2}}
#'   \phi\left(\frac{S\epsilon_i + \mu}{\sqrt{\sigma_u^2 + \sigma_v^2}}\right)
#'   \Phi\left(\frac{\mu_{i*}}{\sigma_*}\right)/\Phi\left(\frac{\mu}{\sigma_u}\right)}
#'
#'   where
#'
#'   \deqn{\mu_{i*}=\frac{\mu\sigma_v^2 - S\epsilon_i\sigma_u^2}{\sigma_u^2 +
#'   \sigma_v^2}}
#'
#'   and
#'
#'   \deqn{\sigma_*^2 = \frac{\sigma_u^2 \sigma_v^2}{\sigma_u^2 + \sigma_v^2}}
#'
#'   In the case of the half normal distribution the convolution is obtained by
#'   setting \eqn{\mu=0}
#'
#' @return \code{sfacross} return a list of class \code{'sfacross'} containing
#'   the following elements: \item{call}{The matched call.}
#'
#'   \item{formula}{The estimated model.}
#'
#'   \item{S}{Integer. Argument \code{'S'} (see above).}
#'
#'   \item{typeSfa}{Character string. \sQuote{Stochastic Production/Profit
#'   Frontier, e = v - u} when \code{S = 1} and \sQuote{Stochastic Cost
#'   Frontier, e = v + u} when \code{S = -1}.}
#'
#'   \item{Nobs}{Numeric. Number of observations used for optimization.}
#'
#'   \item{nXvar}{Numeric. Number of main explanatory variables.}
#'
#'   \item{nmuHvar}{Numeric. If \code{udist = 'tnormal'} or \code{udist =
#'   'lognormal'}, number of variables explaining heterogeneity in the truncated
#'   mean.}
#'
#'   \item{scaling}{Logical. Argument \code{'scaling'} (see above).}
#'
#'   \item{logDepVar}{Logical. Argument \code{'logDepVar'} (see above).}
#'   \item{nuHvar}{Numeric. Number of variables explaining heteroscedasticity in
#'   the one-sided error term.}
#'
#'   \item{nvHvar}{Numeric. Number of variables explaining heteroscedasticity in
#'   the two-sided error term.} \item{nParm}{Numeric. Total number of parameters
#'   estimated.}
#'
#'   \item{udist}{Character string. Argument \code{'udist'} (see above).}
#'
#'   \item{startVal}{Numeric vector. Starting value for M(S)L estimation.}
#'
#'   \item{dataTable}{Tibble. Data frame containing information on the dependent
#'   and explanatory variables use for optimization along with residuals and
#'   fitted values of the OLS and M(S)L estimations, and the individual
#'   observation log likelihood.}
#'
#'   \item{olsParam}{Numeric vector. OLS estimates.}
#'
#'   \item{olsStder}{Numeric vector. Standard errors of OLS estimates.}
#'
#'   \item{olsSigmasq}{Numeric. Estimated variance of OLS random error.}
#'
#'   \item{olsLoglik}{Numeric. Log likelihood value of OLS estimation.}
#'   \item{olsSkew}{Numeric. Skewness of the residuals of the OLS estimation.}
#'
#'   \item{olsM3Okay}{Logical. Indicating if the residuals of the OLS estimation
#'   have the expected skewness.}
#'
#'   \item{CoelliM3Test}{Numeric vector. Coelli's test for OLS residuals
#'   skewness. See Coelli (1995).}
#'
#'   \item{AgostinoTest}{List. D'Agostino's test for OLS residuals skewness. See
#'   D'agostino and Pearson (1973).}
#'
#'   \item{optType}{Character string. Optimization algorithm used.}
#'
#'   \item{nIter}{Numeric. Number of iterations of the ML estimation.}
#'
#'   \item{optStatus}{Character string. Optimization algorithm termination message.}
#'
#'   \item{startLoglik}{Numeric. Log-likelihood at the starting values.}
#'
#'   \item{mleLoglik}{Numeric. Log-likelihood value of the M(S)L estimation.}
#'
#'   \item{mleParam}{Numeric vector. Parameters obtained from M(S)L estimation.}
#'
#'   \item{gradient}{Numeric vector. Each variable gradient of the M(S)L
#'   estimation.}
#'
#'   \item{gradL_OBS}{Matrix. Each variable individual observation gradient of
#'   the M(S)L estimation.}
#'
#'   \item{gradientNorm}{Numeric. Gradient norm of the M(S)L estimation.}
#'
#'   \item{invHessian}{Matrix. Covariance matrix of the parameters obtained from
#'   the M(S)L estimation.}
#'
#'   \item{hessianType}{Integer. Argument \code{'hessianType'} (see above).}
#'
#'   \item{mleDate}{Character string. Date and time of the estimated model.}
#'
#'   \item{simDist}{Character string. Argument \code{'simDist'} (see above).}
#'
#'   \item{Nsim}{Integer. Argument \code{'Nsim'} (see above).}
#'
#'   \item{FiMat}{Matrix. Matrix of random draws used for MSL.}
#'
#' @export
#'
#' @author K Hervé Dakpo, Yann Desjeux and Laure Latruffe
#'
#' @references
#'
#' Aigner, D., Lovell, C. A. K., & Schmidt, P. (1977). Formulation and
#' estimation of stochastic frontier production function models. \emph{Journal
#' of econometrics}, \bold{6}(1), 21--37.
#'
#' Battese, G. E., & Coelli, T. J. (1995). A model for technical inefficiency
#' effects in a stochastic frontier production function for panel data.
#' \emph{Empirical Economics}, \bold{20}(2), 325--332.
#'
#' Caudill, S. B., & Ford, J. M. (1993). Biases in Frontier Estimation Due to
#' Heteroscedasticity. \emph{Economics Letters}, \bold{41}(1), 17--20.
#'
#' Caudill, S. B., Ford, J. M., & Gropper, D. M. (1995). Frontier Estimation and
#' Firm-Specific Inefficiency Measures in the Presence of Heteroscedasticity.
#' \emph{Journal of Business & Economic Statistics}, \bold{13}(1), 105--111.
#'
#' Coelli, T. 1995. Estimators and Hypothesis Tests for a Stochastic Frontier
#' Function - a Monte-Carlo Analysis. \emph{Journal of Productivity Analysis},
#' \bold{6}:247-268.
#'
#' D'Agostino, R., and E.S. Pearson. 1973. Tests for departure from normality.
#' Empirical results for the distributions of \eqn{b_2} and \eqn{\sqrt{b_1}}.
#' \emph{Biometrika}, \bold{60}:613-622.
#'
#' Greene, W. H. (2003). Simulated likelihood estimation of the normal-gamma
#' stochastic frontier function. \emph{Journal of Productivity Analysis},
#' \bold{19}(2-3), 179--190.
#'
#' Hadri, K. (1999). Estimation of a doubly heteroscedastic stochastic frontier
#' cost function. \emph{Journal of Business & Economic Statistics},
#' \bold{17}(3), 359--363.
#'
#' Hajargasht, G. (2015). Stochastic frontiers with a Rayleigh distribution.
#' \emph{Journal of Productivity Analysis}, \bold{44}(2), 199--208.
#'
#' Huang, C. J., & Liu, J.-T. (1994). Estimation of a non-neutral stochastic
#' frontier production function. \emph{Journal of Productivity Analysis},
#' \bold{5}(2), 171--180.
#'
#' Kumbhakar, S. C., Ghosh, S., & McGuckin, J. T. (1991). A Generalized
#' Production Frontier Approach for Estimating Determinants of Inefficiency in
#' U.S. Dairy Farms. \emph{Journal of Business & Economic Statistics},
#' \bold{9}(3), 279--286.
#'
#' Li, Q. (1996). Estimating a stochastic production frontier when the adjusted
#' error is symmetric. \emph{Economics Letters}, \bold{52}(3), 221--228.
#'
#' Meeusen, W., & Vandenbroeck, J. (1977). Efficiency Estimation from
#' Cobb-Douglas Production Functions with Composed Error. \emph{International
#' Economic Review}, \bold{18}(2), 435--445.
#'
#' Migon, H. S., & Medici, E. V. (2001). Bayesian hierarchical models for
#' stochastic production frontier.
#'
#' Nguyen, N. B. (2010). Estimation of technical efficiency in stochastic
#' frontier analysis. Bowling Green State University.
#'
#' Papadopoulos, A. (2020). Stochastic frontier models using the generalized
#' exponential distribution.
#'
#' Reifschneider, D., & Stevenson, R. (1991). Systematic Departures from the
#' Frontier: A Framework for the Analysis of Firm Inefficiency.
#' \emph{International Economic Review}, \bold{32}(3), 715--723.
#'
#' Stevenson, R. E. (1980). Likelihood Functions for Generalized Stochastic
#' Frontier Estimation. \emph{Journal of econometrics}, \bold{13}(1), 57--66.
#'
#' Tsionas, E. G. (2007). Efficiency measurement with the weibull stochastic
#' frontier. \emph{Oxford Bulletin of Economics and Statistics}, \bold{69}(5),
#' 693--706.
#'
#' Wang, K., & Ye, X. (2020). Development of alternative stochastic frontier
#' models for estimating time-space prism vertices. \emph{Transportation}.
#'
#' Wang, H.J., & Schmidt, P. (2002). One-step and two-step estimation of the
#' effects of exogenous variables on technical efficiency levels. \emph{Journal
#' of Productivity Analysis}, \bold{18}:129--144.
#'
#' Wang, J. (2012). A Normal Truncated Skewed-Laplace Model in Stochastic
#' Frontier Analysis. Western Kentucky University.
#'
#' @seealso
#'
#'   \code{\link{summary.sfacross}} for creating and printing summary
#'   results
#'
#'   \code{\link{coef.sfacross}} for extracting coefficients of the estimation
#'
#'   \code{\link{efficiencies.sfacross}} for calculating efficiency estimates
#'
#'   \code{\link{fitted.sfacross}} for obtaining the fitted frontier value
#'
#'   \code{\link{ic.sfacross}} for computing information criterion
#'
#'   \code{\link{logLik.sfacross}} for extracting log-likelihood value of the
#'   estimation
#'
#'   \code{\link{marginal.sfacross}} for calculating marginal effects of
#'   Z-variables on inefficiency
#'
#'   \code{\link{residuals.sfacross}} for extracting residuals of the estimation
#'
#'   \code{\link{vcov.sfacross}} for extracting the variance-covariance matrix
#'   of the coefficients
#'
#' @note For the halton draws, we have adapted the code from
#'   \code{mlogit} package.
#'
#' @examples
#' ## Using data on fossil fuel fired steam electric power generation plants in
#' ## the United States
#'
#' # Translog (cost function) half normal with heteroscedasticity
#'
#' tl_u_h <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
#' log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
#' I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
#' udist = 'hnormal', uhet = ~ regu, data = utility, S = -1, method = 'sparse')
#'
#' summary(tl_u_h)
#'
#' # Translog (cost function) truncated normal with heteroscedasticity
#'
#' tl_u_t <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
#' log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
#' I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
#' udist = 'tnormal', muhet = ~ regu, data = utility, S = -1, method = 'ucminf')
#'
#' summary(tl_u_t)
#'
#' # Translog (cost function) truncated normal with scaling property
#'
#' tl_u_ts <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
#' log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
#' I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
#' udist = 'tnormal', muhet = ~ regu, uhet = ~ regu, data = utility, S = -1,
#' scaling = TRUE, method = 'mla')
#'
#' summary(tl_u_ts)
#'
#' ## Data on Philippine rice producers
#'
#' # Cobb Douglas (production function) generalized exponential, weibull and
#' # lognormal distributions
#'
#' cb_p_ge <- sfacross(formula = log(PROD) ~ log(AREA) + log(LABOR) + log(NPK) +
#' log(OTHER), udist = 'genexponential', data = ricephil, S = 1, method =
#' 'bfgs')
#'
#' summary(cb_p_ge)
#'
#' cb_p_w <- sfacross(formula = log(PROD) ~ log(AREA) + log(LABOR) + log(NPK) +
#' log(OTHER), udist = 'weibull', data = ricephil, S = 1, method = 'bfgs',
#' hessianType = 2, simType = 'halton', Nsim = 150)
#'
#' summary(cb_p_w)
#'
#' cb_p_ln <- sfacross(formula = log(PROD) ~ log(AREA) + log(LABOR) + log(NPK) +
#' log(OTHER), udist = 'lognormal', data = ricephil, S = 1, method = 'bfgs',
#' hessianType = 2, simType = 'halton', Nsim = 150)
#'
#' summary(cb_p_ln)
#'
#' @keywords models optimize cross-section likelihood
sfacross <- function(formula, muhet, uhet, vhet, data, subset,
                     S = 1L, udist = "hnormal", scaling = FALSE, start = NULL,
                     logDepVar = TRUE, method = "bfgs", hessianType = 1L, simType = "halton",
                     Nsim = 100, prime = 2L, burn = 10, antithetics = FALSE, seed = 12345,
                     itermax = 2000, printInfo = FALSE, tol = 1e-12, gradtol = 1e-06,
                     stepmax = 0.1, qac = "marquardt") {
  # u distribution check -------
  udist <- tolower(udist)
  if (!(udist %in% c(
    "hnormal", "exponential", "tnormal", "rayleigh",
    "uniform", "gamma", "lognormal", "weibull", "genexponential",
    "tslaplace"
  ))) {
    stop("Unknown inefficiency distribution: ", paste(udist),
      call. = FALSE
    )
  }
  # Formula manipulation -------
  if (length(Formula(formula))[2] != 1) {
    stop("argument 'formula' must have one RHS part", call. = FALSE)
  }
  cl <- match.call()
  mc <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset"), names(mc), nomatch = 0L)
  mc <- mc[c(1L, m)]
  mc$drop.unused.levels <- TRUE
  formula <- interCheckMain(formula = formula)
  if (!missing(muhet)) {
    muhet <- lhsCheck_mu(formula = muhet, scaling = scaling)
  } else {
    muhet <- ~1
  }
  if (!missing(uhet)) {
    uhet <- lhsCheck_u(formula = uhet, scaling = scaling)
  } else {
    uhet <- ~1
  }
  if (!missing(vhet)) {
    vhet <- lhsCheck_v(formula = vhet)
  } else {
    vhet <- ~1
  }
  formula <- formDist_sfacross(
    udist = udist, formula = formula,
    muhet = muhet, uhet = uhet, vhet = vhet
  )
  # Generate required datasets -------
  if (missing(data)) {
    data <- environment(formula)
  }
  mc$formula <- formula
  mc$na.action <- na.pass
  mc[[1L]] <- quote(model.frame)
  mc <- eval(mc, parent.frame())
  validObs <- rowSums(is.na(mc) | is.infinite.data.frame(mc)) ==
    0
  Yvar <- model.response(mc, "numeric")
  Yvar <- Yvar[validObs]
  mtX <- terms(formula, data = data, rhs = 1)
  Xvar <- model.matrix(mtX, mc)
  Xvar <- Xvar[validObs, , drop = FALSE]
  nXvar <- ncol(Xvar)
  N <- nrow(Xvar)
  if (N == 0L) {
    stop("0 (non-NA) cases", call. = FALSE)
  }
  if (length(Yvar) != nrow(Xvar)) {
    stop(paste("the number of observations of the dependent variable (",
      length(Yvar), ") must be the same to the number of observations of the exogenous variables (",
      nrow(Xvar), ")",
      sep = ""
    ), call. = FALSE)
  }
  if (udist %in% c("tnormal", "lognormal")) {
    mtmuH <- delete.response(terms(formula,
      data = data,
      rhs = 2
    ))
    muHvar <- model.matrix(mtmuH, mc)
    muHvar <- muHvar[validObs, , drop = FALSE]
    nmuHvar <- ncol(muHvar)
    mtuH <- delete.response(terms(formula, data = data, rhs = 3))
    uHvar <- model.matrix(mtuH, mc)
    uHvar <- uHvar[validObs, , drop = FALSE]
    nuHvar <- ncol(uHvar)
    mtvH <- delete.response(terms(formula, data = data, rhs = 4))
    vHvar <- model.matrix(mtvH, mc)
    vHvar <- vHvar[validObs, , drop = FALSE]
    nvHvar <- ncol(vHvar)
  } else {
    mtuH <- delete.response(terms(formula, data = data, rhs = 2))
    uHvar <- model.matrix(mtuH, mc)
    uHvar <- uHvar[validObs, , drop = FALSE]
    nuHvar <- ncol(uHvar)
    mtvH <- delete.response(terms(formula, data = data, rhs = 3))
    vHvar <- model.matrix(mtvH, mc)
    vHvar <- vHvar[validObs, , drop = FALSE]
    nvHvar <- ncol(vHvar)
  }
  # Check other supplied options -------
  if (length(S) != 1 || !(S %in% c(-1L, 1L))) {
    stop("argument 'S' must equal either 1 or -1: 1 for production or profit frontier
   and -1 for cost frontier",
      call. = FALSE
    )
  }
  typeSfa <- if (S == 1L) {
    "Stochastic Production/Profit Frontier, e = v - u"
  } else {
    "Stochastic Cost Frontier, e = v + u"
  }
  if (length(scaling) != 1 || !is.logical(scaling[1])) {
    stop("argument 'scaling' must be a single logical value",
      call. = FALSE
    )
  }
  if (scaling) {
    if (udist != "tnormal") {
      stop("argument 'udist' must be 'tnormal' when scaling option is TRUE",
        call. = FALSE
      )
    }
    if (nuHvar != nmuHvar) {
      stop("argument 'muhet' and 'uhet' must have the same length",
        call. = FALSE
      )
    }
    if (!all(colnames(uHvar) == colnames(muHvar))) {
      stop("argument 'muhet' and 'uhet' must contain the same variables",
        call. = FALSE
      )
    }
    if (nuHvar == 1 || nmuHvar == 1) {
      if (attr(terms(muhet), "intercept") == 1 || attr(
        terms(uhet),
        "intercept"
      ) == 1) {
        stop("at least one exogeneous variable must be provided for the scaling option",
          call. = FALSE
        )
      }
    }
  }
  if (length(logDepVar) != 1 || !is.logical(logDepVar[1])) {
    stop("argument 'logDepVar' must be a single logical value",
      call. = FALSE
    )
  }
  # Number of parameters -------
  nParm <- if (udist == "tnormal") {
    if (scaling) {
      if (attr(terms(muhet), "intercept") == 1 || attr(
        terms(uhet),
        "intercept"
      ) == 1) {
        nXvar + (nmuHvar - 1) + 2 + nvHvar
      } else {
        nXvar + nmuHvar + 2 + nvHvar
      }
    } else {
      nXvar + nmuHvar + nuHvar + nvHvar
    }
  } else {
    if (udist == "lognormal") {
      nXvar + nmuHvar + nuHvar + nvHvar
    } else {
      if (udist %in% c("gamma", "weibull", "tslaplace")) {
        nXvar + nuHvar + nvHvar + 1
      } else {
        nXvar + nuHvar + nvHvar
      }
    }
  }
  # Checking starting values when provided -------
  if (!is.null(start)) {
    if (length(start) != nParm) {
      stop("Wrong number of initial values: model has ",
        nParm, " parameters",
        call. = FALSE
      )
    }
  }
  if (nParm > N) {
    stop("Model has more parameters than observations", call. = FALSE)
  }
  # Check algorithms -------
  method <- tolower(method)
  if (!(method %in% c(
    "ucminf", "bfgs", "bhhh", "nr", "nm",
    "sr1", "mla", "sparse", "nlminb"
  ))) {
    stop("Unknown or non-available optimization algorithm: ",
      paste(method),
      call. = FALSE
    )
  }
  # Check hessian type
  if (length(hessianType) != 1 || !(hessianType %in% c(
    1L,
    2L, 3L
  ))) {
    stop("argument 'hessianType' must equal either 1 or 2 or 3",
      call. = FALSE
    )
  }
  # Draws for SML -------
  if (udist %in% c("gamma", "lognormal", "weibull")) {
    if (!(simType %in% c("halton", "ghalton", "sobol", "uniform"))) {
      stop("Unknown or non-available random draws method",
        call. = FALSE
      )
    }
    if (!is.numeric(Nsim) || length(Nsim) != 1) {
      stop("argument 'Nsim' must be a single numeric scalar",
        call. = FALSE
      )
    }
    if (!is.numeric(burn) || length(burn) != 1) {
      stop("argument 'burn' must be a single numeric scalar",
        call. = FALSE
      )
    }
    if (!is_prime(prime)) {
      stop("argument 'prime' must be a single prime number",
        call. = FALSE
      )
    }
    if (length(antithetics) != 1 || !is.logical(antithetics[1])) {
      stop("argument 'antithetics' must be a single logical value",
        call. = FALSE
      )
    }
    if (antithetics && (Nsim %% 2) != 0) {
      Nsim <- Nsim + 1
    }
    simDist <- if (simType == "halton") {
      "Halton"
    } else {
      if (simType == "ghalton") {
        "Generalized Halton"
      } else {
        if (simType == "sobol") {
          "Sobol"
        } else {
          if (simType == "uniform") {
            "Uniform"
          }
        }
      }
    }
    cat("Initialization of", Nsim, simDist, "draws per observation ...\n")
    FiMat <- drawMat(
      N = N, Nsim = Nsim, simType = simType,
      prime = prime, burn = burn + 1, antithetics = antithetics,
      seed = seed
    )
  }
  # Other optimization options -------
  if (!is.numeric(itermax) || length(itermax) != 1) {
    stop("argument 'itermax' must be a single numeric scalar",
      call. = FALSE
    )
  }
  if (itermax != round(itermax)) {
    stop("argument 'itermax' must be an integer", call. = FALSE)
  }
  if (itermax <= 0) {
    stop("argument 'itermax' must be positive", call. = FALSE)
  }
  itermax <- as.integer(itermax)
  if (length(printInfo) != 1 || !is.logical(printInfo[1])) {
    stop("argument 'printInfo' must be a single logical value",
      call. = FALSE
    )
  }
  if (!is.numeric(tol) || length(tol) != 1) {
    stop("argument 'tol' must be numeric", call. = FALSE)
  }
  if (tol < 0) {
    stop("argument 'tol' must be non-negative", call. = FALSE)
  }
  if (!is.numeric(gradtol) || length(gradtol) != 1) {
    stop("argument 'gradtol' must be numeric", call. = FALSE)
  }
  if (gradtol < 0) {
    stop("argument 'gradtol' must be non-negative", call. = FALSE)
  }
  if (!is.numeric(stepmax) || length(stepmax) != 1) {
    stop("argument 'stepmax' must be numeric", call. = FALSE)
  }
  if (stepmax < 0) {
    stop("argument 'stepmax' must be non-negative", call. = FALSE)
  }
  if (!(qac %in% c("marquardt", "stephalving"))) {
    stop("argument 'qac' must be either 'marquardt' or 'stephalving'",
      call. = FALSE
    )
  }
  # Step 1: OLS -------
  olsRes <- if (colnames(Xvar)[1] == "(Intercept)") {
    lm(Yvar ~ ., data = as.data.frame(Xvar[, -1]))
  } else {
    lm(Yvar ~ -1 + ., data = as.data.frame(Xvar))
  }
  if (any(is.na(olsRes$coefficients))) {
    stop("at least one of the OLS coefficients is NA: ",
      paste(colnames(Xvar)[is.na(olsRes$coefficients)],
        collapse = ", "
      ), "This may be due to a singular matrix
   due to potential perfect multicollinearity",
      call. = FALSE
    )
  }
  olsParam <- c(olsRes$coefficients)
  olsSigmasq <- summary(olsRes)$sigma^2
  olsStder <- sqrt(diag(vcov(olsRes)))
  olsLoglik <- logLik(olsRes)[1]
  if (inherits(data, "plm.dim")) {
    dataTable <- data[validObs, 1:2]
  } else {
    dataTable <- data.frame(IdObs = c(1:sum(validObs)))
  }
  dataTable <- as_tibble(cbind(dataTable, data[validObs, all.vars(terms(formula))]))
  dataTable <- mutate(dataTable,
    olsResiduals = residuals(olsRes),
    olsFitted = fitted(olsRes)
  )
  olsSkew <- skewness(dataTable[["olsResiduals"]])
  olsM3Okay <- if (S * olsSkew < 0) {
    "Residuals have the 'right' skeweness"
  } else {
    "Residuals have the 'wrong' skeweness"
  }
  if (S * olsSkew > 0) {
    warning("The residuals of the OLS are ", if (S == 1) {
      " right"
    } else {
      "left"
    }, "-skewed. This may indicate the absence of inefficiency or
  model misspecification or sample 'bad luck'",
    call. = FALSE
    )
  }
  CoelliM3Test <- c(z = moment(dataTable[["olsResiduals"]],
    order = 3
  ) / sqrt(6 * moment(dataTable[["olsResiduals"]],
    order = 2
  )^3 / N), p.value = 2 * pnorm(-abs(moment(dataTable[["olsResiduals"]],
    order = 3
  ) / sqrt(6 * moment(dataTable[["olsResiduals"]],
    order = 2
  )^3 / N))))
  AgostinoTest <- dagoTest(dataTable[["olsResiduals"]])@test
  # Step 2: MLE arguments -------
  FunArgs <- if (udist == "tnormal") {
    if (scaling) {
      list(
        start = start, olsParam = olsParam, dataTable = dataTable,
        nXvar = nXvar, nuHvar = nuHvar, nvHvar = nvHvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, method = method, printInfo = printInfo,
        itermax = itermax, stepmax = stepmax, tol = tol,
        gradtol = gradtol, hessianType = hessianType,
        qac = qac
      )
    } else {
      list(
        start = start, olsParam = olsParam, dataTable = dataTable,
        nXvar = nXvar, nmuHvar = nmuHvar, nuHvar = nuHvar,
        nvHvar = nvHvar, muHvar = muHvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        method = method, printInfo = printInfo, itermax = itermax,
        stepmax = stepmax, tol = tol, gradtol = gradtol,
        hessianType = hessianType, qac = qac
      )
    }
  } else {
    if (udist == "lognormal") {
      list(
        start = start, olsParam = olsParam, dataTable = dataTable,
        nXvar = nXvar, nmuHvar = nmuHvar, nuHvar = nuHvar,
        nvHvar = nvHvar, muHvar = muHvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        N = N, FiMat = FiMat, method = method, printInfo = printInfo,
        itermax = itermax, stepmax = stepmax, tol = tol,
        gradtol = gradtol, hessianType = hessianType,
        qac = qac
      )
    } else {
      if (udist %in% c("gamma", "weibull")) {
        list(
          start = start, olsParam = olsParam, dataTable = dataTable,
          nXvar = nXvar, nuHvar = nuHvar, nvHvar = nvHvar,
          uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, S = S, N = N, FiMat = FiMat, method = method,
          printInfo = printInfo, itermax = itermax, stepmax = stepmax,
          tol = tol, gradtol = gradtol, hessianType = hessianType,
          qac = qac
        )
      } else {
        list(
          start = start, olsParam = olsParam, dataTable = dataTable,
          nXvar = nXvar, nuHvar = nuHvar, nvHvar = nvHvar,
          uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, S = S, method = method, printInfo = printInfo,
          itermax = itermax, stepmax = stepmax, tol = tol,
          gradtol = gradtol, hessianType = hessianType,
          qac = qac
        )
      }
    }
  }
  ## MLE run -------
  mleList <- tryCatch(switch(udist, hnormal = do.call(
    halfnormAlgOpt,
    FunArgs
  ), exponential = do.call(exponormAlgOpt, FunArgs),
  tnormal = if (scaling) {
    do.call(truncnormscalAlgOpt, FunArgs)
  } else {
    do.call(
      truncnormAlgOpt,
      FunArgs
    )
  }, rayleigh = do.call(raynormAlgOpt, FunArgs),
  gamma = do.call(gammanormAlgOpt, FunArgs), uniform = do.call(
    uninormAlgOpt,
    FunArgs
  ), lognormal = do.call(lognormAlgOpt, FunArgs),
  weibull = do.call(weibullnormAlgOpt, FunArgs), genexponential = do.call(
    genexponormAlgOpt,
    FunArgs
  ), tslaplace = do.call(tslnormAlgOpt, FunArgs)
  ),
  error = function(e) e
  )
  if (inherits(mleList, "error")) {
    stop("The current error occurs during optimization:\n",
      mleList$message,
      call. = FALSE
    )
  }
  # Inverse Hessian + other -------
  mleList$invHessian <- vcovObj(
    mleObj = mleList$mleObj, hessianType = hessianType,
    method = method, nParm = nParm
  )
  mleList <- c(mleList, if (method == "ucminf") {
    list(
      type = "ucminf max.", nIter = unname(mleList$mleObj$info["neval"]),
      status = mleList$mleObj$message, mleLoglik = -mleList$mleObj$value,
      gradient = mleList$mleObj$gradient
    )
  } else {
    if (method %in% c("bfgs", "bhhh", "nr", "nm")) {
      list(
        type = substr(mleList$mleObj$type, 1, 27), nIter = mleList$mleObj$iterations,
        status = mleList$mleObj$message, mleLoglik = mleList$mleObj$maximum,
        gradient = mleList$mleObj$gradient
      )
    } else {
      if (method == "sr1") {
        list(
          type = "SR1 max.", nIter = mleList$mleObj$iterations,
          status = mleList$mleObj$status, mleLoglik = -mleList$mleObj$fval,
          gradient = mleList$mleObj$gradient
        )
      } else {
        if (method == "mla") {
          list(
            type = "Lev. Marquardt max.", nIter = mleList$mleObj$ni,
            status = switch(mleList$mleObj$istop, `1` = "convergence criteria were satisfied",
              `2` = "maximum number of iterations was reached",
              `4` = "algorithm encountered a problem in the function computation"
            ),
            mleLoglik = -mleList$mleObj$fn.value, gradient = mleList$mleObj$grad
          )
        } else {
          if (method == "sparse") {
            list(
              type = "Sparse Hessian max.", nIter = mleList$mleObj$iterations,
              status = mleList$mleObj$status, mleLoglik = -mleList$mleObj$fval,
              gradient = mleList$mleObj$gradient
            )
          } else {
            if (method == "nlminb") {
              list(
                type = "nlminb max.", nIter = mleList$mleObj$iterations,
                status = mleList$mleObj$message, mleLoglik = -mleList$mleObj$objective,
                gradient = mleList$mleObj$gradient
              )
            }
          }
        }
      }
    }
  })
  # quick renaming -------
  if (udist %in% c("tnormal", "lognormal")) {
    names(mleList$startVal) <- fName_mu_sfacross(
      Xvar = Xvar,
      udist = udist, muHvar = muHvar, uHvar = uHvar, vHvar = vHvar,
      scaling = scaling
    )
  } else {
    names(mleList$startVal) <- fName_uv_sfacross(
      Xvar = Xvar,
      udist = udist, uHvar = uHvar, vHvar = vHvar
    )
  }
  names(mleList$mleParam) <- names(mleList$startVal)
  rownames(mleList$invHessian) <- colnames(mleList$invHessian) <- names(mleList$mleParam)
  names(mleList$gradient) <- names(mleList$mleParam)
  colnames(mleList$mleObj$gradL_OBS) <- names(mleList$mleParam)
  # Return object -------
  mleDate <- format(Sys.time(), "Model was estimated on : %b %a %d, %Y at %H:%M")
  dataTable$mleResiduals <- Yvar - as.numeric(crossprod(
    matrix(mleList$mleParam[1:nXvar]),
    t(Xvar)
  ))
  dataTable$mleFitted <- as.numeric(crossprod(
    matrix(mleList$mleParam[1:nXvar]),
    t(Xvar)
  ))
  dataTable$logL_OBS <- mleList$mleObj$logL_OBS
  returnObj <- list()
  returnObj$call <- cl
  returnObj$formula <- formula
  returnObj$S <- S
  returnObj$typeSfa <- typeSfa
  returnObj$Nobs <- N
  returnObj$nXvar <- nXvar
  if (udist %in% c("tnormal", "lognormal")) {
    returnObj$nmuHvar <- nmuHvar
  }
  returnObj$scaling <- scaling
  returnObj$logDepVar <- logDepVar
  returnObj$nuHvar <- nuHvar
  returnObj$nvHvar <- nvHvar
  returnObj$nParm <- nParm
  returnObj$udist <- udist
  returnObj$startVal <- mleList$startVal
  returnObj$dataTable <- dataTable
  returnObj$olsParam <- olsParam
  returnObj$olsStder <- olsStder
  returnObj$olsSigmasq <- olsSigmasq
  returnObj$olsLoglik <- olsLoglik
  returnObj$olsSkew <- olsSkew
  returnObj$olsM3Okay <- olsM3Okay
  returnObj$CoelliM3Test <- CoelliM3Test
  returnObj$AgostinoTest <- AgostinoTest
  returnObj$optType <- mleList$type
  returnObj$nIter <- mleList$nIter
  returnObj$optStatus <- mleList$status
  returnObj$startLoglik <- mleList$startLoglik
  returnObj$mleLoglik <- mleList$mleLoglik
  returnObj$mleParam <- mleList$mleParam
  returnObj$gradient <- mleList$gradient
  returnObj$gradL_OBS <- mleList$mleObj$gradL_OBS
  returnObj$gradientNorm <- sqrt(sum(mleList$gradient^2))
  returnObj$invHessian <- mleList$invHessian
  returnObj$hessianType <- if (hessianType == 1) {
    "Analytic/Numeric Hessian"
  } else {
    if (hessianType == 2) {
      "BHHH Hessian"
    } else {
      if (hessianType == 3) {
        "Robust Hessian"
      }
    }
  }
  returnObj$mleDate <- mleDate
  # returnObj$validObs <- validObs
  if (udist %in% c("gamma", "lognormal", "weibull")) {
    returnObj$simDist <- simDist
    returnObj$Nsim <- Nsim
    returnObj$FiMat <- FiMat
  }
  rm(mleList)
  class(returnObj) <- "sfacross"
  print.sfacross(returnObj)
  return(returnObj)
}

# print for sfacross ----------
#' @rdname sfacross
#' @aliases print.sfacross
#' @export
print.sfacross <- function(x, digits = max(3, getOption("digits") - 2),
                           ...) {
  cat("\nCall:\n")
  cat(deparse(x$call))
  cat("\n\n")
  cat("Likelihood estimates using", x$optType, "\n")
  cat(sfadist(x$udist), "\n\n")
  cat(x$typeSfa, "\n\n")
  print.default(format(x$mleParam, digits = digits),
    print.gap = 2,
    quote = FALSE
  )
  invisible(x)
}
