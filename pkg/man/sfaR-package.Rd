\docType{package}

\name{sfaR-package}
\alias{sfaR-package}
\alias{sfaR}

\title{sfaR: An R package for estimating stochastic frontier models}

\description{The sfaR package provides a set of tools (maximum likelihood -
MLE and maximum simulated likelihood - MSL) for various specifications of
stochastic frontier analysis. Two categories of important functions:
\code{sfacross} and \code{lcmcross} which estimate different
type of frontiers and offering nine alternative optimization algorithms
("bfgs", "bhhh", "nr", "nm", "ucminf", "mla", "sr1", "sparse" and
"nlminb").}

\section{sfacross}{

\code{\link{sfacross}} estimates the frontier for cross/pooled sectional data
and allows for ten different distributions for the one-sided error term.
These distributions comprise the exponential, the gamma, the generalized
exponential, the half normal, the lognormal, the truncated normal, the
truncated skewed laplace, the rayleigh, the uniform and the weibull. In the
case of the gamma, lognormal and weibull distributions maximum simulated
likelihood is used with the possibility of four specific distributions to
construct the draws: halton, Generalized Halton, Sobol and Uniform.
Heteroscedasticity in both error terms can be accommodated in addition to
heterogeneity in the truncated mean parameter in the case of the truncated
normal and lognormal distributions. In the case of the truncated normal
distribution, the scaling property can also be estimated.
}

\section{lcmcross}{

\code{\link{lcmcross}} estimates latent class frontier models, which accounts
for technological heterogeneity by splitting the observations into a
maximum number of five classes. The classification operates based on a logit
functional form that can be specified using some covariates. Only the half
normal distribution is available for the one-sided error term. Moreover,
heteroscedasticity in both error terms is offerered. The choice of the number
of classes can be guided by several information criteria (AIC, BIC or HQIC).
}

\examples{
#######
## NB YD: A mon avis, pas besoin d'exemples ici. Plutôt dans les fonctions concernées.
#######

## Using data on eighty-two countries production (DGP)
# Cobb Douglas (production function) with half normal distribution and
# heteroscedasticity
cb_w_h <- sfacross(formula = ly ~ lk + ll + yr, udist = 'hnormal', uhet = ~
initStat, data = worldprod, S = 1, method = 'bfgs')
  summary(cb_w_h)

# Cobb Douglas (production function) with gamma distribution
cb_w_g <- sfacross(formula = ly ~ lk + ll + yr, udist = 'gamma', uhet = ~ 1,
data = worldprod, S = 1, method = 'bfgs', simType = 'halton', Nsim = 150,
hessianType = 2)
  summary(cb_w_g)

## Data on U.S. electric utility industry.
# Cost frontier with uniform distribution
tl_u_u <- sfacross(formula = log(cost/fprice) ~ log(output) +
I(log(output)^2) + I(log(lprice/fprice)) + I(log(cprice/fprice)), udist =
'uniform', uhet = ~ 1, data = electricity, S = -1, method = 'ucminf')
  summary(tl_u_u)

# Cost frontier with gamma distribution
tl_u_g <- sfacross(formula = log(cost/fprice) ~ log(output) +
I(log(output)^2) + I(log(lprice/fprice)) + I(log(cprice/fprice)), udist =
'gamma', uhet = ~ 1, data = electricity, S = -1, method = 'bfgs', simType =
'halton', Nsim = 200, hessianType = 2)
  summary(tl_u_g)

## Data on Swiss railway.
# Cost frontier with uniform distribution
sr_u_u <- sfacross(formula = LNCT ~ LNQ2 + LNQ3 + LNNET + LNPK + LNPL, udist
= 'uniform', uhet = ~ 1, data = swissrailways, S = -1, method = 'ucminf')
  summary(sr_u_u)
}
  
\references{
Aigner, D., Lovell, C. A. K., & Schmidt, P. (1977). Formulation and
estimation of stochastic frontier production function models. \emph{Journal
of econometrics}, \bold{6}(1):21--37.

Battese, G.E., and T.J. Coelli. 1988. Prediction of Firm-Level Technical
Efficiencies with a Generalized Frontier Production Function and Panel Data.
\emph{Journal of econometrics}, \bold{38}:387--399.

Battese, G. E., & Coelli, T. J. (1995). A model for technical inefficiency
effects in a stochastic frontier production function for panel data.
\emph{Empirical Economics}, \bold{20}(2):325--332.

Bera, A.K., and S.C. Sharma. 1999. Estimating production uncertainty in
stochastic frontier production function models. \emph{Journal of Productivity
Analysis}, \bold{12}:187-- 210.

Caudill, S. B., & Ford, J. M. (1993). Biases in Frontier Estimation Due to
Heteroscedasticity. \emph{Economics Letters}, \bold{41}(1):17--20.

Caudill, S. B., Ford, J. M., & Gropper, D. M. (1995). Frontier Estimation and
Firm-Specific Inefficiency Measures in the Presence of Heteroscedasticity.
\emph{Journal of Business & Economic Statistics}, \bold{13}(1):105--111.

Coelli, T. 1995. Estimators and Hypothesis Tests for a Stochastic Frontier
Function - a Monte-Carlo Analysis. \emph{Journal of Productivity Analysis},
\bold{6}:247--268.

D'Agostino, R., and E.S. Pearson. 1973. Tests for departure from normality.
Empirical results for the distributions of \eqn{b_2} and \eqn{\sqrt{b_1}}.
\emph{Biometrika}, \bold{60}:613--622.

Greene, W. H. (2003). Simulated likelihood estimation of the normal-gamma
stochastic frontier function. \emph{Journal of Productivity Analysis},
\bold{19}(2-3):179--190.

Hadri, K. (1999). Estimation of a doubly heteroscedastic stochastic frontier
cost function. \emph{Journal of Business & Economic Statistics},
\bold{17}(3):359--363.

Hajargasht, G. (2015). Stochastic frontiers with a Rayleigh distribution.
\emph{Journal of Productivity Analysis}, \bold{44}(2):199--208.

Huang, C. J., & Liu, J.-T. (1994). Estimation of a non-neutral stochastic
frontier production function. \emph{Journal of Productivity Analysis},
\bold{5}(2):171--180.

Hjalmarsson, L., S.C. Kumbhakar, and A. Heshmati. 1996. DEA, DFA and SFA: A
comparison. \emph{Journal of Productivity Analysis}, \bold{7}:303--327.

Horrace, W.C., and P. Schmidt. 1996. Confidence statements for efficiency
estimates from stochastic frontier models. \emph{Journal of Productivity
Analysis}, \bold{7}:257--282.

Jondrow, J., C.A.K. Lovell, I.S. Materov, and P. Schmidt. 1982. On the
estimation of technical inefficiency in the stochastic frontier production
function model. \emph{Journal of econometrics}, \bold{19}:233--238.

Kumbhakar, S. C., Ghosh, S., & McGuckin, J. T. (1991). A Generalized
Production Frontier Approach for Estimating Determinants of Inefficiency in
U.S. Dairy Farms. \emph{Journal of Business & Economic Statistics},
\bold{9}(3):279--286.

Li, Q. (1996). Estimating a stochastic production frontier when the adjusted
error is symmetric. \emph{Economics Letters}, \bold{52}(3):221--228.

Meeusen, W., & Vandenbroeck, J. (1977). Efficiency Estimation from
Cobb-Douglas Production Functions with Composed Error. \emph{International
Economic Review}, \bold{18}(2):435--445.

Migon, H. S., & Medici, E. V. (2001). Bayesian hierarchical models for
stochastic production frontier.

Nguyen, N. B. (2010). Estimation of technical efficiency in stochastic
frontier analysis. Bowling Green State University.

Papadopoulos, A. (2020). Stochastic frontier models using the generalized
exponential distribution.

Reifschneider, D., & Stevenson, R. (1991). Systematic Departures from the
Frontier: A Framework for the Analysis of Firm Inefficiency.
\emph{International Economic Review}, \bold{32}(3):715--723.

Stevenson, R. E. (1980). Likelihood Functions for Generalized Stochastic
Frontier Estimation. \emph{Journal of econometrics}, \bold{13}(1):57--66.

Tsionas, E. G. (2007). Efficiency measurement with the weibull stochastic
frontier. \emph{Oxford Bulletin of Economics and Statistics}, \bold{69}(5):
693--706.

Wang, K., & Ye, X. (2020). Development of alternative stochastic frontier
models for estimating time-space prism vertices. \emph{Transportation}.

Wang, H.J., & Schmidt, P. (2002). One-step and two-step estimation of the
effects of exogenous variables on technical efficiency levels. \emph{Journal
of Productivity Analysis}, \bold{18}:129--144.

Wang, J. (2012). A Normal Truncated Skewed-Laplace Model in Stochastic
Frontier Analysis. Western Kentucky University.
}

\author{
K Hervé Dakpo, \href{mailto:k-herve.dakpo@inrae.fr}{k-herve.dakpo@inrae.fr}

Yann Desjeux, \href{mailto:yann.dejeux@inrae.fr}{yann.dejeux@inrae.fr}

Laure Latruffe, \href{mailto:laure.latruffe@inrae.fr}{laure.latruffe@inrae.fr}
}
