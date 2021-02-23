\docType{package}

\name{sfaR-package}
\alias{sfaR-package}
\alias{sfaR}

\title{sfaR: An R package for estimating stochastic frontier models}

\description{The sfaR package provides a set of tools (maximum likelihood -
MLE and maximum simulated likelihood - MSL) for various specifications of
stochastic frontier analysis (SFA).

Two categories of important functions are available:
\code{sfacross} and \code{lcmcross}, which estimate different
type of frontiers and offer nine alternative optimization algorithms
(i.e. "bfgs", "bhhh", "nr", "nm", "ucminf", "mla", "sr1", "sparse" and
"nlminb").}

\section{sfacross}{
\code{\link{sfacross}} estimates the frontier for cross/pooled sectional data
and allows for ten different distributions for the one-sided error term.
These distributions comprise the exponential, the gamma, the generalized
exponential, the half normal, the lognormal, the truncated normal, the
truncated skewed Laplace, the rayleigh, the uniform and the weibull distributions. 
In the case of the gamma, lognormal and weibull distributions, maximum simulated
likelihood is used with the possibility of four specific distributions to
construct the draws: halton, generalized halton, sobol and uniform.
Heteroscedasticity in both error terms can be accommodated, in addition to
heterogeneity in the truncated mean parameter in the case of the truncated
normal and lognormal distributions. In the case of the truncated normal
distribution, the scaling property can also be estimated.
}

\section{lcmcross}{
\code{\link{lcmcross}} estimates latent class stochastic frontier models (LCM), 
which accounts for technological heterogeneity by splitting the observations into a
maximum number of five classes. The classification operates based on a logit
functional form that can be specified using some covariates (namely, the separating 
variables allowing the separation of observations in several classes). Only the half
normal distribution is available for the one-sided error term. Moreover,
heteroscedasticity in both error terms is possible. The choice of the number
of classes can be guided by several information criteria (AIC, BIC or HQIC).
}

\author{K Hervé Dakpo, Yann Desjeux and Laure Latruffe}

\section{bugreport}{Any bug or suggestion can be reported using the \code{sfaR}'\code{s} 
tracker facilities at \url{https://r-forge.r-project.org/projects/sfar/}
}
