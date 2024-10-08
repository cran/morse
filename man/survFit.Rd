% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survFit.R, R/survFit.survDataCstExp.R,
%   R/survFit.survDataVarExp.R
\name{survFit}
\alias{survFit}
\alias{survFit.survDataCstExp}
\alias{survFit.survDataVarExp}
\title{Fits a TKTD model for survival analysis using Bayesian inference}
\usage{
survFit(
  data,
  model_type,
  quiet,
  n.chains,
  n.adapt,
  n.iter,
  n.warmup,
  thin.interval,
  limit.sampling,
  dic.compute,
  dic.type,
  hb_value,
  hb_valueFIXED,
  ...
)

\method{survFit}{survDataCstExp}(
  data,
  model_type = NULL,
  quiet = FALSE,
  n.chains = 3,
  n.adapt = 3000,
  n.iter = NULL,
  n.warmup = NULL,
  thin.interval = NULL,
  limit.sampling = TRUE,
  dic.compute = FALSE,
  dic.type = "pD",
  hb_value = TRUE,
  hb_valueFIXED = NA,
  ...
)

\method{survFit}{survDataVarExp}(
  data,
  model_type = NULL,
  quiet = FALSE,
  n.chains = 3,
  n.adapt = 1000,
  n.iter = NULL,
  n.warmup = NULL,
  thin.interval = NULL,
  limit.sampling = TRUE,
  dic.compute = FALSE,
  dic.type = "pD",
  hb_value = TRUE,
  hb_valueFIXED = NA,
  extend_time = 100,
  ...
)
}
\arguments{
\item{data}{An object of class \code{survDataCstExp} or \code{survDataVarExp}.}

\item{model_type}{Can be \code{"SD"} or \code{"IT"} to choose
between "Stochastic Death" or "Individual Tolerance" models
(resp.). See the modeling vignette for details.}

\item{quiet}{If \code{FALSE}, prints logs and progress bar from
JAGS.}

\item{n.chains}{A positive integer specifying the number of MCMC chains. The minimum required number 
of chains is 2.}

\item{n.adapt}{A positive integer specifying the number of iterations for adaptation. If \code{n.adapt} = 0
then no adaptation takes place.}

\item{n.iter}{A positive integer specifying the number of iterations to monitor for each chain.}

\item{n.warmup}{A positive integer specifying the number of warmup (aka burnin) iterations per chain.}

\item{thin.interval}{A positive integer specifying the period to monitor.}

\item{limit.sampling}{if \code{FALSE} (default is \code{TRUE}), there is no limit to the number of iterations
in MCMC imposed by the \code{raftery.diag} test.}

\item{dic.compute}{if \code{TRUE} (default is \code{FALSE}), it generates penalized deviance samples to compute
the Deviance Information Criterion (DIC) with the \code{rjags} package}

\item{dic.type}{type of penalty to use. A string identifying the type of penalty: \code{pD} or \code{popt}
(see function \code{\link[rjags]{dic.samples}})}

\item{hb_value}{If \code{TRUE}, the background mortality \code{hb} is taken into account.
If \code{FALSE}, parameter \code{hb} is set to 0. The default is \code{TRUE}.}

\item{hb_valueFIXED}{If \code{hb_value} is \code{FALSE}, then \code{hb_valueFiXED} is the value to fix \code{hb}.
If \code{hb_value} is \code{FALSE} and  \code{hb_valueFiXED} is \code{NA}, then \code{hb} is fixed to \code{0}.}

\item{\dots}{Further arguments to be passed to generic methods}

\item{extend_time}{Number of for each replicate used for linear 
interpolation (comprise between time to compute and fitting accuracy)}
}
\value{
an object of class \code{survFit}

The function returns an object of class \code{survFitCstExp}, which is
a list with the following information:
\item{estim.par}{a table of the estimated parameters as medians and 95\%
credible intervals}
\item{mcmc}{an object of class \code{mcmc.list} with the posterior
distribution}
\item{model}{a JAGS model object}
\item{dic}{return the Deviance Information Criterion (DIC) if \code{dic.compute} is \code{TRUE}}
\item{warnings}{a table with warning messages}
\item{parameters}{a list of parameter names used in the model}
\item{n.chains}{an integer value corresponding to the number of chains used
for the MCMC computation}
\item{mcmcInfo}{a table with the number of iterations, chains, adaptation, warmup and the thinning interval.} 
\item{jags.data}{a list of the data passed to the JAGS model}
\item{model_type}{the type of TKTD model used: \code{SD} or \code{IT}}

The function returns an object of class \code{survFitVarExp}, which is
a list with the following information:
\item{estim.par}{a table of the estimated parameters as medians and 95\%
credible intervals}
\item{mcmc}{an object of class \code{mcmc.list} with the posterior
distribution}
\item{model}{a JAGS model object}
\item{dic}{return the Deviance Information Criterion (DIC) if \code{dic.compute} is \code{TRUE}}
\item{warnings}{a table with warning messages}
\item{parameters}{a list of parameter names used in the model}
\item{n.chains}{an integer value corresponding to the number of chains used
for the MCMC computation}
\item{mcmcInfo}{a table with the number of iterations, chains, adaptation, warmup and the thinning interval} 
\item{jags.data}{a list of the data passed to the JAGS model}
\item{model_type}{the type of TKTD model used: \code{SD} or \code{IT}}
}
\description{
This function estimates the parameters of a TKTD model ('SD' or 'IT')
for survival analysis using Bayesian inference. In this model,
the survival rate of individuals is modeled as a function of the chemical compound
concentration with a mechanistic description of the effects on survival over
time.
}
\details{
The function \code{survFit} returns the parameter estimates of Toxicokinetic-toxicodynamic (TKTD) models
\code{SD} for 'Stochastic Death' or \code{IT} fo 'Individual Tolerance'.
TKTD models, and particularly the General Unified Threshold model of
Survival (GUTS), provide a consistent process-based
framework to analyse both time and concentration dependent datasets.
In GUTS-SD, all organisms are assumed to have the same internal concentration 
threshold (denoted \eqn{z}), and, once exceeded, the instantaneous probability
to die increases linearly with the internal concentration.
In GUTS-IT, the threshold concentration is distributed among all the organisms, and once 
exceeded in one individual, this individual dies immediately.


When class of \code{object} is \code{survDataCstExp}, see \link[=survFit.survDataCstExp]{survFit.survDataCstExp} ;
and for a \code{survDataVarExp}, see \link[=survFit.survDataVarExp]{survFit.survDataVarExp}.
}
\references{
Jager, T., Albert, C., Preuss, T. G. and Ashauer, R. (2011) 
General unified threshold model of survival-a toxicokinetic-toxicodynamic
 framework for ecotoxicology, \emph{Environmental Science and Technology}, 45, 2529-2540.
303-314.
}
\keyword{estimation}
