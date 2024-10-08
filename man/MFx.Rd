% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MFx.R, R/MFx.survFit.R, R/MFx_ode.survFit.R
\name{MFx}
\alias{MFx}
\alias{MFx.survFit}
\alias{MFx_ode}
\alias{MFx_ode.survFit}
\title{Predict the Multiplication Factor leading to x\% of reduction in survival
at a specific time.}
\usage{
MFx(object, ...)

\method{MFx}{survFit}(
  object,
  data_predict,
  X = 50,
  time_MFx = NULL,
  MFx_range = c(0, 1000),
  mcmc_size = 1000,
  hb_value = TRUE,
  spaghetti = FALSE,
  accuracy = 0.01,
  quiet = FALSE,
  threshold_iter = 100,
  hb_valueFORCED = 0,
  ode = TRUE,
  interpolate_length = NULL,
  interpolate_method = "linear",
  ...
)

MFx_ode(object, ...)

\method{MFx_ode}{survFit}(
  object,
  data_predict,
  X = 50,
  time_MFx = NULL,
  MFx_range = c(0, 1000),
  mcmc_size = 1000,
  hb_value = TRUE,
  spaghetti = FALSE,
  accuracy = 0.01,
  quiet = FALSE,
  threshold_iter = 100,
  hb_valueFORCED = 0,
  interpolate_length = NULL,
  interpolate_method = "linear",
  ...
)
}
\arguments{
\item{object}{An object of class \code{survFit}.}

\item{\dots}{Further arguments to be passed to generic methods}

\item{data_predict}{A dataframe with two columns \code{time} and \code{conc}.}

\item{X}{Percentage of survival change (e.g., \eqn{50} for survival decrease of 50\%
, or \eqn{-50} for survival increase of 50\%).The default is 50. 
Only time series computed during the adaptation using a binary search in
\eqn{O(log(n))} are returned. However, if \code{NULL}, all time series
computed from the vector \code{MFx_range} are returned.}

\item{time_MFx}{A number giving the time at which  \eqn{MF(x,t)} has to be estimated. 
If NULL, the latest time point of the profile is used.}

\item{MFx_range}{A vector from which lower and upper bound of the range of the
 multiplication factor \code{MFx} are generated. The default is a vector \code{c(0, 1000)}.
If argument \code{X} is \code{NULL}, then all the time series generated with
\code{MFx_range} are returned.}

\item{mcmc_size}{Can be used to reduce the number of MCMC samples in order to speed up
the computation. The default is 1000.}

\item{hb_value}{If \code{TRUE}, the background mortality \code{hb} is taken into
 account from the posterior.
If \code{FALSE}, parameter \code{hb} is set to 0. The default is \code{TRUE}.}

\item{spaghetti}{If \code{TRUE}, return a set of survival curves using
parameters drawn from the posterior distribution.}

\item{accuracy}{Accuracy of the multiplication factor. The default is 0.01.}

\item{quiet}{If \code{FALSE}, print the evolution of accuracy.}

\item{threshold_iter}{Threshold number of iteration.}

\item{hb_valueFORCED}{If \code{hb_value} is \code{FALSE}, it fix \code{hb}.}

\item{ode}{IF \code{ode} is \code{TRUE}, algo use predict_ode rather than predict. Default is \code{TRUE}.}

\item{interpolate_length}{Length of the time sequence for which output is wanted.}

\item{interpolate_method}{The interpolation method for concentration. See package \code{deSolve} for details.
Default is \code{linear}.}
}
\value{
returns an object of class \code{MFx}

The function returns an object of class \code{MFx}, which is a list
 with the following information:
 \item{X_prop}{Survival probability for \code{X} percent of reduction of the initial median 
survival probability at time \code{time_MFx}.}
\item{X_prop_provided}{A number giving the proportion of reduction in survival.}
\item{time_MFx}{A number giving the time at which  \eqn{MF(x,t)} has to be
 estimated as provided in arguments or if NULL, the latest time point of the
  profile is used.}
\item{df_MFx}{A \code{data.frame} with quantiles (median, 2.5\% and 97.5\%)
 of \eqn{MF(x,t)} at time \eqn{t}, \code{time_MFx}, for \eqn{x}\% of survival reduction.}
\item{df_dose}{A \code{data.frame} with quantiles (median, 2.5\% and 97.5\%)
 of survival probability along the computed multiplication factor and at time \code{time_MFx}.}
\item{MFx_tested}{A vector of all multiplication factors computed.} 
\item{ls_predict}{A list of all object of class \code{survFitPredict} obtained
from computing survival probability for every profiles build from the vector of
multiplication factors \code{MFx_tested}.}

The function returns an object of class \code{MFx}, which is a list
 with the following information:
 \item{X_prop}{Survival probability for \code{X} percent of reduction of the initial median 
survival probability at time \code{time_MFx}.}
\item{X_prop_provided}{A number giving the proportion of reduction in survival.}
\item{time_MFx}{A number giving the time at which  \eqn{MF(x,t)} has to be
 estimated as provided in arguments or if NULL, the latest time point of the
  profile is used.}
\item{df_MFx}{A \code{data.frame} with quantiles (median, 2.5\% and 97.5\%)
 of \eqn{MF(x,t)} at time \eqn{t}, \code{time_MFx}, for \eqn{x}\% of survival reduction.}
\item{df_dose}{A \code{data.frame} with quantiles (median, 2.5\% and 97.5\%)
 of survival probability along the computed multiplication factor and at time \code{time_MFx}.}
\item{MFx_tested}{A vector of all multiplication factors computed.} 
\item{ls_predict}{A list of all object of class \code{survFitPredict} obtained
from computing survival probability for every profiles build from the vector of
multiplication factors \code{MFx_tested}.}
}
\description{
Generic method for \code{MFx}, a function denoted \eqn{MF(x,t)} for 
\eqn{x}\% Multiplication Factor at time \eqn{t}.

The function \code{MFx}, \eqn{x}\% Multiplication Factor at time \eqn{t}, (\eqn{MF(x,t)}),
is used to compute the multiplication factor
applied to the concentration exposure profile in order to
reduce by \eqn{x}\% (argument \code{X}) the survival probability at a
 specified test duration \eqn{t} (argument \code{time_MFx}) (default is the maximum
 time point of the experiment).
 
 Mathematical definition of \eqn{x}\% Multiplication Factor at time \eqn{t}
 (at the end of a time series \eqn{T = \{0, \dots, t\}}),
 denoted \eqn{MF(x,t)}, is given by:
 
 \eqn{S(MF(x,t) * C_w(\tau \in T), t) = S( C_w(\tau \in T), t)*(1- x/100)},
 
 where \eqn{C_w(\tau \in T)} is the initial exposure profile without
 multiplication factor. And so the expression \eqn{S(MF(x,t)* C_w(\tau \in T), t)}
  is the survival probability after an exposure profile
   \eqn{MF(x,t)* C_w(\tau \in T)} at time \eqn{t}.

This is a \code{method} to replace function \code{MFx} used on \code{survFit}
object when computing issues happen. \code{MFx_ode} uses the \code{deSolve}
library to improve robustness. However, time to compute may be longer.

The function \code{MFx_ode}, \eqn{x}\% Multiplication Factor at time \eqn{t}, (\eqn{MF(x,t)}),
is used to compute the multiplication factor
applied to the concentration exposure profile in order to
reduce by \eqn{x}\% (argument \code{X}) the survival probability at a
 specified test duration \eqn{t} (argument \code{time_MFx}) (default is the maximum
 time point of the experiment).
 
 Mathematical definition of \eqn{x}\% Multiplication Factor at time \eqn{t}
 (at the end of a time series \eqn{T = \{0, \dots, t\}}),
 denoted \eqn{MF(x,t)}, is given by:
 
 \eqn{S(MF(x,t) * C_w(\tau \in T), t) = S( C_w(\tau \in T), t)*(1- x/100)},
 
 where \eqn{C_w(\tau \in T)} is the initial exposure profile without
 multiplication factor. And so the expression \eqn{S(MF(x,t)* C_w(\tau \in T), t)}
  is the survival probability after an exposure profile
   \eqn{MF(x,t)* C_w(\tau \in T)} at time \eqn{t}.
}
\details{
When class of \code{object} is \code{survFit}, see \link[=MFx.survFit]{MFx.survFit}.
}
