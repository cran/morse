% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LCx.R, R/LCx.survFit.R
\name{LCx}
\alias{LCx}
\alias{LCx.survFit}
\title{Predict \eqn{X}\% Lethal Concentration at the maximum time point (default).}
\usage{
LCx(object, ...)

\method{LCx}{survFit}(object, X, time_LCx = NULL, conc_range = NULL, npoints = 100, ...)
}
\arguments{
\item{object}{An object of class \code{survFit}}

\item{\dots}{Further arguments to be passed to generic methods}

\item{X}{Percentage of individuals dying (e.g., \eqn{50} for \eqn{LC_{50}}, \eqn{10} for \eqn{LC_{10}}, ...)}

\item{time_LCx}{A number giving the time at which  \eqn{LC_{x}} has to be estimated. 
If NULL, the latest time point of the experiment is used.}

\item{conc_range}{A vector of length 2 with minimal and maximal value of the 
range of concentration. If NULL, the range is
define between 0 and the highest tested concentration of the experiment.}

\item{npoints}{Number of time point in \code{conc_range} between 0 and the maximal concentration. 100 by default.}
}
\value{
returns an object of class \code{LCx}.

The function returns an object of class \code{LCx}, which is a list
 with the following information:
\item{X_prop}{Survival probability of individuals surviving considering the median
 of the background mortality (i.e. \eqn{S(0, t)*(1- x/100)})}
\item{X_prop_provided}{Survival probability of individuals surviving as provided in arguments (i.e. \eqn{(100-X)/100)}}
\item{time_LCx}{A number giving the time at which  \eqn{LC_{x}} has to be
 estimated as provided in arguments or if NULL, the latest time point of the
  experiment is used.}
\item{df_LCx}{A \code{data.frame} with quantiles (median, 2.5\% and 97.5\%)
 of \eqn{LC_{X}} at time \code{time_LCx} for \eqn{X}\% of individuals}
\item{df_dose}{A \code{data.frame} with four columns: \code{concentration}, and median \code{q50} and 95\% credible interval
 (\code{qinf95} and \code{qsup95}) of the survival probability at time \code{time_LCx}}
}
\description{
Predict median and 95\% credible interval of the x\% Lethal Concentration.

The function \code{LCx}, \eqn{x}\% Lethal Concentration (\eqn{LC_x}), is use to compute
 the dose required to kill \eqn{x}\% of the members of a tested population
 after a specified test duration (\code{time_LCx}) (default is the maximum
 time point of the experiment).
 
 Mathematical definition of \eqn{x}\% Lethal Concentration at time \eqn{t},
 denoted \eqn{LC(x,t)}, is:
 
 \eqn{S(LC(x,t), t) = S(0, t)*(1- x/100)},
 
 where \eqn{S(LC(x,t), t)} is the survival probability at concentration
 \eqn{LC(x,t)} at time \eqn{t}, and \eqn{S(0,t)} is the survival probability at
 no concentration (i.e. concentration is \eqn{0}) at time \eqn{t} which
 reflect the background mortality \eqn{h_b}:
 
 \eqn{S(0, t) = exp(-hb* t)}.
  
 In the function \code{LCx}, we use the median of \eqn{S(0,t)} to rescale the
 \eqn{x}\% Lethal Concentration at time \eqn{t}.
}
\details{
When class of \code{object} is \code{survFit}, see \link[=LCx.survFit]{LCx.survFit}.
}
