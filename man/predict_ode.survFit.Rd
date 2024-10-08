% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_ode.survFit.R
\name{predict_ode.survFit}
\alias{predict_ode.survFit}
\title{Predict method for \code{survFit} objects}
\usage{
\method{predict_ode}{survFit}(
  object,
  data_predict = NULL,
  spaghetti = FALSE,
  mcmc_size = 1000,
  hb_value = TRUE,
  interpolate_length = 100,
  interpolate_method = "linear",
  hb_valueFORCED = NA,
  ...
)
}
\arguments{
\item{object}{An object of class \code{survFit}.}

\item{data_predict}{A dataframe with three columns \code{time}, \code{conc} and \code{replicate}
used for prediction. If \code{NULL}, prediction is based on \code{x} object of 
class \code{survFit} used for fitting.}

\item{spaghetti}{If \code{TRUE}, return a set of survival curves using
parameters drawn from the posterior distribution.}

\item{mcmc_size}{Can be used to reduce the number of mcmc samples in order to speed up
the computation. \code{mcmc_size} is the number of selected iterations for one chain. Default
is 1000. If all MCMC is wanted, set argument to \code{NULL}.}

\item{hb_value}{If \code{TRUE}, the background mortality \code{hb} is taken into account from the posterior.
If \code{FALSE}, parameter \code{hb} is set to a fixed value. The default is \code{TRUE}.}

\item{interpolate_length}{Length of the time sequence for which output is wanted.}

\item{interpolate_method}{The interpolation method for concentration. See package \code{deSolve} for details.
Default is \code{linear}.}

\item{hb_valueFORCED}{If \code{hb_value} is \code{FALSE}, it fix \code{hb}.}

\item{\dots}{Further arguments to be passed to generic methods}
}
\value{
a \code{list} of \code{data.frame} with the quantiles of outputs in
\code{df_quantiles} or all the MCMC chaines \code{df_spaghetti}
}
\description{
This is the generic \code{predict} S3 method for the \code{survFit} class.
It provides predicted survival rate for "SD" or "IT" models under constant or time-variable exposure.
}
