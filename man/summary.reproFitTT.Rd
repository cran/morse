% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.reproFitTT.R
\name{summary.reproFitTT}
\alias{summary.reproFitTT}
\title{Summary of \code{reproFitTT} object}
\usage{
\method{summary}{reproFitTT}(object, quiet = FALSE, ...)
}
\arguments{
\item{object}{an object of class \code{reproFitTT}}

\item{quiet}{when \code{TRUE}, does not print}

\item{\dots}{Further arguments to be passed to generic methods}
}
\value{
The function returns a list with the following information:
\item{Qpriors}{quantiles of the model priors}
\item{Qposteriors}{quantiles of the model posteriors}
\item{QECx}{quantiles of ECx estimates}
}
\description{
This is the generic \code{summary} S3 method for the \code{reproFitTT} class.
It shows the quantiles of priors and posteriors on parameters
and the quantiles of the posterior on the ECx estimates.
}
\keyword{summary}
