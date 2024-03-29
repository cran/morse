% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_ode.survFit.R
\name{predict_ode}
\alias{predict_ode}
\title{Predict method for \code{survFit} objects}
\usage{
predict_ode(object, ...)
}
\arguments{
\item{object}{an object used to select a method \code{ppc}}

\item{\dots}{Further arguments to be passed to generic methods}
}
\value{
an object of class \code{predict_ode}
}
\description{
This is a \code{method} to replace function \code{predict} used on \code{survFit}
object when computing issues happen. \code{predict_ode} uses the \code{deSolve}
library to improve robustness. However, time to compute may be longer.
}
