% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/priors_distribution.R
\name{priors_distribution}
\alias{priors_distribution}
\title{Density distribution of priors.}
\usage{
priors_distribution(object, ...)
}
\arguments{
\item{object}{An object used to select a method}

\item{\dots}{Further arguments to be passed to generic methods}
}
\value{
an object of class \code{priors_distribution}
}
\description{
Return a \code{data.frame} with prior density distributions of parameters used in
\code{object}.
}
\details{
When the \code{object} is of class \code{survFit}, see \link[=priors_distribution.survFit]{priors_distribution.survFit}
}
