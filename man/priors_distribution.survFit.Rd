% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/priors_distribution.R
\name{priors_distribution.survFit}
\alias{priors_distribution.survFit}
\title{Density distribution of priors from a \code{survFit} object.}
\usage{
\method{priors_distribution}{survFit}(object, size_sample = 1000, EFSA_name = FALSE, ...)
}
\arguments{
\item{object}{An object of class \code{survFit}.}

\item{size_sample}{Size of the random generation of the distribution.
Default is \code{1e3}.}

\item{EFSA_name}{If \code{TRUE}, replace the current terminology by
the one used in the recent EFSA PPR Scientific Opinion (2018).}

\item{\dots}{Further arguments to be passed to generic methods.}
}
\value{
a \code{data.frame} with prio distribution.
}
\description{
Return a \code{data.frame} with priors distribution of parameters used in
\code{object}.
}
\references{
EFSA PPR Scientific Opinion (2018)
\emph{Scientific Opinion on the state of the art of Toxicokinetic/Toxicodynamic (TKTD) effect models for regulatory risk assessment of pesticides for aquatic organisms}
\url{https://www.efsa.europa.eu/en/efsajournal/pub/5377}.
}
