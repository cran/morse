% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survFit.R
\name{priors_survData}
\alias{priors_survData}
\title{Create a list of scalars giving priors to use in Bayesian inference.}
\usage{
priors_survData(x, model_type = NULL)
}
\arguments{
\item{x}{An object of class \code{survData}}

\item{model_type}{TKTD model type ('SD' or 'IT')}
}
\value{
A list for parameterization of priors for Bayesian inference with JAGS.
}
\description{
Create a list of scalars giving priors to use in Bayesian inference.
}
