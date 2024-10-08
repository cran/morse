% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotDoseResponse.survDataCstExp.R
\name{plotDoseResponse.survDataCstExp}
\alias{plotDoseResponse.survDataCstExp}
\title{Plot dose-response from \code{survData} objects}
\usage{
\method{plotDoseResponse}{survDataCstExp}(
  x,
  xlab = "Concentration",
  ylab = "Survival probability",
  main = NULL,
  target.time = NULL,
  style = "ggplot",
  log.scale = FALSE,
  remove.someLabels = FALSE,
  addlegend = TRUE,
  ...
)
}
\arguments{
\item{x}{an object of class \code{survData}}

\item{xlab}{a label for the \eqn{X}-axis, by default \code{Concentration}}

\item{ylab}{a label for the \eqn{Y}-axis, by default \code{Survival probability}}

\item{main}{main title for the plot}

\item{target.time}{a numeric value corresponding to some observed time in \code{data}}

\item{style}{graphical backend, can be \code{'ggplot'} or \code{'generic'}}

\item{log.scale}{if \code{TRUE}, displays \eqn{X}-axis in log-scale}

\item{remove.someLabels}{if \code{TRUE}, removes 75\% of X-axis labels in
\code{'ggplot'} style to avoid the label overlap}

\item{addlegend}{if \code{TRUE}, adds a default legend to the plot}

\item{\dots}{Further arguments to be passed to generic methods}
}
\value{
a plot of class \code{ggplot}
}
\description{
This is the generic \code{plotDoseResponse} S3 method for the \code{survData}
class. It plots the survival probability as a function of concentration at a given
target time.
}
\details{
The function plots the observed values of the survival probability at a given time point
as a function of concentration. The 95 \% binomial confidence interval is added
to each survival probability. It is calculated using function
\code{\link[stats]{binom.test}} from package \code{stats}.
Replicates are systematically pooled in this plot.
}
\note{
When \code{style = "generic"}, the function calls the generic function
\code{\link[graphics]{plot}}

When \code{style = "ggplot"}, the function return an object of class
 \code{ggplot}, see function \code{\link[ggplot2]{ggplot}}
}
\seealso{
\code{\link[stats]{binom.test}}
}
\keyword{plot}
