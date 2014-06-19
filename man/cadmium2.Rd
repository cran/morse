\name{cadmium2}
\alias{cadmium2}
\docType{data}
\encoding{latin1}
\title{
Reproduction and survival datasets for snails exposed to cadmium during 56 days
}
\description{
Reproduction and survival datasets of chronic laboratory bioassays with snails exposed to six concentrations
of one metal contaminant (cadmium) during 56 days.
Six concentrations were tested, with six replicates per concentration. Each replicate contained five
organisms. Reproduction and survival were monitored at 17 time points. 
}
\usage{data(cadmium2)}
\format{
  A data frame with 612 observations of the following five variables:
  \describe{
    \item{\code{replicate}}{A vector of class \code{factor} with the replicate code (\code{A} to \code{F}).}
    \item{\code{conc}}{A vector of class \code{integer} with the cadmium concentrations in \eqn{\mu g.L^{-1}}.}
    \item{\code{time}}{A vector of class \code{integer} with the time points (in days from the beginning of the experiment
    					\eqn{t = 0}).}
    \item{\code{Nsurv}}{A vector of class \code{integer} with the number of alive individuals at each time point for each concentration
    					and each replicate.}
    \item{\code{Nrepro}}{A vector of class \code{integer} with the number of offspring at each time point for each concentration and each
    					replicate.}
  }
}
%\details{
%}
%\source{
%}
%\references{
%}
\examples{
# (1) Load the data
data(cadmium2)

# (2) Plot the number of survivors as a function of time for each concentration
repro.fullsurvplot(cadmium2)

# (3) Plot the number of survivors as a function of concentration
# at the end of the bioassay
repro.survplot(cadmium2, log.scale = TRUE)

# (4) Plot the cumulated number of offspring as a function of concentration   
repro.cumulplot(cadmium2, log.scale = TRUE)

}
\keyword{Datasets}
