\name{cadmium1}
\alias{cadmium1}
\docType{data}
\encoding{latin1}
\title{
Reproduction and survival datasets for \emph{Daphnia magna} exposed to cadmium during 21 days
}
\description{
Reproduction and survival datasets of chronic laboratory bioassays with \emph{Daphnia magna} freshwater invertebrate
exposed to five concentrations of one metal contaminant (cadmium) during 21 days.
Five concentrations were tested, with four replicates per concentration. Each replicate contained 10
organisms. Reproduction and survival were monitored at 10 time points. 
}
\usage{data(cadmium1)}
\format{
  A data frame with 200 observations of the following five variables:
  \describe{
    \item{\code{replicate}}{A vector of class \code{integer} with the replicate code (\code{1} to \code{4}).}
    \item{\code{conc}}{A vector of class \code{numeric} with the cadmium concentrations in \eqn{\mu g.L^{-1}}.}
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
\references{
Billoir, E., Delhaye, H., Forfait, C., Cl\enc{é}{e}ment, B., Triffault-Bouchet, G., Charles, S. and 
Delignette-Muller, M.L. (2012) Comparison of bioassays with different exposure time patterns: The added value of dynamic modelling in predictive ecotoxicology,
\emph{Ecotoxicology and Environmental Safety}, 75, 80-86.
}
\examples{
# (1) Load the data
data(cadmium1)

# (2) Plot the number of survivors as a function of time for each concentration
repro.fullsurvplot(cadmium1)

# (3) Plot the number of survivors as a function of concentration
# at the end of the bioassay
repro.survplot(cadmium1, log.scale = TRUE)

# (4) Plot the cumulated number of offspring as a function of concentration  
repro.cumulplot(cadmium1, log.scale = TRUE)

}
\keyword{Datasets}
