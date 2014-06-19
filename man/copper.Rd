\name{copper}
\alias{copper}
\docType{data}
\encoding{latin1}
\title{
Reproduction and survival datasets for \emph{Daphnia magna} exposed to copper during 21 days
}
\description{
Reproduction and survival datasets of chronic laboratory bioassays with \emph{Daphnia magna} 
freshwater invertebrate exposed to five concentrations of one metal contaminant (copper) during 21 days.
Five concentrations were tested, with three replicates per concentration. Each replicate contained 20
organisms. Reproduction and survival were monitored at 16 time points. 
}
\usage{data(copper)}
\format{
  A data frame with 240 observations of the following five variables:
  \describe{
    \item{\code{replicate}}{A vector of class \code{factor} with the replicate code (\code{A} to \code{C}).}
    \item{\code{conc}}{A vector of class \code{numeric} with the copper concentrations in \eqn{\mu g.L^{-1}}.}
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
Billoir, E., Delignette-Muller, M.L., \enc{Péry}{Pery}, A.R.R. and Charles, S. (2008)
A Bayesian Approach to Analyzing Ecotoxicological Data, \emph{Environmental Science & Technology}, 42 (23), 8978-8984.
}
\examples{
# (1) Load the data
data(copper)

# (2) Plot the number of survivors as a function of time for each concentration
repro.fullsurvplot(copper)

# (3) Plot the number of survivors as a function of concentration
# at the end of the bioassay
repro.survplot(copper, log.scale = TRUE)

# (4) Plot the cumulated number of offspring as a function of concentration  
repro.cumulplot(copper, log.scale = TRUE)
}
\keyword{Datasets}

