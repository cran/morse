\name{chlordan}
\alias{chlordan}
\docType{data}
\encoding{latin1}
\title{
Reproduction and survival datasets for \emph{Daphnia magna} exposed to chlordan during 21 days
}
\description{
Reproduction and survival datasets of chronic laboratory bioassays with \emph{Daphnia magna} freshwater invertebrate
exposed to six concentrations of one organochlorine insecticide during 21 days.
Six concentrations were tested, with 10 replicates per concentration. Each replicate contained one
organism. Reproduction and survival were monitored at 22 time points. 
}
\usage{data(chlordan)}
\format{
  A data frame with 1320 observations of the following five variables:
  \describe{
    \item{\code{replicate}}{A vector of class \code{integer} with the replicate code (\code{1} to \code{10}).}
    \item{\code{conc}}{A vector of class \code{numeric} with the chlordan concentrations in \eqn{\mu g.L^{-1}}.}
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
Manar, R., Bessi, H. and Vasseur, P. (2009) Reproductive effects and bioaccumulation of chlordan in Daphnia magna, 
\emph{Environmental Toxicology and Chemistry}, 28, 2150-2159.
}
\examples{
# (1) Load the data
data(chlordan)

# (2) Plot the number of survivors as a function of time for each concentration
repro.fullsurvplot(chlordan)

# (3) Plot the number of survivors as a function of concentration
# at the end of the bioassay
repro.survplot(chlordan, log.scale = TRUE)

# (4) Plot the cumulated number of offspring as a function of concentration  
repro.cumulplot(chlordan, log.scale = TRUE)
}
\keyword{Datasets}