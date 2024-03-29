% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.survFitTKTD.R
\name{plot.survFitTKTD}
\alias{plot.survFitTKTD}
\title{Plotting method for \code{survFitTKTD} objects}
\usage{
\method{plot}{survFitTKTD}(
  x,
  xlab = "Time",
  ylab = "Survival probablity",
  main = NULL,
  concentration = NULL,
  spaghetti = FALSE,
  one.plot = FALSE,
  adddata = FALSE,
  addlegend = FALSE,
  style = "ggplot",
  ...
)
}
\arguments{
\item{x}{An object of class \code{survFitTKTD}.}

\item{xlab}{A label for the \eqn{X}-axis, by default \code{Time}.}

\item{ylab}{A label for the \eqn{Y}-axis, by default \code{Survival probablity}.}

\item{main}{A main title for the plot.}

\item{concentration}{A numeric value corresponding to some specific concentration in
\code{data}. If \code{concentration = NULL}, draws a plot for each concentration.}

\item{spaghetti}{if \code{TRUE}, draws a set of survival curves using
parameters drawn from the posterior distribution}

\item{one.plot}{if \code{TRUE}, draws all the estimated curves in
one plot instead of one plot per concentration.}

\item{adddata}{if \code{TRUE}, adds the observed data to the plot
with (frequentist binomial) confidence intervals}

\item{addlegend}{if \code{TRUE}, adds a default legend to the plot.}

\item{style}{graphical backend, can be \code{'generic'} or \code{'ggplot'}}

\item{\dots}{Further arguments to be passed to generic methods.}
}
\value{
a plot of class \code{ggplot}
}
\description{
This is the generic \code{plot} S3 method for the
\code{survFitTKTD}.  It plots the fit obtained for each
concentration of chemical compound in the original dataset.
}
\details{
The fitted curves represent the \strong{estimated survival probablity} as a function
of time for each concentration
When \code{adddata = TRUE} the black dots depict the \strong{observed survival
probablity} at each time point. Note that since our model does not take
inter-replicate variability into consideration, replicates are systematically
pooled in this plot.
The function plots both 95\% credible intervals for the estimated survival
probablity (by default the grey area around the fitted curve) and 95\% binomial confidence
intervals for the observed survival probablity (as black error bars if
\code{adddata = TRUE}).
Both types of intervals are taken at the same level. Typically
a good fit is expected to display a large overlap between the two types of intervals.
If \code{spaghetti = TRUE}, the credible intervals are represented by two
dotted lines limiting the credible band, and a spaghetti plot is added to this band.
This spaghetti plot consists of the representation of simulated curves using parameter values
sampled in the posterior distribution (2\% of the MCMC chains are randomly
taken for this sample).
}
\examples{

# (1) Load the survival data
data(propiconazole)

# (2) Create an object of class "survData"
dataset <- survData(propiconazole)

\donttest{
# (3) Run the survFitTKTD function ('SD' model only)
out <- survFitTKTD(dataset)

# (4) Plot the fitted curves in one plot 
plot(out)

# (5) Plot one fitted curve per concentration with credible limits as
# spaghetti, data and confidence intervals
# and with a ggplot style
plot(out, spaghetti = TRUE , adddata = TRUE, one.plot = FALSE,
     style = "ggplot")

# (6) Plot fitted curve for one specific concentration
plot(out, concentration = 36, style = "ggplot")
}

}
\keyword{plot}
