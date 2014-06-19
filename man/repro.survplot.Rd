\name{repro.survplot}
\alias{repro.survplot}
\title{
Plot of survival data
}
\description{
The \code{repro.survplot} plots the number of initially present animals still alive at the
end of the assay for a given concentration of contaminant. There is a single point by concentration and replicate.
The size of points is proportional to the number of overploted replicates.
}
\usage{
repro.survplot(data, xlab, ylab, pch, type = "generic",
				log.scale = FALSE, addlegend = TRUE, \dots)
}

\arguments{
\item{data}{
Raw dataframe with five columns: \code{replicate}, \code{conc}, \code{time}, \code{Nsurv},
\code{Nrepro}. See \code{\link{repro.data}} for details.
}
\item{xlab}{A label for the \eqn{X}-axis, by default \code{Concentrations}.}
\item{ylab}{A label for the \eqn{Y}-axis, by default \code{Number of survivors}.}
\item{pch}{Argument to choose the symbol of replicates.}
\item{type}{Graphical method: \code{generic} or \code{ggplot}.}
\item{log.scale}{If \code{TRUE}, a log-scale is used on the \eqn{X}-axis.}
\item{addlegend}{If \code{TRUE}, a default legend is added to the plot.}
\item{\dots}{Further arguments pass to the generic plot function.}
}
%\details{
%}
%\value{
%}
%\references{
%}
\author{
Marie Laure Delignette-Muller <marielaure.delignettemuller@vetagro-sup.fr>,
Philippe Ruiz <philippe.ruiz@univ-lyon1.fr>
}
\note{
The argument \code{pch} is only used with \code{type = "generic"}.
When \code{type = "ggplot"}, the function calls package \code{\link[ggplot2]{ggplot2}}
and returns an object of class \code{ggplot}.
}
\seealso{
	\code{\link[ggplot2]{ggplot}}
}
\examples{
# (1) Load the data
data(zinc)

# (2) plot the number of survivors depending on the concentration
repro.survplot(zinc, log.scale = TRUE)

# (3) To build a specific legend
repro.survplot(zinc, type = "generic", addlegend = FALSE)
legend("left", legend = "Replicate", pch = 4)

# (4) To build a specific legend with ggplot type
surv <- repro.survplot(zinc, type = "ggplot", addlegend = FALSE)
surv + theme(legend.position = "left")	
}
\keyword{Plot}
