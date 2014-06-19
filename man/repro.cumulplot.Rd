\name{repro.cumulplot}
\alias{repro.cumulplot}
\title{
Cumulative plot of reproduction data
}
\description{
The \code{repro.cumulplot} function plots the cumulative number of offspring (e.g. eggs, clutches\dots) at the end
of the assay for a given concentration of contaminant.
Replicates with and without mortality are plotted with different symbols.
}
\usage{
repro.cumulplot(data, xlab, ylab, type = "generic",
			log.scale = FALSE, addlegend = TRUE, \dots)
}
\arguments{
	\item{data}{
		Raw dataframe with five columns: \code{replicate}, \code{conc}, \code{time}, \code{Nsurv},
		\code{Nrepro}. See \code{\link{repro.data}} for details.
	}
	\item{xlab}{
		A label for the \eqn{X}-axis, by default \code{Concentrations}.
	}
	\item{ylab}{
		A label for the \eqn{Y}-axis, by default \code{Nreprocumul}.
	}
	\item{type}{Graphical method: \code{generic} or \code{ggplot}.}
	\item{log.scale}{
		If \code{TRUE}, a log-scale is used on the \eqn{X}-axis.
	}
	\item{addlegend}{
		If \code{TRUE}, a default legend is added to the plot differentiating
		replicates with and without mortality.
	}
	\item{\dots}{
		Furthers arguments passed to the generic plot function.
	}
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
When \code{type = "ggplot"}, the function calls package \code{\link[ggplot2]{ggplot2}}
and returns an object of class \code{ggplot}.
}
\seealso{
 \code{\link[ggplot2]{ggplot}}
}
\examples{
# (1) Load the data
data(chlordan)

# (2) Plot the cummulative reproduction data
repro.cumulplot(chlordan, log.scale = TRUE)

# (3) Personalize legend 
repro.cumulplot(chlordan, addlegend = FALSE)
legend("left",title = "Legend", pch = c(19,1), bty = "n",
		legend = c("Without mortality","With mortality"))
		
# (4) Personalize legend with ggplot type
cum <- repro.cumulplot(chlordan, type = "ggplot", addlegend = FALSE)
cum + theme(legend.position = "left") + scale_colour_hue("Legend",
breaks = c("1","19"), labels = c("With mortality","Without mortality"))		
}
\keyword{Plot}

