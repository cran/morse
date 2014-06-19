\name{repro.fullsurvplot}
\alias{repro.fullsurvplot}
\title{
Plot of survival data
}
\description{
The \code{repro.fullsurvplot} function plots the number of survivors as a function of
time for each concentration and each replicate. This function is also used by \code{\link{repro.check.data}}.
}
\usage{
repro.fullsurvplot(data, xlab, ylab, type = "generic", addlegend = TRUE)
}
\arguments{
  \item{data}{
     Raw dataframe with five columns: \code{replicate}, \code{conc}, \code{time}, \code{Nsurv},
		\code{Nrepro}. See \code{\link{repro.data}} for details.
	}
  \item{xlab}{A label for the \eqn{X}-axis, by default \code{Time}.}
  \item{ylab}{A label for the \eqn{X}-axis, by default \code{Number of survivors}.}
	\item{type}{Graphical method: \code{generic}, \code{lattice} or \code{ggplot}.}
  \item{addlegend}{If \code{TRUE}, a legend is added to the plot.}
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
When \code{type = "lattice"}, the function returns an object of class \code{trellis}.
}
\seealso{
	\code{\link[ggplot2]{ggplot}}, \code{\link[lattice]{xyplot}}, \code{\link{repro.check.data}}
}
\examples{
# (1) Load the data
data(zinc)

# (2) Plot the survival data
repro.fullsurvplot(zinc, type = "generic", addlegend = TRUE)

# (3) Plot the survival data with a lattice type
repro.fullsurvplot(zinc, type = "lattice", addlegend = TRUE)

# (4) Plot the survival data with a ggplot type
repro.fullsurvplot(zinc, type = "ggplot", addlegend = FALSE)

# (5) To build a specific legend with a ggplot type
fu <- repro.fullsurvplot(zinc, type = "ggplot", addlegend = FALSE)
fu + theme(legend.position = "left") + scale_colour_hue("Replicate")

}
\keyword{Plot}
