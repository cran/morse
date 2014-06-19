\name{repro.data}
\alias{repro.data}
\alias{summary.repro.data}
\alias{print.repro.data}

\title{
Transformed dataset for the \code{repro.fit} function
}
\description{
The \code{repro.data} function creates a \code{repro.data} object needed to run the
\code{\link{repro.fit}} function. A new dataframe called \code{transformed.data} is created by 
adding three new columns at target time for each replicate and each concentration: the initial
number of individuals (\code{Ninit}), the cumulative number of offspring (\code{Nreprocumul})
and the number of individual-days (\code{Nindtime}).
The generic methods are \code{print} and \code{summary}.
The use of the \code{\link{repro.check.data}} function is recommended before using the
\code{\link{repro.data}} function to detect possible structural errors in raw data.
}
\usage{
repro.data(data, target.time)
\method{summary}{repro.data}(object, \dots)
\method{print}{repro.data}(x, \dots)
}

\arguments{
\item{data}{The raw dataframe with five columns passed to the function in argument \code{data}:
 	\describe{
		\item{replicate}{A vector of class \code{integer} or factor for replicate identification.}
		\item{conc}{A vector of class \code{numeric} with tested concentrations (positive values).}
		\item{time}{A vector of class \code{integer} with time points (positive values).
					The first time must be 0.}
		\item{Nsurv}{A vector of class \code{integer} with positive values of the number of alive individuals
					(positive values) at each time point for each concentration and each replicate.}
		\item{Nrepro}{A vector of class \code{integer} with the number of offspring (positive values)
						at each	time point for each concentration and each replicate.}
}}
\item{target.time}{
	The time at which the number of individual-days and the cumulative number of offspring from the beginning
	of the bioassay are calculated. By default the last time point.
}
\item{object}{An object of class \code{repro.data}.}
\item{x}{An object of class \code{repro.data}.}
\item{\dots}{Further arguments to be passed to generic methods.}
}
\details{
The raw dataframe must have exactly five columns arranged in the correct order with the following headings: 
\code{replicate}, \code{conc}, \code{time}, \code{Nsurv} and \code{Nrepro}. It is recommended to use the function
\code{\link{repro.check.data}} in order to check if the dataframe is in the correct format.
The \code{repro.data} function builds a new dataframe called \code{transformed.data}, that is a subset of the raw dataframe
at \code{target.time} with three additional columns: \code{Ninit}, \code{Nreprocumul}, \code{Nindtime}.
The number of individual-days is needed to account for the time-contribution of each individual to the cummulative reproduction
(Delignette-Muller et al., 2014).
See the section \bold{Value}. The new dataframe is automatically reordered by \code{time}, \code{concentration}
and \code{replicate}.
}
\value{
Returns an object of class \code{repro.data}. A list of three objects:
	\item{raw.data}{The raw dataframe with five columns corresponding to the argument passed in the function.}
	\item{transformed.data}{A dataframe with six columns:
	\describe{ 
		\item{replicate}{A vector of class \code{integer} or factor for replicate identification.}
		\item{conc}{A vector of class \code{numeric} with tested concentrations (positive values).}
		\item{Ninit}{A vector of class \code{integer} with the number of individuals
					at the beginning of the bioassay (positive values).}
		\item{Nsurv}{A vector of class \code{integer} with positive values of the number of alive individuals (positive values)
					at the target time for each concentration and each replicate.}
		\item{Nreprocumul}{A vector of class \code{integer} with the cumulative number of offspring (positive values)
							at the target time for each concentration and each replicate.}
		\item{Nindtime}{A vector of class \code{numeric} with the number of individual-days (positive values) at 
						target time for each concentration and each replicate.}
	}}
	\item{target.time}{The time at which the number of individual-days and the cumulative number of offspring from the
						beginning of the bioassay are calculated.}

Generic functions:
\describe{
	\item{\code{summary}}{
		The summary provides information about the structure of the dataset and the experimental design: 
		the number of datapoints per replicate, concentration and time both for the raw dataset and the 
		transformed dataset.
	}
	\item{\code{print}}{
		Print of a \code{repro.data} object with the transformed dataframe
		and the value of the target time.
	}
}}

\references{
Delignette-Muller M.L., Lopes C., Veber P. and Charles S. (2014) \emph{Statistical handling of reproduction data for
exposure-response modelling}. \url{http://arxiv.org/abs/1310.2733}.
}
\author{
Marie Laure Delignette-Muller <marielaure.delignettemuller@vetagro-sup.fr>,
Philippe Ruiz <philippe.ruiz@univ-lyon1.fr>
}
%\note{
%}
\seealso{
\code{\link{repro.check.data}}, \code{\link{repro.fit}}
}
\examples{
# (1) Load the data
data(zinc)

# (2) Create an objet of class 'repro.data'
dat <- repro.data(zinc)

# (3) Print and summarize object dat
print(dat)
summary(dat)
}

\keyword{Model}

