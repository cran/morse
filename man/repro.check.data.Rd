\name{repro.check.data}
\alias{repro.check.data}
\alias{print.repro.check.data}
\title{
Check the consistency of the dataset
}
\description{
The \code{repro.check.dat} function performs several tests on the integrity of the dataset (column headings, type of data\dots) and returns
a list of error messages if data are not in the correct format. The aim of this function is to check the consistency 
of the dataframe before using function \code{\link{repro.data}}. This function highlights possible errors in the data structure
that would disturb or prevent the execution of the function \code{\link{repro.data}}. 
}
\usage{
repro.check.data(data, diagnos.plot = TRUE)

\method{print}{repro.check.data}(x, \dots)
}
\arguments{
  \item{data}{Raw dataframe with five columns. See \code{\link{repro.data}} function for details on the required
  				data format.
  }
  \item{diagnos.plot}{
  	If \code{TRUE}, calls the default \code{\link{repro.fullsurvplot}} function if the number of survivors increases
  	at some time points.
  }
  \item{x}{An object of class repro.check.data.}
  \item{\dots}{Further arguments to be passed to generic methods.}
}
\details{
For a given dataframe, the function checks if: 
	\describe{
		\item{1)}{column headings are correct: \code{replicate} for the column of replicates,
					\code{conc} for the column of concentrations, \code{time} for the column of
					time points, \code{Nsurv} for the column of the number of alive individuals and
					\code{Nrepro} for the column of the number of collected offspring at each time point,}
		\item{2)}{the first time point of the dataset is 0,}	
		\item{3)}{the class of column \code{conc} is \code{numeric},}
		\item{4)}{the classes of columns \code{Nsurv} and \code{Nrepro} are \code{integer},}
		\item{5)}{values of the dataframe are all positive,}	
		\item{6)}{the number of collected offspring is 0 at \eqn{t = 0},}
		\item{7)}{each replicate appears only once per concentration and per time point,}
		\item{8)}{the number of replicates is the same at any concentration and any time point,}
		\item{8)}{the number of alive individuals never increases with time,}
		\item{9)}{at each time \eqn{T}, if the number of alive individuals is null,
					 the number of collected offspring is also null at time \eqn{T+1}.}
	}
}
\value{
Returns an object of class \code{repro.check.data}. A dataframe with two columns of character string, \code{id} and \code{msg}.
The \code{id} is invisible when displaying the function. Print only shows error
messages \code{msg}.
\item{id}{The identifier of the test, equals to: 
	\describe{
		\item{\code{missingcolumn}}{if one or more columns are missing or if the column headings
										are not \code{replicate}, \code{conc}, \code{time},
										 \code{Nsurv} and \code{Nrepro}.}
		\item{\code{firstTime0}}{if the first time point is not 0 at each concentration and each replicate.}
		\item{\code{concNumeric}}{if column \code{conc} does not contain values of class \code{numeric} only.}
		\item{\code{NsurvInteger}}{if column \code{Nsurv} does not contain values of class \code{integer} only.}
		\item{\code{NreproInteger}}{if column \code{Nrepro} does not contain values of class \code{integer} only.}
		\item{\code{tablePositive}}{if there are negative values within the data.}
		\item{\code{Nrepro0T0}}{if \code{Nrepro} is not 0 at time 0 for each concentration and each replicate.}
		\item{\code{onlyReplicate}}{if a replicate is duplicated on different lines for the same time points and the same concentration.}
		\item{\code{missingReplicate}}{if a replicate is missing for at least one time points at one concentration.} 
		\item{\code{NsurvMonotone}}{if \code{Nsurv} increases at some time points compared to the previous one.}
		\item{\code{Nsurvt0Nreprotp1P}}{if at a giving time \eqn{T}, the number of alive individuals is null and 
											the number of collected offspring is not null for the same replicate and the same
								 			concentration at time \eqn{T+1}.}
	}
}		
\item{msg}{One or more user friendly error messages are generated:
	\describe{
		\item{}{\code{The column 'colname' is missing or have a wrong name.}}
		\item{}{\code{Data are required at time 0 for each concentration and each
						replicate.}}
		\item{}{\code{Column 'conc' must contain only numerical values.}}
		\item{}{\code{Column 'Nsurv' must contain only integer values.}}
		\item{}{\code{Column 'Nrepro' must contain only integer values.}}
		\item{}{\code{Data must contain only positive values.}}
		\item{}{\code{'Nrepro' should be 0 at time 0 for each concentration and each
						replicate.}}
		\item{}{\code{Replicate 'replicate' appears
						on different lines for the same time point at concentration 'conc'.}}
		\item{}{\code{Replicate 'replicate' is missing for at least one time points at concentration 'conc'.}}
		\item{}{\code{For replicate 'replicate' and concentration 'conc', 'Nsurv'
						increases at some time points compared to the previous one.}}
		\item{}{\code{For replicate 'replicate' and concentration 'conc', there are some 
						'Nsurv' = 0 followed by 'Nrepro' > 0 at the next time point.}}
}}
}
%\references{
%}
\author{
Marie Laure Delignette-Muller <marielaure.delignettemuller@vetagro-sup.fr>,
Philippe Veber <philippe.veber@univ-lyon1.fr>,
Philippe Ruiz <philippe.ruiz@univ-lyon1.fr>
}
\note{
If an error of type \code{missingcolumn} is detected, the function \code{repro.check.data} is stopped.
When no error is detected the \code{repro.check.data} function returns an empty dataframe.
}
\seealso{
\code{\link{repro.fullsurvplot}}, \code{\link{repro.data}}
}
\examples{
# Run the check data function
data(zinc)
repro.check.data(zinc)

# Example with an error in the dataframe

# (1) Load the data
data(zinc)

# (2) Insert an error (increase the number of survivors at a certain time point compared to its
# value at the previous time point within the same remplicate)
zinc[25,"Nsurv"] <- 20
check <- repro.check.data(zinc, diagnos.plot = TRUE)

# (3) Check for potential errors in the dataframe
check
}
\keyword{Check data}
