\name{repro.fit}
\alias{repro.fit}
\alias{repro.parfit}
\alias{print.repro.fit}
\alias{summary.repro.fit}
\alias{plot.repro.fit}

\title{
Fit an exposure-response model to reproduction data taking mortality into account within the Bayesian framework
}
\description{
The \code{repro.fit} function estimates the parameters of an exposure-response model using Bayesian inference.
One deterministic part is proposed: the log-logistic function.
Two stochastic parts may be chosen by the function in order to take into account the nature of reproduction data
and the inter-replicate variability (Delignette-Muller et al., 2014).
The function calls the \code{rjags} package (Plummer, 2013)  and the \code{dclone} package for parallelization
of chains.
The function returns parameter estimates of the exposure-response model and estimates of \eqn{x} \% effective
concentration for \eqn{x} = 5, 10, 20 and 50.
The \code{repro.parfit} does the same thing as \code{repro.fit}, but chains are run on parallel workers,
so that computations can be faster for long MCMC runs.
Generic methods are \code{print}, \code{plot} and \code{summary}.
}
\usage{
repro.fit(rdata, n.chains = 3, quiet = FALSE)

repro.parfit(rdata, n.chains = 3, quiet = FALSE)

\method{summary}{repro.fit}(object, \dots)

\method{print}{repro.fit}(x, \dots)

\method{plot}{repro.fit}(x, xlab, ylab, fitcol, fitlty, fitlwd,
			ci = FALSE, cicol, cilty, cilwd, addlegend = TRUE, 
			log.scale = FALSE, type = "generic", \dots)
}
\arguments{
	\item{rdata}{An object of class \code{repro.data}.}
	\item{n.chains}{Number of MCMC chains. The minimum required number of chains is 2.}
	\item{quiet}{If \code{TRUE}, make silent all prints and progress bars of JAGS compilation.}
	\item{object}{An object of class \code{repro.fit}.}
	\item{x}{An object of class \code{repro.fit}.}
	\item{xlab}{A label for the \eqn{X}-axis, by default \code{Concentrations}.}
	\item{ylab}{A label for the \eqn{Y}-axis, by default \code{Response}.}
	\item{fitcol}{A single color to plot the fitted curve, by default \code{red}.}
	\item{fitlty}{A single line type to plot the fitted curve, by default \code{1}.}
	\item{fitlwd}{A single numeric which controls the width of the fitted curve, by default \code{1}.}
	\item{ci}{If \code{TRUE}, the 95 \% credible limits of the model are plotted.}
	\item{cicol}{A single color to plot the 95 \% credible limits, by default \code{red}.}
	\item{cilty}{A single line type to plot 95 \% credible limits, by default \code{1}.}
	\item{cilwd}{A single numeric which controls the width of the 95 \% credible limits, by default \code{2}.}
	\item{addlegend}{If \code{TRUE}, a default legend is added to the plot.}
	\item{log.scale}{If \code{TRUE}, a log-scale is used on the \eqn{X}-axis.}
	\item{type}{Graphical method: \code{generic} or \code{ggplot}.}
	\item{\dots}{Further arguments to be passed to generic methods.}
}
\details{
	\describe{
		\item{"log-logistic" deterministic part:}{
			the reproduction rate (expressed in number of offspring per individual-day) at concentration
			\eqn{C_{i}}{C_{i}} is described by:			
			\deqn{f(C_{i}) = \frac{d}{ 1 + (\frac{C_{i}}{e})^b}}{f(C_{i}) = d / (1 + (C_{i} / e)^b)}
			where \eqn{d} stands for the expected number of offspring per individual-day in the control,
			\eqn{e} is the 50 \% effective concentration (\eqn{EC_{50}}{EC50}) and
			\eqn{b} is a slope parameter.
			The number of offspring at concentration \eqn{i} and replicate \eqn{j} is described by a
			Poisson distribution of mean equal to the product of the number of individual-days \eqn{Nindtime_{ij}}{Nindtime ij}
			by a term \eqn{f_{ij}}{f ij} differing between the two stochastic parts.
			\deqn{N_{ij} \sim Poisson(f_{ij} \times Nindtime_{ij})}{N_{ij} ~ Poisson(f_{ij} * Nindtime_{ij})}
			with \eqn{Nindtime_{ij}} the number of individual-days at the target time for each replicate and
			each concentration (Delignette-Muller et al., 2014).
			
		}
		\item{"Poisson" stochastic part :}{
			\eqn{f_{ij}} only depends on the concentration:
			\deqn{f_{ij} = f(C_{i})}{f_{ij} = f(C_{i})}
		}
		\item{"Gamma-Poisson" stochastic part:}{
			\eqn{f_{ij}} is assummed to be variable between replicates at a same concentration and to follow a gamma
			distribution:
			\deqn{f_{ij} \sim Gamma(\frac{f(C_{i})}{\omega},\frac{1}{\omega})}{f_{ij} ~ Gamma(f(C_{i}) / omega, 1 / omega)}
			with \eqn{\omega}{omega} the overdispersion parameter.
		}
	}
	
	DIC:
	The Deviance Information Criterium (DIC) as defined by Spiegelhalter et al. (2002) is provided 
	by the \code{dic.samples} function.
	
	Raftery and Lewis's diagnostic:
	The \code{raftery.diag} is a run length control diagnostic based on a criterion that calculates the appropriate number of
	iterations required to accurately estimate the parameter quantiles. The Raftery and Lewis's diagnostic value
	used in the \code{repro.fit} function is the \code{resmatrix} object. See the \code{raftery.diag} help for more details. 
	
	Model selection:
	The \code{repro.fit} function chooses itself between the Poisson and the Gamma-Poisson model
	depending on the number of MCMC samples and on the DIC values.
	The minimum number of MCMC samples for the pilot run is provided by the Raftery
	and Lewis's diagnostic (Raftery and Lewis 1992). If this number is more than 100 000 or 
	if the DIC difference between models is small (typically less than 1), then the Poisson model is
	selected.
}
\value{
Returns an obbject of class \code{repro.fit}. A list of 13 objects:
  \item{DIC}{DIC value of the selected model.}
  \item{estim.ECx}{A table of the estimated 5, 10, 20 and 50 \% effective concentrations and their 95 \%
  credible intervals.}
  \item{estim.par}{A table of the estimated parameters as medians and 95 \% credible intervals.}
  \item{mcmc}{An object of class \code{mcmc.list} with the posterior distributions.}
  \item{model}{A JAGS model object.}
  \item{model.label}{An undocumented value for internaluse only.}
  \item{n.chains}{An integer value corresponding to the number of chains used for the MCMC computation.}
  \item{n.burnin}{A numerical value corresponding to the number of discarded draws for the burn-in period.}
  \item{n.iter}{A numerical value corresponding to the number of monitored iterations.}
  \item{param.prior}{An undocumented list for internaluse only.}
  \item{n.thin}{A numerical value corresponding to the thinning interval.}
  \item{raw.data}{The raw dataframe with five columns passed to the argument of \code{\link{repro.data}}.}
  \item{transformed.data}{A dataframe with six columns. See \code{\link{repro.data}} for details.}
 
  
Generic functions:
\describe{
	\item{\code{summary}}{
		provides the following information: 
		the type of model used, a summary of the MCMC chains with summary statistics
		for each variable: mean, standard deviation, naive standard error of the mean and time-series standard error 
		based on an estimate of the spectral density at 0, quantiles of the sample distribution using the \code{quantiles}
		argument. See \code{\link[coda]{summary.mcmc}}, prior quantiles, median and 2.5 \% and 97.5 \% quantiles of posterior
		distributions of estimated parameters and ECx estimates (x = 5, 10, 20, 50) as 50 \%, 2.5 \%  and 97.5 \% quantiles. 
	}
	\item{\code{print}}{
		shows information about the estimation method: the full JAGS model, the number
		of chains, the total number of iterations, the number of iterations in the burn-in period, the thin value
		and the DIC.
	}
	\item{\code{plot}}{
		shows the fitted exposure-response curve superimposed to experimental data at target time. The response is here
		expressed as the cumulative number of offspring per individual-day.
		See \code{\link{repro.data}}.
		Two types of output are available: \code{generic} or \code{ggplot}.
	}
}  
}
\references{
Delignette-Muller, M.L., Lopes, C., Veber, P. and Charles, S. (2014) Statistical handling of reproduction data for
exposure-response modelling. \url{http://arxiv.org/abs/1310.2733}.

Plummer, M. (2013) JAGS Version 3.4.0 user manual.
\url{http://sourceforge.net/projects/mcmc-jags/files/Manuals/3.x/jags_user_manual.pdf/download}

Raftery A.E. and Lewis, S.M. (1992) One long run with diagnostics: Implementation strategies for Markov chain Monte Carlo.
\emph{Statistical Science}, 7, 493-497.
 
Spiegelhalter, D., N. Best, B. Carlin, and A. van der Linde (2002) Bayesian measures of model complexity and fit (with discussion). 
\emph{Journal of the Royal Statistical Society}, Series B 64, 583-639. 
}
\author{
Marie Laure Delignette-Muller <marielaure.delignettemuller@vetagro-sup.fr>,
Philippe Ruiz <philippe.ruiz@univ-lyon1.fr>
}
\note{
When the \code{repro.parfit} function is used, the number of clusters is automatically defined by the function. It is equal to
argument \code{n.chains}.
}
\seealso{
\code{\link[rjags]{rjags}}, \code{\link[rjags]{coda.samples}}, \code{\link[rjags]{dic.samples}}, \code{\link[coda]{summary.mcmc}},
\code{\link[dclone]{parJagsModel}}, \code{\link[dclone]{parCodaSamples}}, \code{\link{repro.data}} and \code{\link[coda]{raftery.diag}}
}

\examples{
# (1) Load the data
data(cadmium1)

# (2) Create the transformed dataset
dat <- repro.data(cadmium1)
class(dat)

\dontrun{
# (3) Run the fit
out <- repro.fit(dat, quiet = TRUE)

# (4) Summary
summary(out)

# (5) Plot the fitted curve
plot(out, log.scale = TRUE, c.i. = TRUE)

# (6) Add a specific legend with generic type
plot(out, addlegend = FALSE)
legend("left", legend = c("Without mortality", "With mortality"),
pch = c(19,1)) 
}
}
\keyword{Estimation}