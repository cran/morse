\name{repro.convergence}
\alias{repro.convergence}
\title{
Convergence check of the MCMC chains
}
\description{
The \code{repro.convergence} function checks the convergence of the MCMC chains from the JAGS estimate with the Gelman and Rubin convergence diagnostic (Gelman and Rubin, 1992).
It summarizes the \code{mcmc} or \code{mcmc.list} object with a trace of the sampled output,
a density estimate and an autocorrelation plot for each variable in the chain.}
\usage{
repro.convergence(out, trace = TRUE, density = TRUE,
				autocorr = TRUE, type = "generic")
}
\arguments{
  \item{out}{
An object of class \code{repro.fit}.
}
  \item{trace}{
If \code{TRUE}, the function traces the sampled output
estimate for each variable in the chain.
}
  \item{density}{
If \code{TRUE}, the function plots the density 
estimate for each variable in the chain.
}
  \item{autocorr}{
If \code{TRUE}, the function plots the autocorrelation for
each variable in each chain.
}
	\item{type}{Graphical method: \code{generic} or \code{ggplot}.}
}
%\details{
%}
\value{
A list with the point estimate of the multivariate potential scale reduction factor and the point
estimate of the potential scale reduction factor (Rhat) for each parameter of the Gelman and Rubin
test (Gelman and Rubin, 1992). A value close to 1 is expected when convergence is reached.
See the \code{gelman.diag} help for more details.
}
\references{
Gelman, A. and Rubin, D.B. (1992) \emph{Inference from iterative simulation using multiple sequences},
Statistical Science, 7, 457-511.
}
\author{
Marie Laure Delignette-Muller <marielaure.delignettemuller@vetagro-sup.fr>,
Philippe Ruiz <philippe.ruiz@univ-lyon1.fr>
}
\note{
When \code{type = "ggplot"}, the function calls packages \code{ggmcmc} and \code{gridExtra}
and returns an object of class \code{ggplot}.
}
\seealso{
\code{\link{repro.fit}} and \code{\link[coda]{gelman.diag}}, \code{\link[coda]{plot.mcmc}},
\code{\link[coda]{autocorr.plot}} from the \code{rjags} package and
\code{\link[ggmcmc]{ggs_traceplot}}, \code{\link[ggmcmc]{ggs_density}} and
\code{\link[ggmcmc]{ggs_autocorrelation}} from the \code{ggmcmc}
package (\url{http://xavier-fim.net/packages/ggmcmc})
}
\examples{
# (1) Load the data
data(zinc)

# (2) Create an object of class "repro.data"
dat <- repro.data(zinc)

\dontrun{
# (3) Run the repro.fit function
out <- repro.fit(dat)

# (4) Check the convergence
repro.convergence(out, trace = TRUE, density = FALSE, 
autocorr = TRUE)

# (5) Check the convergence using the "ggmcmc" package
repro.convergence(out, type = "ggplot")
}
}
\keyword{Model}
