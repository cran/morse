\name{morse-package}
\alias{morse-package}
\alias{morse}
\docType{package}
\title{
MOdelling tools for Reproduction and Survival data in Ecotoxicology
}
\description{
The package offers tools for ecotoxicologists and regulators based on advanced and innovative methods for a valuable
quantitative environmental risk assessment.
The package allows the analysis of bioassay reproduction data accounting for mortality all along the bioassay. 
Such data are commonly used to estimate Effective Concentration (\eqn{EC_{x}}) values from chronic toxicity tests. The aim is to
fit an exposure-reponse curve to reproduction data by Bayesian inference while taking into account mortality among
parents without loosing valuable data (Delignette-Muller et al., 2014).
Models are characterized by a deterministic log-logistic part associated with a stochastic part.
Two different stochastic parts can be chosen: Poisson or Gamma-Poisson.
The package uses the \code{rjags} package (Plummer, 2013), an interface from R to the JAGS
library for Bayesian data analysis.
Note that the \code{rjags} package does not include a copy of the JAGS library: you must install it separately.
For instructions on downloading JAGS, see the home page at \url{http://mcmc-jags.sourceforge.net}.
}
\details{
\tabular{ll}{
Package: \tab morse\cr
Type: \tab Package\cr
Version: \tab 1.0.1\cr
Date: \tab 2014-06-19\cr
License: \tab GPL (>=2)\cr
}
}
\author{
Marie Laure Delignette-Muller <marielaure.delignettemuller@vetagro-sup.fr>,
Philippe Ruiz <philippe.ruiz@univ-lyon1.fr>,
Sandrine Charles <sandrine.charles@univ-lyon1.fr>,
Wandrille Duchemin <wandrille.duchemin@insa-lyon.fr>,
Christelle Lopes <christelle.lopes@univ-lyon1.fr>,
Philippe Veber <philippe.veber@univ-lyon1.fr>

Maintainer: Philippe Ruiz <philippe.ruiz@univ-lyon1.fr>
}
\references{
Delignette-Muller, M.L., Lopes, C., Veber, P. and Charles, S. (2014) \emph{Statistical handling of reproduction data for
exposure-response modelling}. \url{http://arxiv.org/abs/1310.2733}.

Plummer, M. (2013) \emph{JAGS Version 3.4.0 user manual}.
\url{http://sourceforge.net/projects/mcmc-jags/files/Manuals/3.x/jags_user_manual.pdf/download}.
}
\keyword{package}
\seealso{
\code{\link[rjags]{rjags}}, \code{\link[lattice]{lattice}}, \code{\link[ggplot2]{ggplot}}
}
\examples{
# (1) Load the data
data(cadmium1)

# (2) Check data
repro.check.data(cadmium1)

# (3) Plot raw data
repro.survplot(cadmium1, log.scale = TRUE)
repro.fullsurvplot(cadmium1)
repro.cumulplot(cadmium1)

\dontrun{
# (4) Fit the log-logistic model
dat <- repro.data(cadmium1)
out <- repro.fit(dat)

# (5) Check the mcmc convergence
repro.convergence(out)

# (6) Summary of the results
plot(out)
summary(out)
print(out)
}
}

\keyword{Package}