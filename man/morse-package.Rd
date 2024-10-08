% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/morse.R
\name{morse-package}
\alias{morse-package}
\alias{morse}
\title{MOdelling tools for Reproduction and Survival data in Ecotoxicology}
\description{
Provides tools for the analysis of survival/reproduction
toxicity test data in quantitative environmental risk assessment. It can be
used to explore/visualize experimental data, and to get estimates
of \eqn{LC_{x}} (\eqn{X}\% Lethal Concentration) or,
\eqn{EC_{x}} (\eqn{X}\% Effective Concentration) by fitting exposure-response
curves. The \eqn{LC_{x}}, \eqn{EC_{x}} and parameters of the curve are
provided along with an indication of the uncertainty of the estimation.
\code{morse} can also be used to get an estimation of the \eqn{NEC} (No Effect Concentration)
by fitting a Toxico-Kinetic Toxico-Dynamic (TKTD) model (GUTS: General Unified Threshold
model of Survival). Within the TKTD-GUTS approach, \eqn{LC(x,t)}, \eqn{EC(x,t)} and \eqn{MF(x,t)}
(\eqn{x}\% Multiplication Factors aka Lethal Profiles) can be explored in proportion \eqn{x} and 
time \eqn{t}.
}
\details{
Estimation procedures in \code{morse} can be used without a deep knowledge of
their underlying probabilistic model or inference methods. Rather, they
were designed to behave as well as possible without requiring a user to
provide values for some obscure parameters. That said, \code{morse} models can also
be used as a first step to tailor new models for more specific situations.

The package currently handles survival and reproduction data. Functions
dedicated to survival (resp. reproduction) analysis start with a
\code{surv} (resp. \code{repro}) prefix. \code{morse} provides a similar
workflow in both cases:
\enumerate{
\item create and validate a data set
\item explore a data set
\item plot a data set
\item fit a model on a data set and output the expected estimates
\item check goodness of fit with posterior preditive check plot (ppc)
}

More specifically, for survival data handles with TKTD `GUTS` model, \code{morse}
provides:
\enumerate{
\item plot \eqn{LC(x,t)} and \eqn{MF(x,t)}.
\item compute goodness-of-fit measures (PPC percent, NRMSE and SPPE)
}

Those steps are presented in more details in the "Tutorial" vignette, while
a more formal description of the estimation procedures are provided in the
vignette called "Models in \code{morse} package". Please refer to these documents
for further introduction to the use of \code{morse}.

This reference manual is a detailed description of the functions exposed in
the package.

\strong{Getting started} The package uses the \code{rjags} package
(Plummer, 2013), an R interface to the JAGS library for Bayesian model
estimation. Note that the \code{rjags} package does not include a copy
of the JAGS library: you need to install it separately. For instructions
on downloading JAGS, see the home page at
\url{https://mcmc-jags.sourceforge.io}. Once done, simply follow the steps
described in the tutorial vignette.

\tabular{ll}{ Package: \tab morse\cr Type: \tab Package\cr Version: \tab
3.2.0\cr Date: \tab 2018-11-15\cr License: \tab GPL (>=2)\cr }
}
\references{
Delignette-Muller, M.L., Ruiz P. and Veber P. (2017)
\emph{Robust fit of toxicokinetic-toxicodynamic models using prior knowledge contained in the design of survival toxicity tests}.

Delignette-Muller, M.L., Lopes, C., Veber, P. and Charles, S. (2014)
\emph{Statistical handling of reproduction data for exposure-response modelling}.

Forfait-Dubuc, C., Charles, S., Billoir, E. and Delignette-Muller, M.L. (2012)
\emph{Survival data analyses in ecotoxicology: critical effect concentrations, methods and models. What should we use?}

Plummer, M. (2013) \emph{JAGS Version 4.0.0 user manual}.
\url{https://sourceforge.net/projects/mcmc-jags/files/Manuals/4.x/jags_user_manual.pdf/download}

Baudrot, V., Preux, S., Ducrot, V., Pavé, A. and Charles, S. (2018)
\emph{New insights to compare and choose TKTD models for survival based on an inter-laboratory study for Lymnaea stagnalis exposed to Cd}.

EFSA PPR Scientific Opinion (2018)
\emph{Scientific Opinion on the state of the art of Toxicokinetic/Toxicodynamic (TKTD) effect models for regulatory risk assessment of pesticides for aquatic organisms}
\url{https://www.efsa.europa.eu/en/efsajournal/pub/5377}.
}
\seealso{
\code{\link[rjags]{rjags}},
\code{\link[ggplot2]{ggplot2}}
}
\author{
Virgile Baudrot  <virgile.baudrot@posteo.net>,
Sandrine Charles <sandrine.charles@univ-lyon1.fr>,
Marie Laure Delignette-Muller <marielaure.delignettemuller@vetagro-sup.fr>,
Wandrille Duchemin <wandrille.duchemin@insa-lyon.fr>,
Benoit Goussen <Benoit.Goussen@ibacon.com>,
Guillaume Kon-Kam-king <guillaume.kon-kam-king@univ-lyon1.fr>,
Christelle Lopes <christelle.lopes@univ-lyon1.fr>,
Philippe Ruiz <philippe.ruiz@univ-lyon1.fr>,
Alexander Singer, <Alexander.Singer@rifcon.de> 
Philippe Veber <philippe.veber@univ-lyon1.fr>

Maintainer: Philippe Veber <philippe.veber@univ-lyon1.fr>
}
