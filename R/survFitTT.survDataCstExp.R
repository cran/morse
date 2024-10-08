#' Fits a Bayesian concentration-response model for target-time survival analysis
#'
#' This function estimates the parameters of an concentration-response
#' model for target-time survival analysis using Bayesian inference. In this model,
#' the survival rate of individuals at a given time point (called target time) is modeled
#' as a function of the chemical compound concentration. The actual number of
#' surviving individuals is then modeled as a stochastic function of the survival
#' rate. Details of the model are presented in the
#' vignette accompanying the package.
#'
#' The function returns
#' parameter estimates of the concentration-response model and estimates of the so-called
#' \eqn{LC_x}, that is the concentration of chemical compound required to get an \eqn{(1 - x/100)} survival rate.
#'
#' @param data an object of class \code{survData}
#' @param target.time the chosen endpoint to evaluate the effect of the chemical compound
#' concentration, by default the last time point available for
#' all concentrations
#' @param lcx desired values of \eqn{x} (in percent) for which to compute
#' \eqn{LC_x}.
#' @param n.chains number of MCMC chains, the minimum required number of chains
#' is 2
#' @param quiet if \code{TRUE}, does not print messages and progress bars from
#' JAGS
#' @param \dots Further arguments to be passed to generic methods
#'
#' @return The function returns an object of class \code{survFitTT}, which is a
#' list with the following information:
#' \item{estim.LCx}{a table of the estimated \eqn{LC_x} along with their 95\%
#' credible intervals}
#' \item{estim.par}{a table of the estimated parameters (medians) and 95\%
#' credible intervals}
#' \item{det.part}{the name of the deterministic part of the used model}
#' \item{mcmc}{an object of class \code{mcmc.list} with the posterior
#' distribution}
#' \item{warnings}{a table with warning messages}
#' \item{model}{a JAGS model object}
#' \item{parameters}{a list of parameter names used in the model}
#' \item{n.chains}{an integer value corresponding to the number of chains used
#' for the MCMC computation}
#' \item{n.iter}{a list of two indices indicating the beginning and the end of
#' monitored iterations}
#' \item{n.thin}{a numerical value corresponding to the thinning interval}
#' \item{jags.data}{a list of the data passed to the JAGS model}
#' \item{transformed.data}{the \code{survData} object passed to the function}
#' \item{dataTT}{the dataset with which the parameters are estimated}
#'
#' @keywords estimation
#
#' @import rjags
#' @importFrom dplyr filter
#'
#' @export
survFitTT.survDataCstExp <- function(data,
                                     target.time = NULL,
                                     lcx = c(5, 10, 20, 50),
                                     n.chains = 3,
                                     quiet = FALSE,
                                     ...) {
  # test class object
  if(! is(data, "survDataCstExp"))
    stop("survFitTT: object of class survDataCstExp expected")

  # select Data at target.time and pool replicates
  dataTT <- selectDataTT(data, target.time)
  
  # Gather replicates according to time and conc
  dataTT <- cbind(aggregate(cbind(Nsurv, Ninit) ~ time + conc, dataTT, sum), replicate = 1)


  # Choose model by testing mortality in the control
  control <- filter(dataTT, conc == 0)
  det.part <-
    if (any(control$Nsurv < control$Ninit)) "loglogisticbinom_3"
  else "loglogisticbinom_2"

  # select model text
  if (det.part == "loglogisticbinom_2") {
    model.text <- llbinom2.model.text
  }
  if (det.part == "loglogisticbinom_3") {
    model.text <- llbinom3.model.text
  }

  # parameters
  parameters <- if (det.part == "loglogisticbinom_2") {
    c("log10b", "log10e")
  } else {
    if (det.part == "loglogisticbinom_3") {
      c("log10b", "d", "log10e")}
  }

  # create priors parameters
  jags.data <- survCreateJagsData(det.part, dataTT)

  # Define model
  model <- survLoadModel(model.program = model.text,
                         data = jags.data, n.chains,
                         Nadapt = 3000, quiet)

  # Determine sampling parameters
  sampling.parameters <- modelSamplingParameters(model,
                                                 parameters, n.chains, quiet)

  if (sampling.parameters$niter > 100000)
    stop("The model needs too many iterations to provide reliable parameter estimates !")

  # Sampling
  prog.b <- ifelse(quiet == TRUE, "none", "text")

  mcmc <- coda.samples(model, parameters,
                       n.iter = sampling.parameters$niter,
                       thin = sampling.parameters$thin,
                       progress.bar = prog.b)

  # summarize estime.par et CIs
  # calculate from the estimated parameters
  estim.par <- survPARAMS(mcmc, det.part)

  # LCx calculation  estimated LCx and their CIs 95%
  # vector of LCX
  estim.LCx <- estimXCX(mcmc, lcx, "LC")

  # check if estimated LC50 lies in the tested concentration range

  warnings <- msgTableCreate()

  LC50 <- log10(estim.par["e", "median"])
  if (!(min(log10(data$conc)) < LC50 & LC50 < max(log10(data$conc)))){
    ##store warning in warnings table
    msg <- "The LC50 estimation (model parameter e) lies outside the range of tested concentrations and may be unreliable as the prior distribution on this parameter is defined from this range !"
    warnings <- msgTableAdd(warnings, "LC50outRange", msg)
    ## print the message
    warning(msg, call. = FALSE)
  }

  # output
  OUT <- list(estim.LCx = estim.LCx,
              estim.par = estim.par,
              det.part = det.part,
              mcmc = mcmc,
              warnings = warnings,
              model = model,
              parameters = parameters,
              n.chains = summary(mcmc)$nchain,
              n.iter = list(start = summary(mcmc)$start,
                            end = summary(mcmc)$end),
              n.thin = summary(mcmc)$thin,
              jags.data = jags.data,
              transformed.data = data,
              dataTT = dataTT)

  class(OUT) <- "survFitTT"
  return(OUT)
}

survCreateJagsData <- function(det.part, data) {
  # Creates the parameters to define the prior of the log-logistic binomial model
  # INPUTS
  # det.part: model name
  # data: object of class survData
  # OUTPUT
  # jags.data : list of data required for the jags.model function

  # Parameter calculation of concentration min and max
  concmin <- min(sort(unique(data$conc))[-1])
  concmax <- max(data$conc)

  # Create prior parameters for the log logistic model

  # Params to define e
  meanlog10e <- (log10(concmin) + log10(concmax)) / 2
  sdlog10e <- (log10(concmax) - log10(concmin)) / 4
  taulog10e <- 1 / sdlog10e^2

  # Params to define b
  log10bmin <- -2
  log10bmax <- 2

  # list of data use by jags
  jags.data <- list(meanlog10e = meanlog10e,
                    Ninit = data$Ninit,
                    Nsurv = data$Nsurv,
                    taulog10e = taulog10e,
                    log10bmin = log10bmin,
                    log10bmax = log10bmax,
                    n = length(data$conc),
                    xconc = data$conc)

  # list of data use by jags
  if (det.part == "loglogisticbinom_3") {
    jags.data <- c(jags.data,
                   dmin = 0,
                   dmax = 1)
  }
  return(jags.data)
}

survPARAMS <- function(mcmc, det.part) {
  # create the table of posterior estimated parameters
  # for the survival analyses
  # INPUT:
  # - mcmc:  list of estimated parameters for the model with each item representing
  # a chain
  # OUTPUT:
  # - data frame with 3 columns (values, CIinf, CIsup) and 3-4rows (the estimated
  # parameters)

  # Retrieving parameters of the model
  res.M <- summary(mcmc)

  if (det.part ==  "loglogisticbinom_3") {
    d <- res.M$quantiles["d", "50%"]
    dinf <- res.M$quantiles["d", "2.5%"]
    dsup <- res.M$quantiles["d", "97.5%"]
  }
  # for loglogisticbinom_2 and 3
  b <- 10^res.M$quantiles["log10b", "50%"]
  e <- 10^res.M$quantiles["log10e", "50%"]
  binf <- 10^res.M$quantiles["log10b", "2.5%"]
  einf <- 10^res.M$quantiles["log10e", "2.5%"]
  bsup <- 10^res.M$quantiles["log10b", "97.5%"]
  esup <- 10^res.M$quantiles["log10e", "97.5%"]

  # Definition of the parameter storage and storage data
  # If Poisson Model

  if (det.part == "loglogisticbinom_3") {
    # if mortality in control
    rownames <- c("b", "d", "e")
    params <- c(b, d, e)
    CIinf <- c(binf, dinf, einf)
    CIsup <- c(bsup, dsup, esup)
  } else {
    # if no mortality in control
    # Definition of the parameter storage and storage data
    rownames <- c("b", "e")
    params <- c(b, e)
    CIinf <- c(binf, einf)
    CIsup <- c(bsup, esup)
  }

  res <- data.frame(median = params, Q2.5 = CIinf, Q97.5 = CIsup,
                    row.names = rownames)

  return(res)
}

llbinom3.model.text <- "\nmodel # Loglogistic binomial model with 3 parameters\n\t\t{\t\nfor (i in 1:n)\n{\np[i] <- d/ (1 + (xconc[i]/e)^b)\nNsurv[i]~ dbin(p[i], Ninit[i])\n}\n\n# specification of priors (may be changed if needed)\nd ~ dunif(dmin, dmax)\nlog10b ~ dunif(log10bmin, log10bmax)\nlog10e ~ dnorm(meanlog10e, taulog10e)\n\nb <- pow(10, log10b)\ne <- pow(10, log10e)\n}\n"

llbinom2.model.text <- "\nmodel # Loglogistic binomial model with 2 parameters\n\t\t{\t\nfor (i in 1:n)\n{\np[i] <- 1/ (1 + (xconc[i]/e)^b)\nNsurv[i]~ dbin(p[i], Ninit[i])\n}\n\n# specification of priors (may be changed if needed)\nlog10b ~ dunif(log10bmin, log10bmax)\nlog10e ~ dnorm(meanlog10e, taulog10e)\n\nb <- pow(10, log10b)\ne <- pow(10, log10e)\n}\n"
