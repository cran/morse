#' Fits a TKTD model for survival analysis using Bayesian inference
#'
#' This function estimates the parameters of a TKTD model ('SD' or 'IT')
#' for survival analysis using Bayesian inference. In this model,
#' the survival rate of individuals is modeled as a function of the chemical compound
#' concentration with a mechanistic description of the effects on survival over
#' time.
#' 
#' The function \code{survFit} returns the parameter estimates of Toxicokinetic-toxicodynamic (TKTD) models
#' \code{SD} for 'Stochastic Death' or \code{IT} fo 'Individual Tolerance'.
#' TKTD models, and particularly the General Unified Threshold model of
#' Survival (GUTS), provide a consistent process-based
#' framework to analyse both time and concentration dependent datasets.
#' In GUTS-SD, all organisms are assumed to have the same internal concentration 
#' threshold (denoted \eqn{z}), and, once exceeded, the instantaneous probability
#' to die increases linearly with the internal concentration.
#' In GUTS-IT, the threshold concentration is distributed among all the organisms, and once 
#' exceeded in one individual, this individual dies immediately.
#' 
#' 
#' When class of \code{object} is \code{survDataCstExp}, see \link[=survFit.survDataCstExp]{survFit.survDataCstExp} ;
#' and for a \code{survDataVarExp}, see \link[=survFit.survDataVarExp]{survFit.survDataVarExp}.
#'
#' @rdname survFit
#'
#' @param data An object of class \code{survDataCstExp} or \code{survDataVarExp}.
#' @param model_type Can be \code{"SD"} or \code{"IT"} to choose
#'   between "Stochastic Death" or "Individual Tolerance" models
#'   (resp.). See the modeling vignette for details.
#' @param quiet If \code{FALSE}, prints logs and progress bar from
#'   JAGS.
#' @param n.chains A positive integer specifying the number of MCMC chains. The minimum required number 
#' of chains is 2.
#' @param n.adapt A positive integer specifying the number of iterations for adaptation. If \code{n.adapt} = 0
#'  then no adaptation takes place.
#' @param n.iter A positive integer specifying the number of iterations to monitor for each chain.
#' @param n.warmup A positive integer specifying the number of warmup (aka burnin) iterations per chain. 
#' @param thin.interval A positive integer specifying the period to monitor.
#' @param limit.sampling if \code{FALSE} (default is \code{TRUE}), there is no limit to the number of iterations
#' in MCMC imposed by the \code{raftery.diag} test.
#' @param dic.compute if \code{TRUE} (default is \code{FALSE}), it generates penalized deviance samples to compute
#' the Deviance Information Criterion (DIC) with the \code{rjags} package
#' @param dic.type type of penalty to use. A string identifying the type of penalty: \code{pD} or \code{popt}
#'  (see function \code{\link[rjags]{dic.samples}})
#' @param hb_value If \code{TRUE}, the background mortality \code{hb} is taken into account.
#' If \code{FALSE}, parameter \code{hb} is set to 0. The default is \code{TRUE}.
#' @param  hb_valueFIXED If \code{hb_value} is \code{FALSE}, then \code{hb_valueFiXED} is the value to fix \code{hb}.
#'   If \code{hb_value} is \code{FALSE} and  \code{hb_valueFiXED} is \code{NA}, then \code{hb} is fixed to \code{0}.
#' @param extend_time Number of for each replicate used for linear 
#' interpolation (comprise between time to compute and fitting accuracy)
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @references Jager, T., Albert, C., Preuss, T. G. and Ashauer, R. (2011) 
#' General unified threshold model of survival-a toxicokinetic-toxicodynamic
#'  framework for ecotoxicology, \emph{Environmental Science and Technology}, 45, 2529-2540.
#' 303-314.
#' 
#' @keywords estimation 
#' 
#' @return an object of class \code{survFit}
#' 
#' @export
#' 
#' 
survFit <- function(data,
                    model_type,
                    quiet,
                    n.chains,
                    n.adapt,
                    n.iter,
                    n.warmup,
                    thin.interval,
                    limit.sampling,
                    dic.compute,
                    dic.type,
                    hb_value,
                    hb_valueFIXED,
                    ...){
  UseMethod("survFit")
}


################################################################################
#
#  PRIORS
#
################################################################################

#' Create a list of scalars giving priors to use in Bayesian inference.
#'
#' @param x An object of class \code{survData}
#' @param model_type TKTD model type ('SD' or 'IT')
#' 
#' @return A list for parameterization of priors for Bayesian inference with JAGS.
#' 
#' @export
#' 
priors_survData <- function(x, model_type = NULL){
  
  data <- filter(x, time != 0)
  
  # Parameter calculation of concentration min and max
  conc_min <- min(data$conc[data$conc != 0], na.rm = TRUE) # to remove 0 and NA
  conc_max <- max(data$conc, na.rm = TRUE)
  
  time_min <- min(data$time)
  time_max <- max(data$time)
  
  conc_unic <- sort(unique(data$conc))
  conc_unicPrec <- dplyr::lag(conc_unic)
  conc_minDelta <- min(conc_unic - conc_unicPrec, na.rm = TRUE)
  
  ##
  ## dominant rate constant: kd
  ##
  
  kd_max <- -log(0.001) / time_min
  kd_min <- -log(0.999) / time_max
  
  ##
  ## background hazard rate
  ##
  
  hb_max <- -log(0.5) / time_min
  hb_min <- -log(0.999) / time_max
  
  ##
  ## killing rate parameter: kk
  ##
  
  kk_max <- -log(0.001) / (time_min * conc_minDelta)
  kk_min <- -log(0.999) / (time_max * (conc_max - conc_min))
  
  ##
  ## beta
  ##
  
  beta_minlog10 <- -2
  beta_maxlog10 <- 2
  
  priorsMinMax <- list(
    conc_min = conc_min,
    conc_max = conc_max,
    
    kd_min = kd_min,
    kd_max = kd_max,
    
    hb_min = hb_min,
    hb_max = hb_max )
  
  ##
  ## Construction of the list of priors
  ##
  
  priorsList <-  list(
    ##
    ## dominant rate constant: kd
    ##
    kd_meanlog10 = (log10(kd_max) + log10(kd_min)) / 2 ,
    kd_sdlog10 = (log10(kd_max) - log10(kd_min)) / 4 ,
    ##
    ## background hazard rate
    ##
    hb_meanlog10 = (log10(hb_max) + log10(hb_min)) / 2 ,
    hb_sdlog10 = (log10(hb_max) - log10(hb_min)) / 4
  )
  
  if(model_type == "IT"){
    
    ## priorsMinMax
    priorsMinMax$beta_min <- beta_minlog10
    priorsMinMax$beta_max <- beta_maxlog10
    
    ## priorsList
    ### non effect threshold: scale parameter & median of a log-logistic distribution
    priorsList$alpha_meanlog10 <- (log10(conc_max) + log10(conc_min)) / 2
    priorsList$alpha_sdlog10 <- (log10(conc_max) - log10(conc_min)) / 4
    
    ### shape parameter of a log-logistic distribution
    priorsList$beta_minlog10 <- beta_minlog10
    priorsList$beta_maxlog10 <- beta_maxlog10
    
  } else if (model_type == "SD"){
    
    ## priorsMinMax
    priorsMinMax$kk_min <- kk_min
    priorsMinMax$kk_max <- kk_max
    
    ## priorsList
    ### killing rate parameter: kk
    priorsList$kk_meanlog10 <- (log10(kk_max) + log10(kk_min)) / 2
    priorsList$kk_sdlog10 <- (log10(kk_max) - log10(kk_min)) / 4
    ### non effect threshold: z
    priorsList$z_meanlog10 <- (log10(conc_max) + log10(conc_min)) / 2
    priorsList$z_sdlog10 <- (log10(conc_max) - log10(conc_min)) / 4
  } else stop("please, provide the 'model_type': 'IT' or 'SD'")
  
  
  return(list(priorsList = priorsList,
              priorsMinMax = priorsMinMax))
}


#############################################################################
#
#    survFit_TKTD_params
#
#############################################################################
  
survFit_TKTD_params <- function(mcmc, model_type, hb_value = TRUE) {
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
    
    kd <- 10^res.M$quantiles["kd_log10", "50%"]
    kd_inf95 <- 10^res.M$quantiles["kd_log10", "2.5%"]
    kd_sup95 <- 10^res.M$quantiles["kd_log10", "97.5%"]
    
    if(hb_value == TRUE){
      hb <- 10^res.M$quantiles["hb_log10", "50%"]
      hb_inf95 <- 10^res.M$quantiles["hb_log10", "2.5%"]
      hb_sup95 <- 10^res.M$quantiles["hb_log10", "97.5%"]
    }
    
    if(model_type == "SD"){
      kk <- 10^res.M$quantiles["kk_log10", "50%"]
      kk_inf95 <- 10^res.M$quantiles["kk_log10", "2.5%"]
      kk_sup95 <- 10^res.M$quantiles["kk_log10", "97.5%"]
      
      z <- 10^res.M$quantiles["z_log10", "50%"]
      z_inf95 <- 10^res.M$quantiles["z_log10", "2.5%"]
      z_sup95 <- 10^res.M$quantiles["z_log10", "97.5%"]
      
      if(hb_value == TRUE){
        res <- data.frame(parameters = c("kd", "hb", "z", "kk"),
                          median = c(kd, hb, z, kk),
                          Q2.5 = c(kd_inf95, hb_inf95, z_inf95, kk_inf95),
                          Q97.5 = c(kd_sup95, hb_sup95, z_sup95, kk_sup95))
      } else{
        res <- data.frame(parameters = c("kd", "z", "kk"),
                          median = c(kd, z, kk),
                          Q2.5 = c(kd_inf95, z_inf95, kk_inf95),
                          Q97.5 = c(kd_sup95, z_sup95, kk_sup95))
      }
      
    } else if (model_type == "IT"){
      alpha <- 10^res.M$quantiles["alpha_log10", "50%"]
      alpha_inf95 <- 10^res.M$quantiles["alpha_log10", "2.5%"]
      alpha_sup95 <- 10^res.M$quantiles["alpha_log10", "97.5%"]
      
      beta <- 10^res.M$quantiles["beta_log10", "50%"]
      beta_inf95 <- 10^res.M$quantiles["beta_log10", "2.5%"]
      beta_sup95 <- 10^res.M$quantiles["beta_log10", "97.5%"]
      
      if(hb_value == TRUE){
        res <- data.frame(parameters = c("kd", "hb", "alpha", "beta"),
                          median = c(kd, hb, alpha, beta),
                          Q2.5 = c(kd_inf95, hb_inf95, alpha_inf95, beta_inf95),
                          Q97.5 = c(kd_sup95, hb_sup95, alpha_sup95, beta_sup95))
      } else{
        res <- data.frame(parameters = c("kd", "alpha", "beta"),
                          median = c(kd, alpha, beta),
                          Q2.5 = c(kd_inf95, alpha_inf95, beta_inf95),
                          Q97.5 = c(kd_sup95, alpha_sup95, beta_sup95))
      }
    } else {
      stop("please, provide the 'model_type': 'IT' or 'SD'")
    }
    
    return(res)
}
