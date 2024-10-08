#' Posterior predictive check plot for \code{survFitCstExp} objects
#'
#' This is the generic \code{ppc} S3 method for the \code{survFitCstExp} class. It
#' plots the predicted values along with 95\% credible intervals
#' versus the observed values for \code{survFit} objects.
#'
#' The black points show the observed number of survivors (pooled
#' replicates, on \eqn{X}-axis) against the corresponding predicted
#' number (\eqn{Y}-axis). Predictions come along with 95\% prediction
#' intervals, which are depicted in green when they contain the
#' observed value and in red otherwise. Samples with equal observed
#' value are shifted on the \eqn{X}-axis. For that reason, the
#' bisecting line (y = x), is represented by steps when observed
#' values are low. That way we ensure green intervals do intersect the
#' bisecting line.
#'
#' @rdname PPC
#' 
#' @param x An object of class \code{survFitCstExp}
#' @param style graphical backend, can be \code{'generic'} or \code{'ggplot'}
#' @param main main title for the plot
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @return a plot of class \code{ggplot}
#'
#' @import ggplot2
#' @import grDevices
#' @importFrom graphics plot
#'
#' @export
#' 
ppc.survFitCstExp <- function(x, style = "ggplot", main = NULL, ...) {

  xlab <- "Observed nb of survivors"
  ylab <- "Predicted nb of survivors"
  
  ppc_gen(EvalsurvTKTDPpc_CstExp(x), style, xlab, ylab, main)
}

#' @importFrom stats rbinom quantile plogis
EvalsurvTKTDPpc_CstExp <- function(x) {
  tot.mcmc <- do.call("rbind", x$mcmc)

  model_type <- x$model_type
  
  kd <- 10^(tot.mcmc[, "kd_log10"])
  # "hb" is not in survFit object of morse <v3.2.0
  if("hb" %in% colnames(tot.mcmc)){
    hb <- tot.mcmc[, "hb"]  
  } else{ hb <- 10^tot.mcmc[, "hb_log10"] }

  if(model_type == "SD"){
    z <- 10^(tot.mcmc[, "z_log10"])
    kk <- 10^(tot.mcmc[, "kk_log10"])
  } else if (model_type == "IT"){
    alpha <- 10^(tot.mcmc[, "alpha_log10"])
    beta <- 10^(tot.mcmc[, "beta_log10"])
  } else{
    stop("'model_type' must be 'SD' or 'IT'")
  }
  
  #NsurvObs <- x$jags.data$y
  NsurvObs <- x$jags.data$Nsurv
  
  #n <- x$jags.data$ndat
  n <- x$jags.data$n_data
  
  #xconc <- x$jags.data$x
  xconc <- x$jags.data$conc
  
  #t <- x$jags.data$t
  time <- x$jags.data$time
  
  Nprec <- x$jags.data$Nprec
  
  if(model_type == "SD"){
    
    niter <- nrow(tot.mcmc)
   
    tprec <- x$jags.data$tprec
    
    NsurvPred <- matrix(NA, nrow = niter, ncol = n)
    psurv = NULL
    
    # bigtime <- x$jags.data$bigtime
    bigtime <- max(time) + 10
    
    for (i in 1:n) {
      for (j in 1:length(kd)) {
        xcor <- ifelse(xconc[i] > 0, xconc[i], 10)
        R <- ifelse(xconc[i] > z[j], z[j]/xcor, 0.1)
        tz <- ifelse(xconc[i] > z[j], -1 / kd[j] * log(1 - R), bigtime)
        tref <- max(tprec[i], tz)
        psurv[j] <- exp(-hb[j] * (time[i] - tprec[i]) +
                          if (time[i] > tz) {
                            -kk[j] * ((xconc[i] - z[j]) * (time[i] - tref) +
                                        xconc[i]/kd[j] * (exp(-kd[j] * time[i]) - exp(-kd[j] * tref)))
                          } else {
                            0
                          })
      }
      NsurvPred[, i] <- rbinom(niter, Nprec[i], psurv)
    }
    NsurvPred <- t(NsurvPred)
  }
  if(model_type == "IT"){
    
    D.max <- matrix(nrow = length(kd), ncol = length(xconc) )
    for(j in 1:n){
      # xconc[j] * (1-exp(-kd * time[j])) is always the max compared to previous time !
      D.max[,j] <- xconc[j] * (1-exp(-kd * time[j])) 
    }
    dtheo.IT <-  exp(-hb %*% t(time)) * (1 - plogis(log(D.max), location = log(alpha), scale = 1/beta))
    
    # transpose dtheo
    dtheo <- t(dtheo.IT)
    # dtheo <- do.call("rbind", lapply(dtheo.IT, t))
    
    i_prec = x$jags.data$i_prec
    ncol_NsurvPred = ncol(dtheo)
    
    NsurvPred = matrix(NA, ncol = ncol_NsurvPred, nrow = nrow(dtheo))
    
    for(i in 1:nrow(dtheo)){
      NsurvPred[i, ] = rbinom(ncol_NsurvPred, size = Nprec[i], prob = as.numeric(dtheo[i,] / dtheo[i_prec[i],]))
    }
    
  }
  
  QNsurvPred <- t(apply(NsurvPred, 1, quantile,
                        probs = c(2.5, 50, 97.5) / 100, na.rm = TRUE))
  
  tab <- data.frame(QNsurvPred,
                    Nprec,
                    NsurvObs,
                    col = ifelse(QNsurvPred[,"2.5%"] > NsurvObs |
                                   QNsurvPred[,"97.5%"] < NsurvObs,
                                 "red", "green"))
  colnames(tab) <- c("P2.5", "P50", "P97.5", "Nprec", "Obs", "col")
  
  return(tab)
}
