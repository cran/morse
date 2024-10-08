#' Plotting method for \code{survFit} objects
#'
#' This is the generic \code{plot} S3 method for the
#' \code{survFit}.  It plots the fit obtained for each
#' concentration of chemical compound in the original dataset.
#'
#' The fitted curves represent the \strong{estimated survival probability} as a function
#' of time for each concentration.
#' The black dots depict the \strong{observed survival
#' probability} at each time point. Note that since our model does not take
#' inter-replicate variability into consideration, replicates are systematically
#' pooled in this plot.
#' The function plots both 95\% credible intervals for the estimated survival
#' probability (by default the grey area around the fitted curve) and 95\%  binomial confidence
#' intervals for the observed survival probability (as black error bars if
#' \code{adddata = TRUE}).
#' Both types of intervals are taken at the same level. Typically
#' a good fit is expected to display a large overlap between the two types of intervals.
#' If \code{spaghetti = TRUE}, the credible intervals are represented by two
#' dotted lines limiting the credible band, and a spaghetti plot is added to this band.
#' This spaghetti plot consists of the representation of simulated curves using parameter values
#' sampled in the posterior distribution (2\% of the MCMC chains are randomly
#' taken for this sample).
#'
#' @param x An object of class \code{survFit}.
#' @param xlab A label for the \eqn{X}-axis, by default \code{Time}.
#' @param ylab A label for the \eqn{Y}-axis, by default \code{Survival probability}.
#' @param main A main title for the plot.
#' @param concentration A numeric value corresponding to some specific concentrations in
#' \code{data}. If \code{concentration = NULL}, draws a plot for each concentration.
#' @param spaghetti if \code{TRUE}, draws a set of survival curves using
#' parameters drawn from the posterior distribution
#' @param one.plot if \code{TRUE}, draws all the estimated curves in
#' one plot instead of one plot per concentration.
#' @param adddata if \code{TRUE}, adds the observed data to the plot
#' with (frequentist binomial) confidence intervals
#' @param addlegend if \code{TRUE}, adds a default legend to the plot.
#' @param style graphical backend, can be \code{'generic'} or \code{'ggplot'}
#' @param \dots Further arguments to be passed to generic methods.
#'
#' @keywords plot
#' 
#' @return a plot of class \code{ggplot}
#'
#' @export
#'
#' @import ggplot2
#' @import grDevices
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr contains
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom grid grid.rect gpar
#' @importFrom graphics plot
#' @importFrom tidyr tibble
#' @importFrom tidyr as_tibble
#' 
plot.survFitCstExp <- function(x,
                         xlab = "Time",
                         ylab = "Survival probability",
                         main = NULL,
                         concentration = NULL,
                         spaghetti = FALSE,
                         one.plot = FALSE,
                         adddata = TRUE,
                         addlegend = FALSE,
                         style = "ggplot", ...) {
  
  if (one.plot && !is.null(concentration)) one.plot <- FALSE
  
  if ((addlegend && is.null(concentration)) ||
      (addlegend && !one.plot))
    warning("The legend is available only if [one.plot] is TRUE or if [concentration] is not NULL !", call. = FALSE)
  
  if (!is.null(concentration) && !any(x$transformed.data$conc == concentration))
    stop("The [concentration] argument is not one of the possible concentration !")
  
  if (one.plot)
    warning("The credible limits and confidence intervals are not drawn when 'one.plot' = TRUE.", call. = FALSE)
  
  conf.int <- survTKTDConfInt_CstExp(x)
  
  dobs <- data.frame(conc = x$transformed.data$conc,
                     time = x$transformed.data$time, 
                     psurv = x$transformed.data$Nsurv / x$transformed.data$Ninit,
                     Points = "Observed values",
                     #color = as.numeric(as.factor(x$transformed.data$conc)),
                     color = x$transformed.data$replicate,
                     conf.int)
  
  
  # remove time 0 in dobs
  dobs <- dplyr::filter(dobs, time != 0)
  
  data.credInt <- survFitPlotCITKTD_CstExp(x)
  
  if (style == "generic") {
    survFitPlotTKTDGeneric(data.credInt, dobs, xlab, ylab, main, concentration,
                           one.plot, spaghetti,
                           adddata, addlegend)
  } else if (style == "ggplot") {
    survFitPlotTKTDGG(data.credInt, dobs, xlab, ylab, main, concentration,
                      one.plot, spaghetti, adddata, addlegend)
  } else stop("Unknown style")
}

#' @importFrom stats aggregate binom.test
survTKTDConfInt_CstExp <- function(x) {
  # create confidente interval on observed data
  # binomial model by a binomial test
  # INPUT:
  # - x : object of class survFitTT
  # OUTPUT:
  # - ci : confidente interval

  ci <- apply(x$transformed.data, 1, function(x) {
    binom.test(as.numeric(x["Nsurv"]), as.numeric(x["Ninit"]))$conf.int
  })
  ci <- as.data.frame(t(ci))
  colnames(ci) <- c("qinf95", "qsup95")
  ci$Conf.Int <- "Confidence interval"
  
  return(ci)
}

Surv_SD <- function(Cw, time, kk, kd, z, hb) {
  S <- exp(-hb*time)
  x <- ifelse(Cw > z, 1 - z/Cw, NA)
  tz <- -(1/kd)*log(x)
  y <- ifelse(time > tz,
              exp( kk/kd*Cw*(exp(-kd*tz) - exp(-kd*time))
                   - kk*(Cw - z)*(time - tz)),
              NA)
  return(ifelse(!is.na(x) & !is.na(y), S * y, S))
}

Surv_IT <- function(Cw, time, kd, hb, alpha, beta){
  D <- Cw*(1 - exp(-kd %*% t(time)))
  D.max <- t(apply(D, 1, cummax))
  S <- exp(-hb %*% t(time)) * (1 - plogis(log(D.max), location = log(alpha), scale = 1/beta))
  return(S)
}

survFitPlotCITKTD_CstExp <- function(x) {
  # INPUT
  # x : An object of class survFitTKTD
  # OUTPUT
  # A list of - dobs : observed values
  #           - dtheo : estimated values
  npoints <- 100
  
  concobs <- unique(x$transformed.data$conc)
  ## tfin <- seq(0, max(x$jags.data$t), length.out = npoints)
  tfin <- seq(0, max(x$jags.data$time), length.out = npoints)
  
  # prameters
  mctot <- do.call("rbind", x$mcmc)
  kd <- 10^mctot[, "kd_log10"]
  # "hb" is not in survFit object of morse <v3.2.0
  if("hb" %in% colnames(mctot)){
    hb <- mctot[, "hb"]
  } else{ hb <- 10^mctot[, "hb_log10"] }
  
  # all theoretical
  k <- 1:length(concobs)
  j <- 1:npoints
  
  model_type = x$model_type
  if(model_type == "SD"){
    z <- 10^mctot[, "z_log10"]
    kk <- 10^mctot[, "kk_log10"]
    
    dtheo <- lapply(k, function(kit) { # conc
      sapply(j, function(jit) { # time
        Surv_SD(Cw = concobs[kit],
                time = tfin[jit],
                kk = kk,
                kd = kd,
                z = z,
                hb = hb)
      })
    })
  }
  if(model_type == "IT"){
    alpha <- 10^mctot[, "alpha_log10"]
    beta <- 10^mctot[, "beta_log10"]
    
    dtheo <- lapply(k, function(kit) { # concentration pour chaque concentration
      Surv_IT(Cw = concobs[kit],
              time = tfin,
              kd = kd,
              hb = hb,
              alpha = alpha,
              beta = beta)
    })
  }

  # transpose dtheo
  dtheo <- do.call("rbind", lapply(dtheo, t))
  
  # quantile
  qinf95 <- apply(dtheo, 1, quantile, probs = 0.025, na.rm = TRUE)
  qsup95 <- apply(dtheo, 1, quantile, probs = 0.975, na.rm = TRUE)
  q50 <- apply(dtheo, 1, quantile, probs = 0.5, na.rm = TRUE)
  
  
  dtheo <- as.data.frame(dtheo)
  colnames(dtheo) <- paste0("X", 1:length(kd))
  # divide number of mcmc by 50
  sel <- sample(ncol(dtheo))[1:ceiling(ncol(dtheo) / 50)]
  dtheo <- dtheo[, sel]
  
  dtheo$conc <- rep(concobs, rep(npoints, length(concobs)))
  dtheo$time <- rep(tfin, length(concobs))
  
  # add credible limits
  dtheo$qinf95 <- qinf95
  dtheo$qsup95 <- qsup95
  dtheo$q50 <- q50
  
  # names for legend plots
  dtheo$Cred.Lim <- "Credible limits"
  dtheo$Mean.C <- "Mean curve"
  # vector color
  dtheo$color <- as.numeric(as.factor(dtheo$conc))


  return(dtheo)
}