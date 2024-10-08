#' Posterior predictive check plot for \code{reproFitTT} objects
#'
#' This is the generic \code{ppc} S3 method for the \code{reproFitTT} class.
#' It plots the predicted values with 95\% credible intervals versus the observed
#' values.
#' 
#' The coordinates of black points are the observed values of the cumulated number
#' of reproduction outputs for a given concentration (\eqn{X}-scale) and the corresponding 
#' predicted values (\eqn{Y}-scale). 95\% prediction intervals are added to each predicted
#' value, colored in green if this interval contains the observed value and in red
#' in the other case. As replicates are not pooled in this plot, overlapped points
#' are shifted on the \eqn{X-}axis to help the visualization of replicates. The bisecting
#' line (y = x) is added to the plot in order to see if each prediction interval
#' contains each observed value. As replicates are shifted on the \eqn{X}-axis, this
#' line may be represented by steps.
#'
#' @rdname PPC
#' 
#' @param x An object of class \code{reproFitTT}
#' @param xlab A label for the \eqn{X}-axis, by default \code{Observed Cumul. Nbr. of offspring}
#' @param ylab A label for the \eqn{Y}-axis, by default \code{Predicted Cumul. Nbr. of offspring}
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
ppc.reproFitTT <- function(x,
                           style = "ggplot", 
                           xlab = "Observed Cumul. Nbr. of offspring",
                           ylab = "Predicted Cumul. Nbr. of offspring",
                           main = NULL, ...) {
  if (!is(x, "reproFitTT"))
    stop("x is not of class 'reproFitTT'!")

  ppc_gen(EvalreproPpc(x), style, xlab, ylab, main)
}


#' @importFrom stats rgamma rpois quantile
EvalreproPpc <- function(x) {
  tot.mcmc <- do.call("rbind", x$mcmc)

  if (x$model.label == "GP") {
    omega <- 10^tot.mcmc[, "log10omega"]
  }
  b <- 10^tot.mcmc[, "log10b"]
  d <- tot.mcmc[, "d"]
  e <- 10^tot.mcmc[, "log10e"]

  niter <- nrow(tot.mcmc)
  n <- x$jags.data$n
  xconc <- x$jags.data$xconc
  Nindtime <- x$jags.data$Nindtime
  NcumulObs <- x$jags.data$Ncumul
  NcumulPred <- matrix(NA, nrow = niter, ncol = n)

  if (x$model.label == "GP") {
    for (i in 1:n) {
      popmean <- d / (1 + (xconc[i]/e)^b)
      indmean <- rgamma(n = niter, shape = popmean / omega, rate = 1 / omega)
      NcumulPred[, i] <- rpois(niter, indmean * Nindtime[i])
    }

  }
  if (x$model.label == "P") {
    for (i in 1:n) {
      ytheo <- d / (1 + (xconc[i]/e)^b)
      nbtheo <- ytheo * Nindtime[i]
      NcumulPred[, i] <- rpois(niter, nbtheo)
    }
  }
  QNreproPred <- t(apply(NcumulPred, 2, quantile,
                         probs = c(2.5, 50, 97.5) / 100))
  tab <- data.frame(QNreproPred,
                    Nindtime, NcumulObs,
                    col = ifelse(QNreproPred[,"2.5%"] > NcumulObs | QNreproPred[,"97.5%"] < NcumulObs,
                                 "red", "green"))
  colnames(tab) <- c("P2.5", "P50", "P97.5", "Nindtime", "Obs", "col")

  return(tab)
}
