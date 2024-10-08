#' Summary of \code{reproFitTT} object
#'
#' This is the generic \code{summary} S3 method for the \code{reproFitTT} class.
#' It shows the quantiles of priors and posteriors on parameters
#' and the quantiles of the posterior on the ECx estimates.
#'
#' @param object an object of class \code{reproFitTT}
#' @param quiet when \code{TRUE}, does not print
#' @param \dots Further arguments to be passed to generic methods
#'
#' @return The function returns a list with the following information:
#' \item{Qpriors}{quantiles of the model priors}
#' \item{Qposteriors}{quantiles of the model posteriors}
#' \item{QECx}{quantiles of ECx estimates}
#'
#' @keywords summary
#' 
#' @importFrom stats qnorm qunif
#' 
#' @export
summary.reproFitTT <- function(object, quiet = FALSE, ...) {

  # quantiles of priors parameters
  n.iter <- object$n.iter$end - object$n.iter$start

  # b
  log10b <- qunif(p = c(0.5, 0.025, 0.975),
                  min = object$jags.data$log10bmin,
                  max = object$jags.data$log10bmax)

  b <- 10^log10b

  # d
  d <- qnorm(p = c(0.5, 0.025, 0.975),
             mean = object$jags.data$meand,
             sd = 1 / sqrt(object$jags.data$taud))

  # e
  log10e <- qnorm(p = c(0.5, 0.025, 0.975),
                  mean = object$jags.data$meanlog10e,
                  sd = 1 / sqrt(object$jags.data$taulog10e))

  e <- 10^log10e

  if (object$model.label == "P") {
    res <- rbind(b, d, e)
  }
  if (object$model.label == "GP") {
    # omega
    log10omega <- qunif(p = c(0.5, 0.025, 0.975),
                        min = object$jags.data$log10omegamin,
                        max = object$jags.data$log10omegamax)

    omega <- 10^log10omega

    res <- rbind(b, d, e, omega)
  }

  ans1 <-  format(data.frame(res), scientific = TRUE, digits = 4)
  colnames(ans1) <- c("50%", "2.5%", "97.5%")

  # quantiles of estimated model parameters
  ans2 <- format(object$estim.par, scientific = TRUE, digits = 4)
  colnames(ans2) <- c("50%", "2.5%", "97.5%")

  # estimated ECx and their CIs 95%
  ans3 <- format(object$estim.ECx, scientific = TRUE, digits = 4)
  colnames(ans3) <- c("50%", "2.5%", "97.5%")

  if (! quiet) {
    cat("Summary: \n\n")
    if (object$model.label == "GP")
      cat("The ", object$det.part, " model with a Gamma Poisson stochastic part was used !\n\n")
    if(object$model.label == "P")
      cat("The ", object$det.part, " model with a Poisson stochastic part was used !\n\n")
    cat("Priors on parameters (quantiles):\n\n")
    print(ans1)
    cat("\nPosteriors of the parameters (quantiles):\n\n")
    print(ans2)
    cat("\nPosteriors of the ECx (quantiles):\n\n")
    print(ans3)
  }

  invisible(list(Qpriors = ans1,
                 Qposteriors = ans2,
                 QECx = ans3))
}

