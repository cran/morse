#' Summary of \code{survFitTT} object
#'
#' This is the generic \code{summary} S3 method for the \code{survFitTT} class.
#' It shows the quantiles of priors and posteriors on parameters and the quantiles
#' of the posteriors on the LCx estimates.
#'
#' @param object an object of class \code{survFitTT}
#' @param quiet when \code{TRUE}, does not print
#' @param \dots Further arguments to be passed to generic methods
#'
#' @return The function returns a list with the following information:
#' \item{Qpriors}{quantiles of the model priors}
#' \item{Qposteriors}{quantiles of the model posteriors}
#' \item{QLCx}{quantiles of LCx estimates}
#'
#'
#' @keywords summary
#'
#' @importFrom stats qnorm qunif
#' 
#' @export
summary.survFitTT <- function(object, quiet = FALSE, ...) {

  # quantiles of priors parameters
  n.iter <- object$n.iter$end - object$n.iter$start

  # b
  log10b <- qunif(p = c(0.5, 0.025, 0.975),
                  min = object$jags.data$log10bmin,
                  max = object$jags.data$log10bmax)

  b <- 10^log10b

  # e
  log10e <- qnorm(p = c(0.5, 0.025, 0.975),
                  mean = object$jags.data$meanlog10e,
                  sd = 1 / sqrt(object$jags.data$taulog10e))

  e <- 10^log10e

  # d
  if (object$det.part == "loglogisticbinom_3") {

    d <- qunif(p = c(0.5, 0.025, 0.975),
               min = object$jags.data$dmin,
               max = object$jags.data$dmax)

    res <- rbind(b, d, e)
  } else {
    res <- rbind(b, e)
  }

  ans1 <- format(data.frame(res), scientific = TRUE, digits = 4)
  colnames(ans1) <- c("50%", "2.5%", "97.5%")

  # quantiles of estimated model parameters
  ans2 <- format(object$estim.par, scientific = TRUE, digits = 4)
  colnames(ans2) <- c("50%", "2.5%", "97.5%")

  # estimated ECx and their CIs 95%
  ans3 <- format(object$estim.LCx, scientific = TRUE, digits = 4)
  colnames(ans3) <- c("50%", "2.5%", "97.5%")

  # print
  if (! quiet) {
    cat("Summary: \n\n")
    cat("The ", object$det.part, " model with a binomial stochastic part was used !\n\n")
    cat("Priors on parameters (quantiles):\n\n")
    print(ans1)
    cat("\nPosteriors of the parameters (quantiles):\n\n")
    print(ans2)
    cat("\nPosteriors of the LCx (quantiles):\n\n")
    print(ans3)
  }

  invisible(list(Qpriors = ans1,
                 Qpost = ans2,
                 QLCx = ans3))
}

