#' Print of \code{survFitVarExp} object
#'
#' This is the generic \code{print} S3 method for the \code{survFitVarExp} class.
#' It prints the underlying JAGS model and some information on the Bayesian
#' inference procedure.
#'
#' @param x An object of class \code{survFitVarExp}
#' @param \dots Further arguments to be passed to generic methods.
#'
#' @keywords print
#' 
#' @return print the model text and the Jags Computing information
#'
#' @export
print.survFitVarExp <- function(x, ...) {
  # print the model text and the Jags Computing information
  # for an object of class survFitTKTD
  
  mcmcInfo = x$mcmcInfo
  
  # M.C.M.C. informations
  nbr.thin = mcmcInfo$nbr.thin
  mcmc_info =
  cat("Model:\n")
  print(x$model)
  cat("\nComputing information:\n\n")
  cat("Number of iterations per chain = ", mcmcInfo$n.iter, "\n")
  cat("Thinning interval =", mcmcInfo$thin.interval, "\n")
  cat("Number of chains =", mcmcInfo$n.chains, "\n")
  cat("Number iterations in warmup per chain =", mcmcInfo$n.warmup, "\n")
  cat("Sample size per chain =", mcmcInfo$n.iter / mcmcInfo$thin.interval , "\n")
}