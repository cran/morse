#' Print of \code{reproFitTT} object
#' 
#' This is the generic \code{print} S3 method for the \code{reproFitTT} class.
#' It prints the underlying JAGS model and some information on the Bayesian 
#' inference procedure.
#' 
#' @param x An object of class \code{reproFitTT}
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @return print the model text and the Jags Computing information
#' 
#' @keywords print
#' 
#' @export
print.reproFitTT <- function(x, ...) {
  # print the model text and the Jags Computing information
  # for an object of class reproFitTT
  
  # M.C.M.C. informations
  cat("Model:\n")
  print(x$model)
  cat("\nComputing information:\n\n")
  cat("\n", "Iterations = ", x$n.iter[["start"]], ":",
      x$n.iter[["end"]], "\n", sep = "")
  cat("Thinning interval =", x$n.thin, "\n")
  cat("Number of chains =", x$n.chains, "\n")
  cat("Sample size per chain =",
      (x$n.iter[["end"]] - x$n.iter[["start"]]) / x$n.thin + 1, "\n")
}
