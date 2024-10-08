#' Summary of \code{survDataCstExp} object
#' 
#' The generic \code{summary} S3 method for the \code{survDataCstExp} class provides
#' information about the structure of the data set and the experimental design.
#' 
#' @param object an object of class \code{survDataCstExp}
#' @param quiet when \code{TRUE}, does not print
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @return The function returns a list with the following information:
#' \item{NbrepTimeConc}{nb of replicates for all concentrations and time points}
#' \item{NbsurvTimeConc}{nb of survivors. for all concentrations and time points}
#' 
#' @keywords summary
#' 
#' @export
summary.survDataCstExp <- function(object, quiet = FALSE, ...) {
  # matrix of number of replicate by time / conc
  ans1 <- table(object[, c("conc", "time")])
  
  # matrix of number of survival (sum of all replicate) by time / conc
  ans2 <- tapply(object$Nsurv, list(as.factor(object$conc),
                                    as.factor(object$time)), sum)
  
  if (! quiet) {
    cat("\nNumber of replicates per time and concentration: \n")
    print(ans1)
    cat("\nNumber of survivors (sum of replicates) per time and concentration: \n")
    print(ans2)
  }
  
  invisible(list(NbrepTimeConc = ans1,
                 NbsurvTimeConc = ans2))
}
