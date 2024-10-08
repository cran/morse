#' Summary of \code{reproData} object
#' 
#' This is the generic \code{summary} S3 method for the \code{reproData} class.
#' It provides information about the structure of the data set and the experimental
#' design.
#' 
#' @param object an object of class \code{reproData}
#' @param quiet if \code{TRUE}, does not print
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @return The function returns a list with the same information than 
#' \code{\link{summary.survDataCstExp}} plus an additional one:
#' \item{NboffTimeConc}{nb of offspring for all concentrations and time points}
#' 
#' @keywords summary
#' 
#' @export
summary.reproData <- function(object, quiet = FALSE, ...) {
  res <- summary.survDataCstExp(object, quiet = quiet)
  
  # matrix of number of offspring (sum of all replicate) by time / conc
  ans3 <- tapply(object$Nrepro,
                 list(as.factor(object$conc), as.factor(object$time)), sum)
  
  if (! quiet) {
    cat("\nNumber of offspring (sum of replicates) per time and concentration: \n")
    print(ans3)
  }
  
  invisible(c(res, list(NboffTimeConc = ans3)))
}
