#' Fits a Bayesian concentration-response model for target-time survival analysis
#' 
#' @param data an object used to select a method 'survFitTT'
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @return an object of class \code{survFitTT}

#' @export
survFitTT <- function(data, ...){
  UseMethod("survFitTT")
}
