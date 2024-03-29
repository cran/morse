#' Predict the Multiplication Factor leading to x\% of reduction in survival
#' at a specific time.
#' 
#' Generic method for \code{MFx}, a function denoted \eqn{MF(x,t)} for 
#' \eqn{x}\% Multiplication Factor at time \eqn{t}.
#' 
#' When class of \code{object} is \code{survFit}, see \link[=MFx.survFit]{MFx.survFit}.
#' 
#' @param object An object used to select a method
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @return returns an object of class \code{MFx}
#' 
#' @export
#' 
MFx <- function(object, ...){
  UseMethod("MFx")
}

