#' @title Create a list giving data to use in Bayesian inference.
#' 
#' @name modelData
#' 
#' @param x An object of class \code{survData}
#' @param model_type TKTD GUTS model type ('SD' or 'IT')
#' @param extend_time Number of for each replicate used for linear 
#' interpolation (comprise between time to compute and fitting accuracy)
#' @param \dots Further arguments to be passed to generic methods
#' 
#' 
#' @return A list for parameterization of priors for Bayesian inference.
#' 
modelData <- function(x, ...){
  UseMethod("modelData")
}
