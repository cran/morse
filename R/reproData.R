#' Creates a dataset for reproduction toxicity analysis
#'
#' This function creates a \code{reproData} object from experimental
#' data provided as a \code{data.frame}. The resulting object can then be used
#' for plotting and model fitting. The \code{reproData} class is a sub-class
#' of \code{survData}, meaning that all functions and method available for
#' survival analysis can be used with \code{reproData} objects.
#'
#' The \code{x} argument contains the experimental data, and should have
#' the same structure than the argument of \code{survData}, plus a single
#' additional column providing the total number of offspring observed since the
#' last time point. The function fails if \code{x} does not meet the
#' expected requirements. Please run \code{\link{reproDataCheck}} to ensure
#' \code{x} is well-formed.
#'
#' Note that experimental data with time-variable exposure are not supported.
#'
#' @aliases reproData
#'
#' @param x a dataframe as expected by \code{survData} containing one
#' additional \code{Nrepro} column of class \code{integer} with positive
#' values only. This column should
#' provide the number of offspring produced since the last observation.
#'
#' @return An object of class \code{reproData}.
#'
#' @keywords transformation
#'
#' @export
reproData <- function(x) {

  # test the integrity of the data with reproDataCheck
  if (dim(reproDataCheck(x, diagnosis.plot = FALSE))[1] > 0)
    stop("The [x] argument is not well-formed, please use [reproDataCheck] for details.")

  if (!is_exposure_constant(x)) {
    stop("[x] displays variable exposure, which is not supported for reproduction analysis")
  }

  x <- survData(x)
  x$Nindtime <- Nindtime(x)
  x$Nreprocumul <- Nreprocumul(x)

  # force concentration as type double
  x$conc <- as.double(x$conc)
  class(x) <- c("reproData", class(x))
  return(x)
}

# Computes cumulated offspring for each line of a
# \code{reproData} object
#
# @param x an object of class \code{reproData}
# @return an integer vector
#
# @importFrom dplyr left_join
#
# @export
#
Nreprocumul <- function(x) {
  T <- sort(unique(x$time)) # observation times
  Nreprocumul <- x$Nrepro
  for (i in 2:length(T)) {
    now <- x$time == T[i]
    before <- x$time == T[i - 1]
    Nreprocumul[now] <- Nreprocumul[before] + x$Nrepro[now]
  }
  return(Nreprocumul)
}
