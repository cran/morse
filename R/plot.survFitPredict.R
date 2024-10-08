#' Plotting method for \code{survFitPredict} objects
#'
#' This is the generic \code{plot} S3 method for the
#' \code{survFitPredict}.  It plots the predicted survival probability for each
#' concentration of the chemical compound in the provided dataset.
#'
#' The fitted curves represent the \strong{predicted survival probability} as a function
#' of time for each concentration.
#' The function plots both the 95\% credible band and the predicted survival
#' probability over time.
#' If \code{spaghetti = TRUE}, the credible intervals are represented by two
#' dotted lines limiting the credible band, and a spaghetti plot is added to this band.
#' This spaghetti plot consists of the representation of simulated curves using parameter values
#' sampled in the posterior distribution (10\% of the MCMC chains are randomly
#' taken for this sample).
#'
#' @param x An object of class \code{survFitPredict}.
#' @param xlab A label for the \eqn{X}-axis, by default \code{Time}.
#' @param ylab A label for the \eqn{Y}-axis, by default \code{Survival probability}.
#' @param main A main title for the plot.
#' @param spaghetti If \code{TRUE}, draws a set of survival curves using
#' parameters drawn from the posterior distribution
#' @param one.plot if \code{TRUE}, draws all the estimated curves in
#' one plot instead of one plot per concentration.
#' @param mcmc_size A numerical value refering by default to the size of the mcmc in object \code{survFitPredict}.
#'  This option is specific to \code{survFitPredict} objects for which computing time may be long.
#'  \code{mcmc_size} can be used to reduce the number of mcmc samples in order to speed up
#'  the computation.
#'  
#' @param \dots Further arguments to be passed to generic methods.
#'
#' @keywords plot
#' 
#' @return a plot of class \code{ggplot}
#' 
#' @export
#' 
#' @importFrom tidyr gather
#'
plot.survFitPredict <- function(x,
                               xlab = "Time",
                               ylab = "Survival probability",
                               main = NULL,
                               spaghetti = FALSE,
                               one.plot = FALSE,
                               mcmc_size = NULL,
                               ...) {

  df_prediction <-  x$df_quantile
  df_spaghetti <-  x$df_spaghetti

  # Plot
  plt <- ggplot() +
    theme_minimal() +
    scale_x_continuous(name = xlab) +
    scale_y_continuous(name = ylab,
                       limits = c(0,1)) +
    theme(legend.position = "top")
  
  # spaghetti
  if(spaghetti == TRUE){
    
    df_spaghetti_gather <- df_spaghetti %>%
      tidyr::gather(survRate_key, survRate_value, -c(time,conc,replicate))
    
    plt <- plt +
      geom_line(data = df_spaghetti_gather,
                aes(x = time, y = survRate_value, group = interaction(survRate_key, replicate)),
                alpha = 0.2) +
      geom_line(data = df_prediction,
                aes(x = time, y= qinf95, group = replicate),
                color = "orange", linetype = 2) +
      geom_line(data = df_prediction,
                aes(x = time, y = qsup95, group = replicate),
                color = "orange", linetype = 2)
  }
  if(spaghetti != TRUE){
    plt <- plt + 
      geom_ribbon(data = df_prediction,
                  aes(x = time, ymin = qinf95,ymax = qsup95, group = replicate),
                  fill = "grey70", alpha = 0.4)
  }
  # Prediction
  plt <- plt +
    geom_line(data = df_prediction,
              aes(x = time, y = q50, group = replicate),
              col="orange", size = 1)

  # facetting
  if(one.plot == FALSE){
    plt <- plt + facet_wrap(~ replicate)
  }  
  
  return(plt)
}

