#' @title Plotting method for \code{LCx} objects
#'
#' @description
#' This is the generic \code{plot} S3 method for the LCx class. 
#' It plots the survival probability as a function of concentration.
#'
#'
#' @param x An object of class \code{LCx}.
#' @param xlab A label for the \eqn{X}-axis, by default \code{Concentration}.
#' @param ylab A label for the \eqn{Y}-axis, by default \code{Survival probability median and 95 CI}.
#' @param main A main title for the plot.
#' @param subtitle A subtitle for the plot
#' @param \dots Further arguments to be passed to generic methods.
#'
#' @keywords plot
#' 
#' @return a plot of class \code{ggplot}
#'
#' @export
#'
plot.LCx <- function(x,
                     xlab = "Concentration",
                     ylab = "Survival probability \n median and 95 CI",
                     main = NULL,
                     subtitle = NULL,
                     ...){
  
  df_dose <- x$df_dose
  df_LCx <- x$df_LCx
  X_prop <- x$X_prop
  X_prop_provided <- x$X_prop_provided
  time_LCx <- x$time_LCx
  
  # Check if df_LCx are not all NA:
  if (all(is.na(df_LCx$LCx))) {
    warning(paste0("No LCx can be computed at X=", 100-X_prop_provided*100, " and time_LCx=", time_LCx, 
                   "\nSee the grey dotted line on the graph does not cross zone of credible interval.",
                   "\nUse LCx function 'LCx' with other values for arguments 'time_LCx' (default is the maximum time of the experimental data),
and/or other 'X'."))
  } else{
    legend.point = data.frame(
      x.pts = df_LCx$LCx,
      y.pts = rep(X_prop,3),
      pts.leg = c(paste("median: ", round(df_LCx$LCx[1],digits = 2)),
                  paste("quantile 2.5%: ", round(df_LCx$LCx[2],digits = 2)),
                  paste("quantile 97.5%: ", round(df_LCx$LCx[3],digits = 2)))
    )
  }
  
  
  if(is.null(main)){
    main <- paste("Concentration-response curve: LC", 100 - X_prop_provided*100, " at time", time_LCx)
  } 
  
  LCx_plt <- ggplot() + theme_minimal() +
    theme(legend.position = "top",
          legend.title = element_blank())+
    scale_y_continuous(limits = c(0,1)) +
    labs(title = main,
         subtitle = subtitle,
         x = xlab,
         y = ylab) +
    geom_ribbon(data = df_dose,
                aes(x = concentration, ymin = qinf95, ymax = qsup95), fill = "lightgrey") + 
    geom_line(data = df_dose,
              aes( x = concentration, y = q50), color ="orange") +
    geom_hline(yintercept = X_prop, col="grey70", linetype=2)
  
    # LCx points
    if(!all(is.na(df_LCx$LCx))){
      LCx_plt <- LCx_plt  +
        geom_point(data = legend.point,
                   aes(x = x.pts, y=y.pts,
                       color= pts.leg))+
        
        scale_color_manual(values=c("orange", "black", "black"))
    }

  return(LCx_plt)
  
}
