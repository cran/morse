#' Plotting method for \code{reproFitTT} objects
#' 
#' This is the generic \code{plot} S3 method for the \code{reproFitTT} class.
#' It plots the concentration-effect fit under target time reproduction
#' analysis.
#' 
#' The fitted curve represents the \strong{estimated reproduction rate} at the target time
#'  as a function of the chemical compound concentration.
#' The function plots 95\% credible intervals for the estimated reproduction
#' rate (by default the grey area around the fitted curve). Typically
#' a good fit is expected to display a large overlap between the two types of intervals.
#' If spaghetti = TRUE, the credible intervals are represented by two dotted
#' lines limiting the credible band, and a spaghetti plot is added to this band.
#' It consists of the representation of simulated curves using parameter values
#' sampled in the posterior distribution (10\% of the MCMC chains are randomly
#' taken for this sample).
#'
#' @param x an object of class \code{reproFitTT}
#' @param xlab a label for the \eqn{X}-axis, by default \code{Concentration}
#' @param ylab a label for the \eqn{Y}-axis, by default \code{Nb of offspring per ind/day}
#' @param main main title for the plot
#' @param fitcol color of the fitted curve
#' @param fitlty line type of the fitted curve
#' @param fitlwd width of the fitted curve
#' @param spaghetti if \code{TRUE}, the credible interval is represented by 
#' multiple curves
#' @param cicol color of the 95 \% credible limits
#' @param cilty line type of the 95 \% credible limits
#' @param cilwd width of the 95 \% credible limits
#' @param ribcol color of the ribbon between lower and upper credible limits.
#' Transparent if \code{NULL}
#' @param addlegend if \code{TRUE}, adds a default legend to the plot
#' @param log.scale if \code{TRUE}, displays \eqn{X}-axis in log-scale
#' @param style graphical backend, can be \code{'ggplot'} or \code{'generic'}
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @note When \code{style = "generic"}, the function calls the generic function
#' \code{\link[graphics]{plot}}
#' @note When \code{style = "ggplot"}, the function return an object of class
#'  \code{ggplot}, see function \code{\link[ggplot2]{ggplot}} 
#'
#' 
#' @import ggplot2
#' @import grDevices
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom grid grid.rect gpar
#' @importFrom graphics plot axis legend lines par points polygon
#' segments title
#' @importFrom reshape2 melt
#'
#' @keywords plot
#' 
#' @return a plot of class \code{ggplot}
#' 
#' @export
plot.reproFitTT <- function(x,
                            xlab = "Concentration",
                            ylab = "Nb of offspring per ind/day",
                            main = NULL,
                            fitcol = "orange",
                            fitlty = 1,
                            fitlwd = 1,
                            spaghetti = FALSE,
                            cicol = "orange",
                            cilty = 2,
                            cilwd = 1,
                            ribcol = "grey70",
                            addlegend = FALSE,
                            log.scale = FALSE,
                            style = "ggplot", ...) {
  # plot the fitted curve estimated by reproFitTT
  # INPUTS
  # - x:  reproFitTT object
  # - xlab : label x
  # - ylab : label y
  # - main : main title
  # - fitcol : color fitted curve
  # - fitlty : type line fitted curve
  # - fitlwd : width line fitted curve
  # - cicol : color ci
  # - cilty : type line ci
  # - cilwd : width line ci
  # - addlegend : boolean
  # - log.scale : x log option
  # - style : generic ou ggplot
  # OUTPUT:
  # - plot of fitted regression
  
  # Selection of datapoints that can be displayed given the type of scale
  sel <- if (log.scale) x$dataTT$conc > 0 else TRUE
  
  dataTT <- x$dataTT[sel, ]
  dataTT$resp <- dataTT$Nreprocumul / dataTT$Nindtime
  transf_data_conc <- optLogTransform(log.scale, dataTT$conc)
  
  # Concentration values used for display in linear scale
  display.conc <- (function() {
    x <- optLogTransform(log.scale, dataTT$conc)
    s <- seq(min(x),max(x), length = 100)
    if(log.scale) exp(s) else s
  })()
  
  # Possibly log transformed concentration values for display
  curv_conc <- optLogTransform(log.scale, display.conc)
  
  cred.int <- reproMeanCredInt(x, display.conc)
  
  spaghetti.CI <- if (spaghetti) { reproSpaghetti(x, display.conc) } else NULL
  dataCIm <- if (spaghetti) { melt(cbind(curv_conc, spaghetti.CI),
                                   id.vars = c("curv_conc", "conc"))} else NULL
  
  curv_resp <- data.frame(conc = curv_conc, resp = cred.int[["q50"]],
                          Line = "loglogistic")
  
  # ylim
  ylim_CI <- if (spaghetti) { max(dataCIm$value, cred.int$qsup95)
  } else {
    max(cred.int$qsup95)
  }

  if (style == "generic") {
    reproFitPlotGenericCredInt(x, dataTT$conc, transf_data_conc, dataTT$resp,
                               curv_conc, curv_resp,
                               cred.int, spaghetti.CI, dataCIm,
                               xlab, ylab, fitcol, fitlty, fitlwd,
                               main, addlegend,
                               cicol, cilty, cilwd, ribcol, log.scale, ylim_CI)
  }
  else if (style == "ggplot") {
    reproFitPlotGG(x, dataTT$conc, transf_data_conc, dataTT$resp,
                   curv_conc, curv_resp,
                   cred.int, spaghetti.CI, dataCIm,
                   xlab, ylab, fitcol, fitlty, fitlwd,
                   main, addlegend,
                   cicol, cilty, cilwd, ribcol, log.scale, ylim_CI)
  }
  else stop("Unknown style")
}

#' @importFrom stats quantile rgamma
reproMeanCredInt <- function(fit, x) {
  # create the parameters for credible interval for the log logistic model
  # INPUT:
  # - fit : object of class reproFitTT
  # - x : vector of concentrations values (x axis)
  # OUTPUT:
  # - ci : credible limit
  
  mctot <- do.call("rbind", fit$mcmc)
  k <- nrow(mctot)
  # parameters
  d2 <- mctot[, "d"]
  log10b2 <- mctot[, "log10b"]
  b2 <- 10^log10b2
  log10e2 <- mctot[, "log10e"]
  e2 <- 10^log10e2
  
  # quantiles
  qinf95 = NULL
  q50 = NULL
  qsup95 = NULL
  
  # poisson
  if (fit$model.label == "P") {
    for (i in 1:length(x)) {
      theomean <- d2 / (1 + (x[i] / e2)^(b2)) # mean curve
      # IC 95%
      qinf95[i] <- quantile(theomean, probs = 0.025, na.rm = TRUE)
      qsup95[i] <- quantile(theomean, probs = 0.975, na.rm = TRUE)
      q50[i] <- quantile(theomean, probs = 0.5, na.rm = TRUE)
    }
  }
  
  # gamma poisson
  else if (fit$model.label == "GP") {
    # parameters
    log10omega2 <- mctot[, "log10omega"]
    omega2 <- 10^(log10omega2)
    
    for (i in 1:length(x)) {
      theomean <- d2 / (1 + (x[i] / e2)^(b2)) # mean curve
      theo <- rgamma(n = k, shape = theomean / omega2, rate = 1 / omega2)
      # IC 95%
      qinf95[i] <- quantile(theo, probs = 0.025, na.rm = TRUE)
      qsup95[i] <- quantile(theo, probs = 0.975, na.rm = TRUE)
      q50[i] <- quantile(theo, probs = 0.5, na.rm = TRUE)
    }
  }
  # values for cred.int
  ci <- data.frame(qinf95 = qinf95,
                   q50 = q50,
                   qsup95 = qsup95)
  
  return(ci)
}

reproSpaghetti <- function(fit, x) {
  mctot <- do.call("rbind", fit$mcmc)
  sel <- sample(nrow(mctot))[1:ceiling(nrow(mctot) / 10)]
  k <- nrow(mctot[sel,])
  # parameters
  d2 <- mctot[, "d"][sel]
  log10b2 <- mctot[, "log10b"][sel]
  b2 <- 10^log10b2
  log10e2 <- mctot[, "log10e"][sel]
  e2 <- 10^log10e2
  if (fit$model.label == "GP") {
    log10omega2 <- mctot[, "log10omega"][sel]
    omega2 <- 10^(log10omega2)
  }
  
  # all theorical
  dtheo <- array(data = NA, dim = c(length(x), length(e2)))
  if (fit$model.label == "GP") dtheotemp <- dtheo
  for (i in 1:length(e2)) {
    if (fit$model.label == "P") {
      dtheo[, i] <- d2[i] / (1 + (x / e2[i])^(b2[i])) # mean curve
    }
    else if (fit$model.label == "GP") {
      dtheotemp[, i] <- d2[i] / (1 + (x / e2[i])^(b2[i])) # mean curve
      dtheo[, i] <- rgamma(n = length(x), shape = dtheotemp[, i] / omega2[i],
                           rate = 1 / omega2[i])
    }
  }
  dtheof <- as.data.frame(cbind(x, dtheo))
  names(dtheof) <- c("conc", paste0("X", 1:length(sel)))
  
  return(dtheof)
}

#' @importFrom epitools pois.exact
reproFitPlotGenericCredInt <- function(x, data_conc, transf_data_conc, data_resp,
                                       curv_conc, curv_resp,
                                       cred.int, spaghetti.CI, dataCIm,
                                       xlab, ylab, fitcol, fitlty, fitlwd,
                                       main, addlegend,
                                       cicol, cilty, cilwd, ribcol, log.scale, ylim_CI) {

  # plot the fitted curve estimated by reproFitTT
  # with generic style with credible interval
  
  plot(transf_data_conc, data_resp,
       xlab = xlab,
       ylab = ylab,
       main = main,
       xaxt = "n",
       yaxt = "n",
       ylim = c(0, ylim_CI + 0.01),
       type = "n")
  
  # axis
  axis(side = 2, at = pretty(c(0, ylim_CI)))
  axis(side = 1,
       at = transf_data_conc,
       labels = data_conc)
  
  # Plotting the theoretical curve
  # cred.int ribbon + lines
  if (!is.null(spaghetti.CI)) {
    color <- "gray"
    color_transparent <- adjustcolor(color, alpha.f = 0.05)
    by(dataCIm, dataCIm$variable, function(x) {
      lines(x$curv_conc, x$value, col = color_transparent)
    })
  } else {
    polygon(c(curv_conc, rev(curv_conc)), c(cred.int[["qinf95"]],
                                            rev(cred.int[["qsup95"]])),
            col = ribcol, border = NA)
  }
  
  lines(curv_conc, cred.int[["qsup95"]], type = "l", col = cicol, lty = cilty,
        lwd = cilwd)
  lines(curv_conc, cred.int[["qinf95"]], type = "l", col = cicol, lty = cilty,
        lwd = cilwd)
  
  # fitted curve
  lines(curv_conc, curv_resp[, "resp"], col = fitcol,
        lty = fitlty, lwd = fitlwd, type = "l")
  
  # legend
  if(addlegend)  {
    legend("bottomleft",
           lty = c(cilty, fitlty),
           lwd = c(cilwd, fitlwd),
           col = c(cicol, fitcol),
           legend = c("Credible limits", "loglogistic"),
           bty = "n")
  }
}

reproFitPlotGGCredInt <- function(curv_resp, cred.int, spaghetti.CI, dataCIm,
                                  cicol, cilty, cilwd, valCols, fitlty, fitlwd, ribcol,
                                  xlab, ylab, main, ylim_CI) {
  # IC
  data.three <- data.frame(conc = curv_resp$conc,
                           qinf95 = cred.int[["qinf95"]],
                           qsup95 = cred.int[["qsup95"]],
                           Cred.Lim = "Credible limits")
  
  plt_31 <- if (!is.null(spaghetti.CI)) {
    ggplot(data.three) + geom_line(data = dataCIm, aes(x = curv_conc, y = value,
                                                       group = variable),
                                   col = "gray", alpha = 0.05)
  } else {
    ggplot(data.three) + geom_ribbon(data = data.three, aes(x = conc,
                                                            ymin = qinf95,
                                                            ymax = qsup95),
                                     fill = ribcol, col = NA,
                                     alpha = 0.4)
  }
  
  plt_3 <- plt_31 +
    geom_line(data = data.three, aes(conc, qinf95, color = Cred.Lim),
              linetype = cilty, size = cilwd) +
    geom_line(data = data.three, aes(conc, qsup95, color = Cred.Lim),
              linetype = cilty, size = cilwd) +
    scale_color_manual("", values = valCols$cols4) +
    theme_minimal()
  
  # plot IC
  # final plot
  
  if (!is.null(spaghetti.CI)) {
    plt_40 <- ggplot(data.three) +
      geom_line(data = dataCIm, aes(x = curv_conc, y = value, group = variable),
                col = "gray", alpha = 0.05)
  } else {
    plt_40 <- ggplot(data.three) + geom_ribbon(data = data.three,
                                               aes(x = conc,
                                                   ymin = qinf95,
                                                   ymax = qsup95),
                                               fill = ribcol,
                                               col = NA, alpha = 0.4)
  }
  
  plt_4 <- plt_40 +
    geom_line(data = data.three, aes(conc, qinf95),
              linetype = cilty, size = cilwd, color = valCols$cols4) +
    geom_line(data = data.three, aes(conc, qsup95),
              linetype = cilty, size = cilwd, color = valCols$cols4) +
    geom_line(aes(conc, resp), curv_resp,
              linetype = fitlty, size = fitlwd, color = valCols$cols2) +
    ylim(0, ylim_CI + 0.2) +
    labs(x = xlab, y = ylab) +
    ggtitle(main) + theme_minimal()
  
  return(list(plt_3 = plt_3,
              plt_4 = plt_4))
}

reproFitPlotGG <- function(x, data_conc, transf_data_conc, data_resp,
                           curv_conc, curv_resp,
                           cred.int, spaghetti.CI, dataCIm,
                           xlab, ylab, fitcol, fitlty, fitlwd,
                           main, addlegend,
                           cicol, cilty, cilwd, ribcol, log.scale, ylim_CI) {
  
  if (Sys.getenv("RSTUDIO") == "") {
    dev.new() # create a new page plot
    # when not use RStudio
  }
  
  # dataframes points (data) and curve (curv)
  # colors
  valCols <- fCols(curv_resp, fitcol, cicol)

  plt_4 <-
    reproFitPlotGGCredInt(curv_resp, cred.int, spaghetti.CI, dataCIm,
                          cicol, cilty, cilwd, valCols, fitlty, fitlwd, ribcol, xlab,
                          ylab, main, ylim_CI)$plt_4
  
  if (addlegend) {
    
    # create legends
    
    # curve (to create the legend)
    plt_2 <- ggplot(curv_resp) +
      geom_line(data = curv_resp, aes(conc, resp, colour = Line),
                linetype = fitlty, size = fitlwd) +
      scale_color_manual("", values = valCols$cols2) +
      theme_minimal()
    
    mylegend_2 <- legendGgplotFit(plt_2) # mean line legend
    
    plt_5 <- plt_4 + scale_x_continuous(breaks = transf_data_conc,
                                        labels = data_conc)
    
    plt_3 <- reproFitPlotGGCredInt(curv_resp, cred.int, spaghetti.CI, dataCIm,
                                   cicol, cilty, cilwd, valCols, fitlty,
                                   fitlwd, ribcol, xlab, ylab, main, ylim_CI)$plt_3
    
    mylegend_3 <- legendGgplotFit(plt_3)
    
    grid.arrange(plt_5, arrangeGrob(mylegend_2, mylegend_3,
                                    nrow = 6), ncol = 2,
                 widths = c(6, 2))
  }
  else { # no legend
    plt_5 <- plt_4 + scale_x_continuous(breaks = transf_data_conc,
                                        labels = data_conc)
    return(plt_5)
  }
}

