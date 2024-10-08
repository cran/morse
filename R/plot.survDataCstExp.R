#' Plotting method for \code{survData} objects
#'
#' This is the generic \code{plot} S3 method for the \code{survData} class.
#' It plots the number of survivors as a function of time.
#'
#' @param x an object of class \code{survData}
#' @param xlab a label for the \eqn{X}-axis, by default \code{Time}
#' @param ylab a label for the \eqn{Y}-axis, by default \code{Number of survivors}
#' @param main main title for the plot
#' @param concentration a numeric value corresponding to some concentration(s) in
#' \code{data}. If \code{concentration = NULL}, draws a plot for each concentration
#' @param style graphical backend, can be \code{'generic'} or \code{'ggplot'}
#' @param pool.replicate if \code{TRUE}, the datapoints of each replicate are
#' summed for a same concentration
#' @param addlegend if \code{TRUE}, adds a default legend to the plot
#' @param remove.someLabels if \code{TRUE}, removes 3/4 of \eqn{X}-axis labels in
#' \code{'ggplot'} style to avoid label overlap
#' @param \dots Further arguments to be passed to generic methods
#'
#' @note When \code{style = "ggplot"} (default), the function calls function
#' \code{\link[ggplot2]{ggplot}} and returns an object of class \code{ggplot}.
#'
#' @keywords plot
#' 
#' @return a plot of class \code{ggplot}
#'
#' @import ggplot2
#' @import grDevices
#' @import dplyr
#' @importFrom graphics plot axis legend lines par points title
#' @importFrom methods is
#' @importFrom stats aggregate
#'
#' @export
plot.survDataCstExp <- function(x,
                                xlab = "Time",
                                ylab = "Number of survivors",
                                main = NULL,
                                concentration = NULL,
                                style = "ggplot",
                                pool.replicate = FALSE,
                                addlegend = FALSE,
                                remove.someLabels = FALSE, ...) {

  if (!is(x,"survDataCstExp"))
    stop("plot.survData: object of class survData expected")

  if (style == "generic" && remove.someLabels)
    warning("'remove.someLabels' argument is valid only in 'ggplot' style.",
            call. = FALSE)

  if (is.null(concentration) && addlegend)
    warning("'addlegend' argument is valid only when 'concentration' is not null.",
            call. = FALSE)

  if (pool.replicate) {
    # agregate by sum of replicate
    x <- cbind(aggregate(Nsurv ~ time + conc, x, sum),
               replicate = 1)
  }

  x <- as.data.frame(x)

  if (is.null(concentration)) {
    survDataPlotFull(x, xlab, ylab, style, remove.someLabels)
  }  else {
    survDataPlotFixedConc(x, xlab, ylab, main, concentration,
                          style, addlegend, remove.someLabels)
  }
}


# [ReplicateIndex(data)] builds a list of indices, each one named after
# a replicate of [data], thus providing a dictionary from replicate names to
# integer keys.
ReplicateIndex <- function(data) {
  replicate <- unique(data$replicate)
  r <- as.list(seq(1, length(replicate)))
  names(r) <- as.character(replicate)
  return(r)
}


# General full plot: one subplot for each concentration, and one color for
# each replicate (for generic graphics)
dataPlotFullGeneric <- function(data, xlab, ylab, resp) {
  replicate.index <- ReplicateIndex(data)

  # creation of a vector of colors
  colors <- rainbow(length(unique(data$replicate)))
  pchs <- as.numeric(unique(data$replicate))
  # split of the graphical window in subplots
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar)) 
  par(mfrow = plotMatrixGeometry(length(unique(data$conc))))

  by(data, data$conc, function(x) {
    x <- as.data.frame(x)
    # background
    plot(x$time, rep(0, length(x$time)),
         xlab = xlab,
         ylab = ylab,
         ylim = c(0, max(x[, resp])),
         type = "n",
         col = 'white',
         xaxt = "n",
         yaxt = "n")

    # axis
    axis(side = 1, at = sort(unique(as.data.frame(x)[, "time"])))
    axis(side = 2, at = unique(round(pretty(c(0, max(x[, resp]))))))


    # lines and points
    by(x, x$replicate, function(y) {
      index <- replicate.index[[y$replicate[1]]]
      lines(y$time, y[, resp],
            type = "l",
            col = colors[index])
      points(y$time, y[, resp],
             pch = pchs[index],
             col = colors[index])
    })

    # title
    title(paste("Conc: ", unique(x$conc), sep = ""))
  })

  par(mfrow = c(1, 1))
}

# general full plot (ggplot variant): one subplot for each concentration,
# and one color for each replicate
#' @import ggplot2
dataPlotFullGG <- function(data, xlab, ylab, resp, remove.someLabels) {

  data <- as.data.frame(data)
  time = NULL
  Nsurv = NULL

  data$response <- data[, resp]

  # create ggplot object Nsurv / time / replicate / conc
  fg <- ggplot(data = data, aes(time, response, colour = factor(replicate))) +
        geom_point() +
        geom_line()  +
        labs(x = xlab, y = ylab) +
        facet_wrap(~conc, ncol = 2) +
        scale_x_continuous(breaks = unique(data$time),
                           labels = if (remove.someLabels) {
                             exclude_labels(unique(data$time))
                           } else {
                             unique(data$time)
                           }
        ) +
        scale_y_continuous(breaks = unique(round(pretty(c(0, max(data$response)))))) +
        expand_limits(x = 0, y = 0) +
        theme_minimal()

   fd <- fg + theme(legend.position = "none") # remove legend

   return(fd)

}

dataPlotFull <- function(data, xlab, ylab, resp, style = "generic",
                         remove.someLabels = FALSE) {

  if (missing(xlab)) xlab <- "Time"

  if (style == "generic")
    dataPlotFullGeneric(data, xlab, ylab, resp)
  else if (style == "ggplot")
    dataPlotFullGG(data, xlab, ylab, resp, remove.someLabels)
  else stop("Unknown plot style")
}

survDataPlotFull <- function(data, xlab, ylab,
                             style = "ggplot",
                             remove.someLabels = FALSE) {
  dataPlotFull(data, xlab, ylab, "Nsurv", style, remove.someLabels)
}

dataPlotFixedConc <- function(x,
                              xlab,
                              ylab,
                              main,
                              resp,
                              concentration,
                              style = "generic",
                              addlegend = FALSE,
                              remove.someLabels = FALSE) {

  x <- as.data.frame(x) # x is interpreted as a tibble

    if (missing(xlab)) xlab <- "Time"

  legend.position <- ifelse(resp == "Nsurv", "bottomleft", "topleft")

  # check concentration value
  if (!concentration %in% x$conc)
    stop("The argument [concentration] should correspond to one of the tested concentrations")

  # select the concentration
  x <- filter(x, x$conc == concentration)

  # vector color
  x$color <- as.numeric(as.factor(x$replicate))

  if (style == "generic") {
    plot(x$time, x[, resp],
         type = "n",
         xaxt = "n",
         yaxt = "n",
         main = main,
         xlim = range(x$time),
         ylim = c(0, max(x[, resp])),
         xlab = xlab,
         ylab = ylab)

    # one line by replicate
    by(x, list(x$replicate),
       function(x) {
         lines(x$time, x[,resp], # lines
               col = x$color)
         points(x$time, x[,resp], # points
                pch = 16,
                col = x$color)
       })

    # axis
    axis(side = 1, at = sort(unique(x[, "time"])))
    axis(side = 2, at = unique(round(pretty(c(0, max(x[, resp]))))))

    if (addlegend && !length(unique(x$replicate)) == 1) {
      legend(legend.position, legend = unique(x$replicate) ,
             col = unique(x$color),
             pch = 16,
             lty = 1)
    }
  }
  else if (style == "ggplot") {
    x$response <- x[,resp]

    if (length(unique(x$replicate)) == 1) {
      df <- ggplot(x, aes(x = time, y = response))
    } else {
      df <- ggplot(x, aes(x = time, y = response,
                          color = factor(replicate),
                          group = replicate))
    }
    fd <- df + geom_line() + geom_point() + ggtitle(main) +
      theme_minimal() +
      labs(x = xlab,
           y = ylab) +
      scale_color_hue("Replicate") +
      scale_x_continuous(breaks = unique(x$time),
                         labels = if (remove.someLabels) {
                           exclude_labels(unique(x$time))
                         } else {
                           unique(x$time)
                         }) +
      scale_y_continuous(breaks = unique(round(pretty(c(0, max(x$response)))))) +
      expand_limits(x = 0, y = 0)

    if (addlegend) {# only if pool.replicate == FALSE
      fd
    } else {
      fd + theme(legend.position = "none") # remove legend
    }
  }
  else stop("Unknown plot style")
}

survDataPlotFixedConc <- function(x,
                                  xlab,
                                  ylab,
                                  main,
                                  concentration,
                                  style = "generic",
                                  addlegend = FALSE,
                                  remove.someLabels = FALSE) {

  dataPlotFixedConc(x, xlab, ylab, main, "Nsurv", concentration,
                    style, addlegend, remove.someLabels)
}
