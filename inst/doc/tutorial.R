## ----include=FALSE, echo=FALSE-------------------------------------------

knitr::opts_chunk$set(fig.width = 7, 
                      fig.height = 4, 
                      cache = TRUE)


## ---- echo=FALSE, cache=TRUE, results='hide'-----------------------------
library(morse)
library(ggplot2)

## ---- cache=TRUE---------------------------------------------------------
data(cadmium2)
survDataCheck(cadmium2, diagnosis.plot = TRUE)

## ---- cache=TRUE---------------------------------------------------------
dat <- survData(cadmium2)
head(dat)

## ---- cache=TRUE---------------------------------------------------------
plot(dat, style = "ggplot", pool.replicate = FALSE)

## ---- cache=TRUE---------------------------------------------------------
plot(dat, concentration = 124, addlegend = TRUE,
     pool.replicate = FALSE, style = "ggplot")

## ---- cache=TRUE---------------------------------------------------------
plot(dat, target.time = 21, pool.replicate = FALSE, style = "ggplot",
     addlegend = TRUE)

## ---- cache=TRUE---------------------------------------------------------
plot(dat, concentration = 232, target.time = 21)

## ---- cache=TRUE---------------------------------------------------------
summary(dat)

## ---- results="hide", cache=TRUE-----------------------------------------
fit <- survFitTT(dat,
                lcx = c(10, 20, 30, 40, 50))

## ---- cache=TRUE---------------------------------------------------------
summary(fit)

## ---- cache=TRUE---------------------------------------------------------
plot(fit, log.scale = TRUE, ci = TRUE, style = "ggplot",
     addlegend = TRUE)

## ----cache=TRUE, eval=FALSE----------------------------------------------
#  ppc(out2, style = "ggplot")

## ---- results="hide", cache=TRUE-----------------------------------------
data("cadmium1")
fit <- survFitTT(survData(cadmium1),
                lcx = c(10, 20, 30, 40, 50))
plot(fit, log.scale = TRUE, ci = TRUE, style = "ggplot",
     addlegend = TRUE)

## ---- cache=TRUE---------------------------------------------------------
# (1) load dataset
data(cadmium2)

# (2) check structure and integrity of the dataset
reproDataCheck(cadmium2, diagnosis.plot = TRUE)

# (3) create a `reproData` object
dat <- reproData(cadmium2)

# (4) represent the cumulated number of offspring as a function time
plot(dat, style = "ggplot", pool.replicate = FALSE)
plot(dat, target.time = 21, addlegend = TRUE, style = "ggplot",
     pool.replicate = FALSE)
plot(dat, concentration = 124, addlegend = TRUE, style = "ggplot",
     pool.replicate = FALSE)
plot(dat, concentration = 124, target.time = 21, style = "ggplot")

# (5) check information on the experimental design
summary(dat)

# (6) fit an exposure-response model at target-time
fit <- reproFitTT(dat, stoc.part = "bestfit",
                  ecx = c(10, 20, 30, 40, 50),
                  quiet = TRUE)
plot(fit, log.scale = TRUE, ci = TRUE, 
     style = "ggplot", addlegend = TRUE)

## ---- cache=TRUE---------------------------------------------------------
summary(fit)

## ---- cache=TRUE---------------------------------------------------------
dat <- reproData(cadmium2)
plot(as.survData(dat))

## ---- cache=TRUE, eval=FALSE---------------------------------------------
#  ppc(out, style = "ggplot")

