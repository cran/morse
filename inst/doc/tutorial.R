## ----include=FALSE, echo=FALSE-------------------------------------------
knitr::opts_chunk$set(fig.width = 7,
                      fig.height = 4,
                      cache = TRUE)

## ---- echo=TRUE, cache=TRUE----------------------------------------------
# (6) fit the TK-TD model SD
fit_varSD <- survFit(dat, quiet = TRUE, model_type = "SD")

## ---- echo=TRUE, cache=TRUE----------------------------------------------
# fit a TK-TD model IT
fit_varIT <- survFit(dat, quiet = TRUE, model_type = "IT")

