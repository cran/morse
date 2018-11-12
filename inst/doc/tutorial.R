## ----include=FALSE, echo=FALSE-------------------------------------------
knitr::opts_chunk$set(fig.width = 7,
                      fig.height = 4,
                      cache = TRUE)

## ----checkNsurvPRED, cache=TRUE------------------------------------------
predict_Nsurv_check(predict_Nsurv_PRZ_SD_cstTOvar)

## ----plotPredict_Nsurv, cache=TRUE---------------------------------------
plot(predict_Nsurv_PRZ_SD_cstTOvar)

## ----ppcPredict_Nsurv, cache=TRUE----------------------------------------
ppc(predict_Nsurv_PRZ_SD_cstTOvar)

