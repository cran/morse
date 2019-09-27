
# object  = survFit_SD
# 
# listMCMC = object$mcmc
# subListMCMC_N = listMCMC[, grep("Nsurv_sim", varnames(listMCMC))]
# subListMCMC_p = listMCMC[, grep("psurv", varnames(listMCMC))]
# 
# DFsubListMCMC_N = as.data.frame(do.call("cbind", subListMCMC_N))
# DFsubListMCMC_p = as.data.frame(do.call("cbind", subListMCMC_p))
# 
# vecNiter = table(object$jags.data$replicate)
# 
# tprec = object$jags.data$tprec
# time = object$jags.data$time
# 
# 
# prodInit = c()
# diff_N = list()
# diff_p = list()
# for(i in 1:length(vecNiter)){
#    prodInit[i] = factorial(DFsubListMCMC_N[, 1]) *  factorial(DFsubListMCMC_N[, vecNiter[i]])^factorial(DFsubListMCMC_N[, vecNiter[i]]) / factorial(DFsubListMCMC_N[, vecNiter[i]])
#   for(j in 1:vecNiter[i]){
#     diff_N[i][j] = DFsubListMCMC_N[, tprec+1] - DFsubListMCMC_N[, time+1]
#     diff_p[i][j] = DFsubListMCMC_p[, tprec+1] - DFsubListMCMC_p[, time+1]
#   }
# }
