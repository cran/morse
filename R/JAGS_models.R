#
# jags_TKTD_cstSD is the JAGS model used in `MORSE` 2.Y.Z packages in the function `survFitTKTD()`
#
jags_TKTD_cstSD <-
  "model {
    ##--------------------- priors

    kd_taulog10 <- 1 / kd_sdlog10^2
    hb_taulog10 <- 1 / hb_sdlog10^2
    kk_taulog10 <- 1 / kk_sdlog10^2
    z_taulog10 <- 1 / z_sdlog10^2

    kk_log10 ~ dnorm(kk_meanlog10, kk_taulog10)
    z_log10  ~ dnorm(z_meanlog10 , z_taulog10)
    kd_log10 ~ dnorm(kd_meanlog10, kd_taulog10)
    hb_log10 ~ dnorm(hb_meanlog10, hb_taulog10)

    ##-------------------- parameter transformation
    kd <- 10**kd_log10
    hb <- 10**hb_log10
    kk <- 10**kk_log10
    z  <- 10**z_log10

    bigtime <- max(time[1:n_data]) + 1

    ##------------------- Computation of the likelihood

    for (i in 1:n_data){

      tz[i] <- ifelse(conc[i] > z, -1/kd * log( 1- R[i]), bigtime)
      R[i] <- ifelse(conc[i] > z, z/xcor[i], 0.1)
      xcor[i] <- ifelse(conc[i] > 0, conc[i], 10)

      tref[i] <- max(tprec[i], tz[i])
      psurv[i] <- exp(-hb * (time[i] - tprec[i]) + ifelse(time[i] > tz[i], -kk * ((conc[i] - z) * (time[i] - tref[i]) + conc[i]/kd * ( exp(-kd * time[i]) - exp(-kd * tref[i]))), 0))
     
      Nsurv[i] ~ dbin(psurv[i] , Nprec[i])
  }

}"


#
# IT model with log-logistic function.
#

jags_TKTD_cstIT <-"
model {

  ##------------------------------------------------------ priors

  kd_taulog10 <- 1 / kd_sdlog10^2
  hb_taulog10 <- 1 / hb_sdlog10^2
  alpha_taulog10 <- 1 / alpha_sdlog10^2

  kd_log10 ~ dnorm(kd_meanlog10, kd_taulog10)
  hb_log10 ~ dnorm(hb_meanlog10, hb_taulog10)

  alpha_log10 ~ dnorm(alpha_meanlog10, alpha_taulog10)
  beta_log10 ~ dunif(beta_minlog10, beta_maxlog10) # unif -2, 2

  ##---------------------------------------- parameter transformation

  kd <- 10**kd_log10
  hb <- 10**hb_log10

  alpha <- 10**alpha_log10
  beta <- 10**beta_log10

  ##------------------------------------ model

  for( i in 1:n_data){

    D[replicate_ID[i], time_ID[i]] <- conc[i] * ( 1 - exp( - kd * time[i] ))

    D_max[replicate_ID[i], time_ID[i]] <- max(D[replicate_ID[i], 1:time_ID[i]])

    F[i]  <- D_max[replicate_ID[i], time_ID[i]]**beta / ( D_max[replicate_ID[i], time_ID[i]]**beta + alpha**beta )

    psurv[i] <-  exp(-hb * time[i]) * (1- F[i])
    
    Nsurv[i] ~ dbin(psurv[i]/psurv[i_prec[i]] , Nprec[i])

  }

}"


#
# SD with non-constant concentration exposure
#


jags_TKTD_varSD <-
  "model {

  ##------------------------------------------ parameter transformation
  kd_taulog10 <- 1 / kd_sdlog10^2
  hb_taulog10 <- 1 / hb_sdlog10^2
  kk_taulog10 <- 1 / kk_sdlog10^2
  z_taulog10 <- 1 / z_sdlog10^2

  ##------- priors
  kd_log10 ~ dnorm(kd_meanlog10, kd_taulog10)
  hb_log10 ~ dnorm(hb_meanlog10, hb_taulog10)
  kk_log10 ~ dnorm(kk_meanlog10, kk_taulog10)
  z_log10  ~ dnorm(z_meanlog10 , z_taulog10)

  ##----- parameter transformation
  kd <- 10**kd_log10
  hb <- 10**hb_log10
  kk <- 10**kk_log10
  z  <- 10**z_log10

  ##------ Computation of the likelihood

  for( i in 1:n_data_long){

    ###-------------- midpoint method

    conc_midpoint[replicate_ID_long[i], time_ID_long[i]] <-  (exp(kd * time_long[i]) * conc_long[i] + exp(kd * tprec_long[i])
    * concprec_long[i]) / 2 * (time_long[i] - tprec_long[i])

    D_int[replicate_ID_long[i], time_ID_long[i]] <- kd * exp(-kd * time_long[i]) *
    sum( conc_midpoint[replicate_ID_long[i], 1:time_ID_long[i]] )

    h[replicate_ID_long[i], time_ID_long[i]] <- kk * max(0, D_int[replicate_ID_long[i], time_ID_long[i]] - z) + hb

    h_midPoint[replicate_ID_long[i], time_ID_long[i]] <- (h[replicate_ID_long[i], time_ID_long[i]] +
    h[replicate_ID_long[i], tprec_ID_long[i]])/2 * (time_long[i] - tprec_long[i])

    H_int[replicate_ID_long[i], time_ID_long[i]] <- sum( h_midPoint[replicate_ID_long[i], 1:time_ID_long[i]] )

  }
  for( i in 1:n_data_red){

    H[replicate_ID[i], time_ID_red[i]]  <- H_int[replicate_ID[i], time_ID_long_red[i]]

    psurv[i] = exp( - H[replicate_ID[i], time_ID_red[i]])

    Nsurv[i] ~ dbin(psurv[i]/psurv[i_prec[i]] , Nprec[i])

 ## ---------------------- generated data 
 
    Nsurv_ppc[i] ~ dbin(psurv[i]/psurv[i_prec[i]] , Nprec[i]) 
 
  } 
 
}"


#
# IT model with log-logistic function, and non-constant exposure concentration
#

jags_TKTD_varIT <-"model {
  #------------------------------------------ parameter transformation
  kd_taulog10 <- 1 / kd_sdlog10^2
  hb_taulog10 <- 1 / hb_sdlog10^2
  alpha_taulog10 <- 1 / alpha_sdlog10^2

  #-------------------------------------------------------- priors
  kd_log10 ~ dnorm(kd_meanlog10, kd_taulog10)
  hb_log10 ~ dnorm(hb_meanlog10, hb_taulog10)

  alpha_log10 ~ dnorm(alpha_meanlog10, alpha_taulog10)
  beta_log10 ~ dunif(beta_minlog10, beta_maxlog10)

  #------------------------------------------ parameter transformation

  kd <- 10**kd_log10
  hb <- 10**hb_log10

  alpha <- 10**alpha_log10
  beta <- 10**beta_log10

  ##------------------------------------ model

  for( i in 1:n_data_long){

    ###------------ midpoint method

    diff.int[replicate_ID_long[i], time_ID_long[i]] <-  (exp(kd * time_long[i]) * conc_long[i] + exp(kd * tprec_long[i])
    * concprec_long[i]) / 2 * (time_long[i] - tprec_long[i])

    D_int[replicate_ID_long[i], time_ID_long[i]] <- kd * exp(-kd * time_long[i]) *
    sum( diff.int[replicate_ID_long[i], 1:time_ID_long[i]] )

  }

  for( i in 1:n_data_red){

    D[replicate_ID[i], time_ID_red[i]]  <- D_int[replicate_ID[i], time_ID_long_red[i]]

    D_max[replicate_ID[i], time_ID_red[i]] <- max(D[replicate_ID[i],1:time_ID_red[i]])

    F[i]  <- D_max[replicate_ID[i], time_ID_red[i]]**beta / ( D_max[replicate_ID[i], time_ID_red[i]]**beta + alpha**beta )

    psurv[i] <-  exp(-hb * time[i]) * (1- F[i])

    Nsurv[i] ~ dbin(psurv[i]/psurv[i_prec[i]] , Nprec[i])


  ## ---------------------- generated data 
 
    Nsurv_ppc[i] ~ dbin(psurv[i]/psurv[i_prec[i]] , Nprec[i]) 
 
  } 
 
  # ### initialization is requires to use 'Nsurv_sim[i-1]' require in JAGS language (avoid auto-loop issue). 
  # Nsurv_sim[1] ~ dbin(psurv[1]/psurv[1], Nprec[1]) 
  # for( i in 2:n_data_red){ 
  #   Nsurv_sim[i] ~ dbin(psurv[i]/psurv[i_prec[i]], ifelse( i == i_prec[i], Nprec[i], Nsurv_sim[i-1])) 
  # }

}"