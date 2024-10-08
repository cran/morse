#' Predict method for \code{survFit} objects
#' 
#' This is a \code{method} to replace function \code{predict_Nsurv} used on \code{survFit}
#' object when computing issues happen. \code{predict_nsurv_ode} uses the \code{deSolve}
#' library to improve robustness. However, time to compute may be longer.
#' 
#' @name predict
#' 
#' @param object An object of class \code{survFit}.
#' @param data_predict A dataframe with three columns \code{time}, \code{conc} and \code{replicate}
#'  used for prediction. If \code{NULL}, prediction is based on \code{x} object of 
#'  class \code{survFit} used for fitting.
#' @param spaghetti If \code{TRUE}, return a set of survival curves using
#' parameters drawn from the posterior distribution.
#' @param mcmc_size Can be used to reduce the number of mcmc samples in order to speed up
#'  the computation. \code{mcmc_size} is the number of selected iterations for one chain. Default
#'  is 1000. If all MCMC is wanted, set argument to \code{NULL}.
#' @param hb_value If \code{TRUE}, the background mortality \code{hb} is taken into account from the posterior.
#' If \code{FALSE}, parameter \code{hb} is set to 0. The default is \code{TRUE}.
#' @param  hb_valueFORCED If \code{hb_value} is \code{FALSE}, it fix \code{hb}.
#' @param extend_time Length of time points interpolated with variable exposure profiles.
#' @param interpolate_length Length of the time sequence for which output is wanted.
#' @param interpolate_method The interpolation method for concentration. See package \code{deSolve} for details.
#' Default is \code{linear}.
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @return an object of class \code{predict_Nsurv_ode}.
#' 
#' @export
#' 
predict_Nsurv_ode <- function(object,
                               data_predict,
                               spaghetti,
                               mcmc_size,
                               hb_value,
                               hb_valueFORCED,
                               extend_time,
                               interpolate_length,
                               interpolate_method,
                               ...){
  UseMethod("predict_Nsurv_ode")
}

#' @import deSolve
#' @importFrom stats approxfun
#' 
#' @rdname predict
#' @return a \code{list} of \code{data.frame} with the quantiles of outputs in
#' \code{df_quantiles} or all the MCMC chaines \code{df_spaghetti}
#' 
#' @export
#'
predict_Nsurv_ode.survFit <- function(object,
                                  data_predict = NULL,
                                  spaghetti = FALSE,
                                  mcmc_size = 1000,
                                  hb_value = TRUE,
                                  hb_valueFORCED = NA,
                                  extend_time = 100,
                                  interpolate_length = NULL,
                                  interpolate_method = "linear",
                                  ...) {
  x <- object # Renaming to satisfy CRAN checks on S3 methods
  # arguments should be named the same when declaring a
  # method and its instantiations
  
  if(!("Nsurv" %in% colnames(data_predict))){
    warning("Please provide a column 'Nsurv' in the 'data_predict' argument to have
            prediction on the Number of survivor.")
  }
  
  message("Note that computing can be quite long (several minutes).
  Tips: To reduce that time you can reduce Number of MCMC chains (default mcmc_size is set to 1000).")
  
  # Initialisation
  mcmc <- x$mcmc
  model_type <- x$model_type
  extend_time <- extend_time
  
  if(is.null(data_predict)){
    if("survFitVarExp" %in% class(x)){
      x_interpolate = data.frame(
        time = x$jags.data$time_long,
        conc = x$jags.data$conc_long,
        replicate = x$jags.data$replicate_long)
    } else{
      data_predict = data.frame(
        time = x$jags.data$time,
        conc = x$jags.data$conc,
        replicate = x$jags.data$replicate,
        Nsurv = x$jags.data$Nsurv)
      
      x_interpolate <- predict_interpolate(data_predict,  extend_time = extend_time) %>%
        dplyr::arrange(replicate, time)
    }
  }
  if(!is.null(data_predict)){
    x_interpolate <- predict_interpolate(data_predict,  extend_time = extend_time) %>%
      dplyr::arrange(replicate, time)
  }
  
  df <- data.frame(
    time = x_interpolate$time,
    conc = x_interpolate$conc,
    replicate = x_interpolate$replicate)
  
  unique_replicate <- unique(df$replicate)
  
  ls_time <- list()
  ls_conc <- list()
  
  for(i in 1:length(unique_replicate)){
    
    ls_time[[i]] <- dplyr::filter(df, replicate == unique_replicate[i])$time
    ls_conc[[i]] <- dplyr::filter(df, replicate == unique_replicate[i])$conc
    
  }
  
  # ------- Computing
  mcmc.samples = mcmc
  
  if(!is.null(mcmc_size)){
    reduc_tab = lapply(mcmc.samples, "[", 
                       seq(1, nrow(mcmc.samples[[1]]), length = mcmc_size),
                       1:ncol(mcmc.samples[[1]]))
    mcmc.samples = reduc_tab
  }
  
  mctot = do.call("rbind", mcmc.samples)
  mcmc_size = nrow(mctot)
  
  kd = 10^mctot[, "kd_log10"]
  
  if(hb_value == TRUE){
    # "hb" is not in survFit object of morse <v3.2.0
    if("hb" %in% colnames(mctot)){
      hb <- mctot[, "hb"]  
    } else{ hb <- 10^mctot[, "hb_log10"] }
  } else if(hb_value == FALSE){
    if(is.na(hb_valueFORCED)){
      if(is.na(x$hb_valueFIXED)){
        stop("Please provide value for `hb` using `hb_valueFORCED`.")
      } else{
        hb <- rep(x$hb_valueFIXED, nrow(mctot))
      } 
    } else{
      hb <- rep(hb_valueFORCED, nrow(mctot))
    }
  }
  
  k = 1:length(unique_replicate)
  
  if (model_type == "SD") {
    kk <- 10^mctot[, "kk_log10"]
    z <- 10^mctot[, "z_log10"]
    
    dtheo = lapply(k, function(kit) { # For each replicate
      SurvSD_ode(Cw = ls_conc[[kit]],
                 time = ls_time[[kit]],
                 replicate = unique_replicate[kit],
                 kk = kk,
                 kd = kd,
                 hb = hb,
                 z = z,
                 mcmc_size = mcmc_size,
                 interpolate_length = interpolate_length,
                 interpolate_method = interpolate_method)
    })
    
  }
  if (model_type == "IT") {
    
    alpha <- 10^mctot[, "alpha_log10"]
    beta <- 10^mctot[, "beta_log10"]
    
    dtheo = lapply(k, function(kit) { # For each replicate
      SurvIT_ode(Cw = ls_conc[[kit]],
                 time = ls_time[[kit]],
                 replicate = unique_replicate[kit],
                 kd = kd,
                 hb = hb,
                 alpha = alpha,
                 beta = beta,
                 mcmc_size = mcmc_size,
                 interpolate_length = interpolate_length,
                 interpolate_method = interpolate_method)
    })
  }
  # Transpose
  df_theo <- do.call("rbind", dtheo)
  # dtheo <- do.call("rbind", lapply(dtheo, t))
  
  # Computing Nsurv
  df_mcmc <- as_tibble(do.call("rbind", x$mcmc))
  NsurvPred_valid <- select(df_mcmc, contains("Nsurv_sim"))
  NsurvPred_check <- select(df_mcmc, contains("Nsurv_ppc"))
  
  if (is.null(data_predict) &
     # The following condition are always true for survFit done after morse v3.2.0 !
     ncol(NsurvPred_valid) > 0 &
     ncol(NsurvPred_check) > 0) {
    
    df_quantile <- data.frame(
      time = data_predict$time,
      conc = data_predict$conc,
      replicate = data_predict$replicate,
      Nsurv = data_predict$Nsurv,
      Nsurv_q50_check = apply(NsurvPred_check, 1, quantile, probs = 0.5, na.rm = TRUE),
      Nsurv_qinf95_check = apply(NsurvPred_check, 1, quantile, probs = 0.025, na.rm = TRUE),
      Nsurv_qsup95_check = apply(NsurvPred_check, 1, quantile, probs = 0.975, na.rm = TRUE),
      Nsurv_q50_valid = apply(NsurvPred_valid, 1, quantile, probs = 0.5, na.rm = TRUE),
      Nsurv_qinf95_valid = apply(NsurvPred_valid, 1, quantile, probs = 0.025, na.rm = TRUE),
      Nsurv_qsup95_valid = apply(NsurvPred_valid, 1, quantile, probs = 0.975, na.rm = TRUE))
    
  } else{
    # --------------------

    df_psurv <- as_tibble(df_theo) %>%
      select(-conc) %>%
      mutate(time = df$time,
             replicate = df$replicate)
    
    df_filter <- dplyr::inner_join(df_psurv, data_predict, by = c("replicate", "time")) %>%
      filter(!is.na(Nsurv)) %>%
      group_by(replicate) %>%
      arrange(replicate, time) %>%
      mutate(Nprec = ifelse(time == min(time), Nsurv, lag(Nsurv)),
             iter = row_number(),
             iter_prec = ifelse(time == min(time), iter, lag(iter))) %>%
      ungroup()
    
    mat_psurv <- df_filter %>%
      select(-c("time", "conc", "replicate",
                 "q50", "qinf95", "qsup95", 
                 "Nsurv", "Nprec", "iter", "iter_prec")) %>%
      as.matrix()
    
    ncol_NsurvPred <- ncol(mat_psurv)
    nrow_NsurvPred <- nrow(mat_psurv)
    iter = df_filter$iter
    iter_prec = df_filter$iter_prec
    
    NsurvPred_valid <- matrix(ncol = ncol_NsurvPred, nrow = nrow(mat_psurv))

    Nprec <- cbind(df_filter$Nprec)[, rep(1,ncol_NsurvPred)]
    
    mat_psurv_prec = matrix(ncol = ncol_NsurvPred, nrow = nrow_NsurvPred)
    for (i in 1:nrow_NsurvPred) {
      if (iter[i] == iter_prec[i]) {
        mat_psurv_prec[i,] = mat_psurv[i,]
      } else{
        mat_psurv_prec[i,] = mat_psurv[i-1,]
      }
    }
    mat_pSurv_ratio = mat_psurv / mat_psurv_prec
    
    NsurvPred_check_vector = rbinom(ncol_NsurvPred*nrow_NsurvPred,
                                    size = Nprec,
                                    prob =  mat_pSurv_ratio)
    NsurvPred_check = matrix(NsurvPred_check_vector, byrow = FALSE, nrow = nrow_NsurvPred)
    

    NsurvPred_valid[1, ] = rep(Nprec[1], ncol_NsurvPred)
    for (i in 2:nrow(mat_psurv)) {
      if (iter[i] == iter_prec[i]) {
        NsurvPred_valid[i,] = NsurvPred_check[i,]
      } else{
        NsurvPred_valid[i,] = rbinom(ncol_NsurvPred,
                                     size = NsurvPred_valid[i - 1,],
                                     prob = mat_pSurv_ratio[i,])
      }
    }
    
    
    df_quantile <- data.frame(time = df_filter$time,
                              conc = df_filter$conc,
                              replicate = df_filter$replicate,
                              Nsurv = df_filter$Nsurv,
                              Nsurv_q50_check = apply(NsurvPred_check, 1, quantile, probs = 0.5, na.rm = TRUE),
                              Nsurv_qinf95_check = apply(NsurvPred_check, 1, quantile, probs = 0.025, na.rm = TRUE),
                              Nsurv_qsup95_check = apply(NsurvPred_check, 1, quantile, probs = 0.975, na.rm = TRUE),
                              Nsurv_q50_valid = apply(NsurvPred_valid, 1, quantile, probs = 0.5, na.rm = TRUE),
                              Nsurv_qinf95_valid = apply(NsurvPred_valid, 1, quantile, probs = 0.025, na.rm = TRUE),
                              Nsurv_qsup95_valid = apply(NsurvPred_valid, 1, quantile, probs = 0.975, na.rm = TRUE))
    
  } 
  
  if (spaghetti == TRUE) {
    random_column <- sample(1:ncol(NsurvPred_valid), size = round(10/100 * ncol(NsurvPred_valid)))
    df_spaghetti <- as_tibble(NsurvPred_valid[, random_column]) %>%
      mutate(time = data_predict$time,
             conc = data_predict$conc,
             replicate = data_predict$replicate,
             Nsurv = data_predict$Nsurv)
  } else df_spaghetti <- NULL
  
  #ls_check_on_Nsurv <- check_on_Nsurv(df_quantile)
  
  return_object <- list(df_quantile = df_quantile,
                        df_spaghetti = df_spaghetti)
  
  class(return_object) <- c(class(return_object), "survFitPredict_Nsurv")
  
  return(return_object)
  
  }



