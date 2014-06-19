summary.repro.fit <-
		function(object,...)

{
	set.seed(1234)
	
	# transformed data at time target.time
	#	distribution
	# d
	# normal
	sumdistd <- rnorm(n = object$n.iter,
			mean = object$param.prior$meand,
			sd = 1/sqrt(object$param.prior$taud))				
	
	# b
	# uniform
	sumdistlog10b <- runif(n = object$n.iter,
			min = object$param.prior$log10bmin,
			max = object$param.prior$log10bmax)
	
	# e
	# normal
	sumdistlog10e <- rnorm(n = object$n.iter,
			mean = object$param.prior$meanlog10e,
			sd = 1/sqrt(object$param.prior$taulog10e))
	
	if(object$model.label == "P"){
		res <- cbind(sumdistd, sumdistlog10b, sumdistlog10e)
		res.names <- c("d", "log10b", "log10e")
	}
	
	if(object$model.label == "GP"){	
		# omega
		# uniform
		sumdistlog10omega <- runif(n = object$n.iter,
				min = object$param.prior$log10omegamin,
				max = object$param.prior$log10omegamax)
		res <- cbind(sumdistd, sumdistlog10b, sumdistlog10e, sumdistlog10omega)
		res.names <- c("d", "log10b", "log10e", "log10omega")
	}
	
	#	quantile
	med <- apply(res, 2, function(x) quantile(x, probs = 0.5))
	Q2.5 <- apply(res, 2, function(x) quantile(x, probs = 0.025))
	Q97.5 <- apply(res, 2, function(x) quantile(x, probs = 0.975))
	ans1 <-  round(data.frame(med, Q2.5, Q97.5,
					row.names = res.names), digits = 3)
	colnames(ans1) <- c("50%", "2.5%", "97.5%")
	
	# parameters median and quantile obtained from the mcmc 
	ans2 <- round(object$estim.par, digits = 3)
	colnames(ans2) <- c("50%", "2.5%", "97.5%")
	
	# estimated ECx and their CIs 95%
	ans3 <- round(object$estim.ECx, digits = 3)
	colnames(ans3) <- c("50%", "2.5%", "97.5%")
	
	# summary of mcmc posteriors from coda file
	ans4 <- summary(object$mcmc, quantiles = c(0.5, 0.025, 0.975))
	
	cat("Summary: \n\n")
	if(object$model.label == "GP")
		cat("The log-logistic model with a Gamma Poisson stochastic part
						was used.\n\n")
	if(object$model.label == "P")
		cat("The log-logistic model with a Poisson stochastic part 
						was used.\n\n")
	
	cat("Quantiles of priors on parameters: \n\n")
	print(ans1)
	cat("\nEstimated parameters: \n\n")
	print(ans2)
	cat("\nEstimated ECx:\n\n")
	print(ans3)
	cat("\nMCMC estimation: \n")	
	print(ans4)
}
