repro.convergence <-
		function(out, 
				trace = TRUE,
				density = TRUE,
				autocorr = TRUE,
				type = "generic")

#	user function	#
####### Convergence test
{
	# INPUT:
	# - out$mcmc: list of estimated parameters for the model with each item representing 
	#	a chains
	#   -each items: size number of parameters columns and Niter rows
	# - trace: option to plot the trace  of mcmc
	# - density: option to plot the density of mcmc
	# - autocorr: plot autocorrelation of mcmc 
	# - type: generic ou ggplot	
	# OUTPUT:
	#	- Gelman and Rubin test A list with 
	#		psrf: The point estimates of the potential scale reduction factor
	# 		mpsrf:  The point estimate of the multivariate potential scale reduction factor
	
	#	test class object
	if (class(out) != "repro.fit")
		stop("The object passed in argument 'out' is not of class 'repro.fit'!\n")
	
	#	MCMC object
	mcmc <- out$mcmc
	#	diagnostic
	
	GelRubmulti <- gelman.diag(mcmc)$mpsrf
	GelRubesti <- gelman.diag(mcmc)$psrf[,"Point est."]
	
	####	output	
	
	#	return the psrf and mprsf value	
	cat("Gelman and Rubin:\n")
	cat("Potential scale reduction factor for each parameter:\n")
	print(GelRubesti)
	cat("\nMultivariate potential scale reduction factor:\n", GelRubmulti,"\n")
	
	
	#	generic plot
	if (type == "generic"){
		#	trace and density
		if (trace == TRUE || density == TRUE)
			plot(mcmc, trace = trace, density = density)
		#	autocorrelation
		if (autocorr == TRUE){
			dev.new()
			autocorr.plot(mcmc)
		}
	}
	
	#	ggplot loading
	
	if (type == "ggplot"){
		require(ggmcmc, quietly = TRUE)
		require(gridExtra, quietly = TRUE)
		
		# creat ggs objects
		D <- ggs(mcmc)
		if (trace == TRUE) trp <- ggs_traceplot(D)
		if (density == TRUE) dns <- ggs_density(D)
		if (autocorr == TRUE) atc <- ggs_autocorrelation(D)
		if(trace == FALSE || density == FALSE || autocorr == FALSE)		
			blank <- grid.rect(gp=gpar(col="white"))
		
		do.call(grid.arrange, list(if(trace == TRUE) trp else blank, 
						if(density == TRUE) dns else blank, 
						if(autocorr == TRUE) atc else blank, ncol = 2))
		
	}
	return(invisible(list(mpsrf = GelRubmulti,
							psrf = GelRubesti)))
}
