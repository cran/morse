repro.fit <-
		function(rdata,  
				n.chains = 3, 
				quiet = FALSE) {
	
	#	user function	#
	# INPUT
	# - rdata from TrasnformData function
	#   - output of repro.data function
	# - n.chain : number of mcmc chains
	# - quiet : silent option
	# OUTPUT
	# - raw.data
	# - transformed.data 
	# - mortality
	# - model.label
	# - model 
	# - n.chains
	# - n.iter
	# - n.burnin
	# - n.thin 
	# - param.prior
	# - mcmc
	# - DIC
	# - estim.par
	# - estim.ECx
	
	#	test class object
	if (class(rdata) != "repro.data")
		stop("The object passed in argument 'rdata' is not of class 'repro.data'!\n")
	
	#	load rjags
	require(rjags, quietly = TRUE)
	
	####### Variables
	#	variables names for the jags model
	Variables <- list(poisson = c("d", "log10b", "log10e"),	
			gammapoisson = c("d", "log10b","log10e", "log10omega"))	
	
	#	Data from repro.data
	
	rdata <- .repro.data2fit(rdata)
	
	#	Priors
	
	#	Priors distributions
	
	ddistr <- "dnorm(meand, taud)T(0,)"
	log10bdistr <- "dunif(log10bmin, log10bmax)"
	log10edistr <- "dnorm(meanlog10e , taulog10e)"
	log10omegadistr <- "dunif(log10omegamin, log10omegamax)"
	
	distr.prior <- list(ddistr = "normal",
			log10bdistr = "uniform",
			log10edistr = "normal",
			log10omegadistr = "uniform")
	
	#creation model.text
	
	#	poisson
	poisson.model.text <- .repro.loglogistic.poisson.text(.poisson.model.text.p1,
			ddistr,
			.poisson.model.text.p2,
			log10bdistr,
			.poisson.model.text.p3,
			log10edistr,
			.poisson.model.text.p4)
	
	#	gammapoisson
	gammapoisson.model.text <- .repro.loglogistic.gammapoisson.text(.gammapoisson.model.text.p1,
			ddistr,
			.gammapoisson.model.text.p2,
			log10bdistr,
			.gammapoisson.model.text.p3,
			log10edistr,
			.gammapoisson.model.text.p4,
			log10omegadistr,
			.gammapoisson.model.text.p5)
	
	#	Priors parameters
	
	paramdefault <- .repro.loglogistic.auto.param.prior(rdata$concmin,
			rdata$concmax,
			rdata$NreprocumulIndtime0,
			rdata$tab0)
	
	####### list of parameters values used to define Priors
	param.prior <- list(n = rdata$n,
			xconc = rdata$conc,
			Nindtime = rdata$Nindtime,
			meanlog10e = paramdefault$meanlog10e,					  
			taulog10e = paramdefault$taulog10e,					  
			meand = paramdefault$meand,
			taud = paramdefault$taud,
			yNcumul = rdata$Ncumul,
			log10omegamin = paramdefault$log10omegamin,
			log10omegamax = paramdefault$log10omegamax,								
			log10bmin = paramdefault$log10bmin,
			log10bmax = paramdefault$log10bmax)
	
	
	#	Model computing
	
	####### Define model
	poisson.model <- .repro.load.poisson.model(model.program = poisson.model.text,
			data = param.prior,
			n.chains,
			quiet)
	####### Determine sampling parameters
	poisson.sampling.parameters <- .repro.model.sampling.parameters(poisson.model,
			Variables$poisson, n.chains, quiet)
	
	if(poisson.sampling.parameters$niter > 100000){
		error <- 'The models need too many iterations to provide reliable parameter estimates'
		return(error)
	}
	
	####### Model Selection by the DIC
	
	poisson.DIC <- .repro.DIC(poisson.model, poisson.sampling.parameters, quiet)
	####### Define model		
	gammapoisson.model <- .repro.load.gammapoisson.model(model.program = gammapoisson.model.text,
			data = param.prior,
			n.chains,
			quiet)
	####### Determine sampling parameters		
	gammapoisson.sampling.parameters <- .repro.model.sampling.parameters(gammapoisson.model,
			Variables$gammapoisson, n.chains, quiet)
	
	####### parameters
	if(gammapoisson.sampling.parameters$niter > 100000){
		listparam <- .repro.listllparameters(poisson.model,
				"P",
				poisson.sampling.parameters,
				Variables$poisson,
				poisson.DIC)	
		
	}
	else {
		gammapoisson.DIC <- .repro.DIC(gammapoisson.model, gammapoisson.sampling.parameters,
				quiet)
		
		if(poisson.DIC <= gammapoisson.DIC) {
			listparam <- .repro.listllparameters(poisson.model,
					"P",						
					poisson.sampling.parameters,
					Variables$poisson,
					poisson.DIC)	
		}
		else {
			listparam <- .repro.listllparameters(gammapoisson.model,
					"GP",						
					gammapoisson.sampling.parameters,
					Variables$gammapoisson,
					gammapoisson.DIC)		
		}
	}
	
	####### Sampling
	prog.b <- ifelse(quiet == TRUE, "none", "text")
	mcmc <- coda.samples(listparam$model, 
			listparam$variables, 
			n.iter = listparam$niter, 
			thin = listparam$thin, 
			progress.bar = prog.b)
	
	####### summarize estime.par et CIs
	# calculate from the estimated parameters
	estim.par <- .repro.PARAMS(mcmc, listparam$model.label)#, det.part)
	
	####### ECx calculation  estimated ECx and their CIs 95%
	
	estim.ECx <- .repro.ECX(mcmc, c(5,10,20,50))
	
	#### output
	OUT <- list(DIC = listparam$DIC,
			estim.ECx = estim.ECx,
			estim.par = estim.par, 
			mcmc = mcmc,
			model = listparam$model,
			model.label = listparam$model.label,
			n.chains = length(mcmc),
			n.burnin = listparam$nburnin,
			n.iter = listparam$niter,
			n.thin = listparam$thin,
			param.prior = param.prior,
			raw.data = rdata$raw.data,
			transformed.data = rdata$transformed.data) 
	
	class(OUT) <- "repro.fit"
	return(OUT)
	
}
