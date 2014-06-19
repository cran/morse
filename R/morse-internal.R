.repro.keep.only.ints <-
#	function to select only integer for genric plot axis
		function(xs) xs[floor(xs) == xs]
.repro.calcul.Ninit <- 
		function(t1, t2, type)
#	calcul the correct number Ninit
#	for each replicate con and time
{
	#	INPUTS
	#	- x list of splited dataframe
	#	- tabletime0
	#	OUTPUTS
	#	- list of splited dataframe with the news column Ninit
	
	join(t1, t2, by = c("replicate","conc"), type = type)[,c("replicate", "conc",
					"time", "Nsurv", "Nrepro", "Ninit")]
}	
.repro.TransformData <-
		function(data, target.time) 

# internal function #
# original data to data for control data creation
{
	# INPUTS
	# - data: reordered original dataframe with 5 columns (original data):
	#   - replicate: replicate indentification
	#   - conc: tested concentration
	#   - time: observation time
	#   - Nsurv: number of alive individuals at concentration "conc" and at time "time"
	#   - Nrepro: number of collected offspring at concentration "conc" and at time "time"
	#  !!!!!! data are supposed to be sorted by replicate, conc and time !!!!!!
	# - time: the time we want to consider as the last time for the analysis. By default, it is the last time of original data.
	# OUTPUT: a dataframe with 6 columns corresponding to survival, Nindtime and repro data at the time given in the input of the function
	#   - replicate: replicate indentification
	#   - conc: tested concentration
	#   - Ninit: number of individuals at the beginning of the bioassay
	#   - Nsurv: number of alive individuals at the time specified in the input of the function
	#   - Nreprocumul: cumulative number of offspring at this time
	#   - Nindtime: number of individual-days from the beginning of the bioassay to this time
	
	temp <- split(data, data$time)
	### add column Ninit
	tabletime0 <- data[data$time==0,]
	tabletime0[,"Ninit"] <- tabletime0[,"Nsurv"]
	
	# if there are all replicate for all time and concentration
	if(length(which((sapply(temp,dim)[1,]) == sapply(temp,dim)[1,1])) == length(temp))
		res <- lapply(temp, function(x) .repro.calcul.Ninit(x, tabletime0, "right"))
	
	# if one or more replicate is missing	
	# NA gestion
	if(length(which((sapply(temp,dim)[1,]) < sapply(temp,dim)[1,1])) != 0){
		
		res <- lapply(temp, function(x) .repro.calcul.Ninit(x, tabletime0, "right"))
				
		# replace NA time by correct time value
		namenatime <- names(lapply(res,
						function(x){match(NA, as.matrix(x))})[lapply(res,
							function(x){match(NA, as.matrix(x))}) != "NA"])

		for(i in 1:length(namenatime)){
			res[namenatime][[namenatime[i]]]$time[which(is.na(res[namenatime][[namenatime[i]]]$time))] <- namenatime[i]
		
			# replace Ninit value by NA
			res[namenatime][[namenatime[i]]]$Ninit[which(is.na(res[namenatime][[namenatime[i]]]$Nsurv))] <- NA
		}
	}
	
	# if one or more replicate is added
#	if(length(which((sapply(temp,dim)[1,]) > sapply(temp,dim)[1,1])) != 0){
#
#		res <- lapply(temp, function(x) .repro.calcul.Ninit(x, tabletime0, "left"))
#		
#		# replace NA init by correct init value
#		namenainit <- names(lapply(res,
#						function(x){match(NA, as.matrix(x))})[lapply(res,
#								function(x){match(NA, as.matrix(x))}) != "NA"])
#		
#		for(i in 1:length(namenainit)){
#			res[namenainit][[namenainit[i]]]$Ninit[which(is.na(res[namenainit][[namenainit[i]]]$Ninit))] <- tabletime0[1, "Ninit"] 			
#		}
#	}
		
	res2 <- do.call(rbind, res)
	rownames(res2) <- 1:dim(res2)[1]
	res2$time <- as.numeric(res2$time)	
	### reorder dataset
	data <- res2[order(res2$replicate, res2$conc, res2$time),]

		
	T <-unique(data$time) # times of obs without repetitions
	finalnbr <- match(target.time, T) # index of the time at which we want the results in vector T
	if(finalnbr==1){stop("!!!! It isn't possible to use the first observation time as the last observation time !!!!")}
	
	#### calculation of cumulative number of repro and individual-days
	for(i in 2:finalnbr)
	{
		# original dataset at T[i-1]
		dataTim1 <- subset(data, data$time==T[i-1])
		# original dataset at T[i]
		dataTi <- subset(data, data$time==T[i])
		
		# check if data have been properly captured
		if(any(dataTim1$replicate!=dataTi$replicate) || any(dataTim1$conc!=dataTi$conc))
#		{stop("!!!! BE CAREFUL concentrations and/or replicates are not identical at each time !!!!")}
		{warning("!!!! BE CAREFUL concentrations and/or replicates are not identical at each time !!!!")}
		
		# build the output date at time T[i]
		if(i==2)
		{
			# (number of survivors at T[i]) * T[i]
			# + (number of dead organisms between T[i] and T[i-1]) * (T[i]+T[i-1])/2
			NindtimeTi <- ( dataTi$Nsurv*dataTi$time + 
						(dataTim1$Nsurv - dataTi$Nsurv)*(dataTi$time + dataTim1$time)/2 )
			
			NreprocumulTi <- dataTim1$Nrepro + dataTi$Nrepro
		}
		else
		{
			# (number of survivors at T[i]) * T[i]
			# + (number of dead organisms between t and tm1) * (t+tm1)/2
			# + number of individual-time at T[i-1] (accounting for dead organisms 
			# at time before T[i-1])
			# - number of individual-time at T[i-1] that are staying alive at T[i-1]
			NindtimeTi <- ( dataTi$Nsurv*dataTi$time + 
						(dataTim1$Nsurv - dataTi$Nsurv)*(dataTi$time + dataTim1$time)/2 +		
						tableTim1$Nindtime -						
						dataTim1$Nsurv * dataTim1$time )
			
			NreprocumulTi <- tableTim1$Nreprocumul + dataTi$Nrepro
		}
		
		tableTi <- data.frame(
				replicate = dataTi$replicate,
				conc = dataTi$conc,
				Ninit = dataTi$Ninit,
				Nsurv = dataTi$Nsurv,
				Nreprocumul = NreprocumulTi,  # cumulative number of offspring
				Nindtime = NindtimeTi)
		
		# tabelTi stored as tableTim1 for next iteration
		tableTim1 <- tableTi
	}
	
	tablefinale <- tableTi
	return(tablefinale)
}
.repro.data2fit <-
		function(rdata)
#	internal function	#
#	calcul data for repro.fit from repro.data data
#	INPUTS:
#	- rdata: repro.data object
#	OUTPUTS:
#	- raw.data
#	- transformed.data
#	- tab0: data at conc = 0
#	- tab: data at conc != 0
#	- Nindtime
#	- NreprocumulIndtime0
#	- conc
#	- Ncumul
#	- n
{
	####### Data from repro.data
	raw.data <- rdata$raw.data
	transformed.data <- rdata$transformed.data
	
	# control data
	tab0 <- rdata$transformed.data[rdata$transformed.data$conc == min(rdata$transformed.data$conc),]
	tab <- rdata$transformed.data[rdata$transformed.data$conc != min(rdata$transformed.data$conc),]
	
	Nindtime <- tab$Nindtime
	NreprocumulIndtime0 <- tab0$Nreprocumul/tab0$Nindtime
	conc <- tab$conc
	Ncumul <- tab$Nreprocumul
	n <- nrow(tab)
	
	## Parameter calculation of concentration min and max
	concmin <- min(sort(unique(conc))[-1])
	concmax <- max(conc)
	
	return(list(raw.data = raw.data,
					transformed.data = transformed.data, 
					tab = tab,
					tab0 = tab0,
					Nindtime = Nindtime,
					NreprocumulIndtime0 = NreprocumulIndtime0,
					conc = conc,
					Ncumul = Ncumul,
					n = n,
					concmin = concmin,
					concmax = concmax))
}
.repro.loglogistic.manu.distr.prior <-
		function(distr.prior){
	#	internal function	#
	#	create part of model.text with choosen distribution of priors variables
	#	INPUTS :
	#	- distr.prior : list of distribution for priors parameters
	#					loglogistic model
	#					- log10edistr : distribution of log10e prior "normal" or "uniform"
	#					- ddistr : distribution of d prior "normal" or "uniform"
	#					- log10bdistr : distribution of log10b prior "normal" or "uniform"
	#					- log10rdistr : distribution of log10r prior "normal" or "uniform"
	#	OUTPUTS :
	#	- distr.prior : list of characters chains for each part of model.text
	#					- ddistr
	#					- log10bdistr
	#					- log10edistr
	#					- log10rdistr
	
	#	distribution verification
	namedistr <- c("normal", "uniform")
	if(!any(namedistr == distr.prior[["ddistr"]]))
		stop("distribution of d not recognized !")
	
	if(!any(namedistr == distr.prior[["log10bdistr"]]))
		stop("distribution of log10b not recognized !")
	
	if(!any(namedistr == distr.prior[["log10edistr"]]))
		stop("distribution of log10e not recognized !")
	
	if(!is.null(distr.prior[["log10omegadistr"]]) && !any(namedistr == distr.prior[["log10omegadistr"]]))
		stop("distribution of log10omega not recognized !")	
	
	#	adapt part of JAGS model character string
	distr.prior2 = list()	
	
	distr.prior2[["ddistr"]] <- ifelse(distr.prior[["ddistr"]] == "normal",
			"dnorm(meand, taud)T(0,)",
			ifelse(distr.prior[["ddistr"]] == "uniform",
					"dunif(dmin, dmax)"))
	
	distr.prior2[["log10bdistr"]] <- ifelse(distr.prior[["log10bdistr"]] == "normal",
			"dnorm(meanlog10b, taulog10b)",
			ifelse(distr.prior[["log10bdistr"]] == "uniform",
					"dunif(log10bmin, log10bmax)"))
	
	distr.prior2[["log10edistr"]] <- ifelse(distr.prior[["log10edistr"]] == "normal",
			"dnorm(meanlog10e, taulog10e)",
			ifelse(distr.prior[["log10edistr"]] == "uniform",
					"dunif(log10emin, log10emax)"))
	
	distr.prior2[["log10omegadistr"]] <- ifelse(distr.prior[["log10omegadistr"]] == "normal",
			"dnorm(meanlog10omega, taulog10omega)",
			ifelse(distr.prior[["log10omegadistr"]] == "uniform",
					"dunif(log10omegamin, log10omegamax)"))		
	return(distr.prior2)	
}
.repro.loglogistic.poisson.text <-
		#	internal function	#
		#	create the poisson jags model text depending on the choosed distribution
		#	for each parameters
		#	INPUTS :
		#		- string part of the text
		#		- string part who define distribution in bugs language
		#	OUTPUTS :
		#		- complete bugs model in character string
		function(.poisson.model.text.p1, ddistr, .poisson.model.text.p2, log10bdistr,
				.poisson.model.text.p3,	log10edistr, .poisson.model.text.p4){
	return(paste(.poisson.model.text.p1,
					ddistr,
					.poisson.model.text.p2,
					log10bdistr,
					.poisson.model.text.p3,
					log10edistr,
					.poisson.model.text.p4, sep = ""))
}
.repro.loglogistic.gammapoisson.text <-
		#	internal function	#
		#	create the gammapoisson jags model text depending on the choosed distribution
		#	for each parameters
		#	INPUTS :
		#		- string part of the text
		#		- string part who define distribution in bugs language
		#	OUTPUTS :
		#		- complete bugs model in character string
		function(.gammapoisson.model.text.p1, ddistr, .gammapoisson.model.text.p2, log10bdistr,
				.gammapoisson.model.text.p3, log10edistr, .gammapoisson.model.text.p4,
				log10omegadistr, .gammapoisson.model.text.p5){
	return(paste(.gammapoisson.model.text.p1,
					ddistr,
					.gammapoisson.model.text.p2,
					log10bdistr,
					.gammapoisson.model.text.p3,
					log10edistr,
					.gammapoisson.model.text.p4,
					log10omegadistr,
					.gammapoisson.model.text.p5, sep = ""))
}
.repro.loglogistic.auto.param.prior <-
		function(concmin, concmax, NreprocumulIndtime0, tab0)
#	internal function	#
#	create default prior parameters for the log logistic model
{
	#	INPUTS:
	#	- Parameter calculation of concentration min and max
	#	- NreprocumulIndtime0: tab0$Nreprocumul/tab0$Nindtime
	#	- tab0: transformed data at concentration = 0 
	
	## Params to define log10EC50
	meanlog10e <- (log10(concmin)+log10(concmax))/2
	sdlog10e <- (log10(concmax)-log10(concmin))/4
	taulog10e <- 1/sdlog10e^2
	
	## Params to define d
	meand <- mean(NreprocumulIndtime0)
	SEd <- sd(NreprocumulIndtime0)/sqrt(length(unique(tab0$replicate)))
	taud <- 1/(SEd)^2
	
	## Params to define rate
	log10omegamin <- -4
	log10omegamax <- 4
	
	## Params to define b
	log10bmin <- -2
	log10bmax <- 2
	
	return(list(meanlog10e = meanlog10e,
					taulog10e = taulog10e,
					meand = meand,
					taud = taud,
					log10omegamin = log10omegamin,
					log10omegamax = log10omegamax,
					log10bmin = log10bmin,
					log10bmax = log10bmax))
}
.repro.loglogistic.manu.param.prior <-
		function(distr.prior, param.prior)
#	internal function	#
#	transform quantile to mean and tau for log10e and d
#	and quantile to boundary for b and r
{
	#	INPUTS:
	# - param.prior: list of:
	# 	- Q25log10e
	# 	- Q975log10e
	#	- Q25d
	#	- Q975d
	#	- Q975log10b
	#	- Q25log10b
	#	- Q975log10r
	#	- Q25log10r
	#	OUTPUTS :
	# - param.prior2 : list of :
	#	- meanlog10e
	#	- sdlog10e
	#	- meand
	#	- SEd
	#	- taulog10e
	#	- taud
	#	- log10bmax
	#	- log10bmin
	#	- log10rmax
	#	- log10rmin
	param.prior2 <- param.prior
	#	calcul of parameters for distribution
	#	e
	#	normal
	if(distr.prior[["log10edistr"]] == "normal"){
		param.prior2[["meanlog10e"]] <- (param.prior2[["Q25log10e"]] +
					param.prior2[["Q975log10e"]])/2
		
		param.prior2[["sdlog10e"]] <- (param.prior2[["Q975log10e"]] -
					param.prior2[["Q25log10e"]])/
				(2 * qnorm(0.975))
		
		#	calcul of the precision log10e
		param.prior2[["taulog10e"]] <- 1/(param.prior2[["sdlog10e"]])^2
		param.prior2[["sdlog10e"]] <- NULL
	}
	#	uniform
	if(distr.prior[["log10edistr"]] == "uniform"){
		param.prior2[["log10emax"]] <- ((param.prior2[["Q975log10e"]] +
						param.prior2[["Q25log10e"]]) +
					(param.prior2[["Q975log10e"]] -
						param.prior2[["Q25log10e"]])/0.95)/2
		
		param.prior2[["log10emin"]] <- (param.prior2[["Q975log10e"]] +
					param.prior2[["Q25log10e"]]) -
				param.prior2[["log10emax"]]
	}
	
	#	d
	#	normal
	if(distr.prior[["ddistr"]] == "normal"){
		param.prior2[["meand"]] <- (param.prior2[["Q25d"]] +
					param.prior2[["Q975d"]])/2
		
		param.prior2[["SEd"]] <- (param.prior2[["Q975d"]] -
					param.prior2[["Q25d"]])/
				(2 * qnorm(0.975))
		
		#	calcul of the precision d
		param.prior2[["taud"]] <- 1/(param.prior2[["SEd"]])^2
		param.prior2[["SEd"]] <- NULL
	}
	#	uniform
	if(distr.prior[["ddistr"]] == "uniform"){
		param.prior2[["dmax"]] <- ((param.prior2[["Q975d"]] +
						param.prior2[["Q25d"]]) +
					(param.prior2[["Q975d"]] -
						param.prior2[["Q25d"]])/0.95)/2
		
		param.prior2[["dmin"]] <- (param.prior2[["Q975d"]] +
					param.prior2[["Q25d"]]) -
				param.prior2[["dmax"]]
	}
	
	#	b
	#	normal
	if(distr.prior[["log10bdistr"]] == "normal"){
		param.prior2[["meanlog10b"]] <- (param.prior2[["Q25log10b"]] +
					param.prior2[["Q975log10b"]])/2
		
		param.prior2[["sdlog10b"]] <- (param.prior2[["Q975log10b"]] -
					param.prior2[["Q25log10b"]])/
				(2 * qnorm(0.975))
		
		#	calcul of the precision log10b
		param.prior2[["taulog10b"]] <- 1/(param.prior2[["sdlog10b"]])^2
		param.prior2[["sdlog10b"]] <- NULL
	}
	#	uniform
	if(distr.prior[["log10bdistr"]] == "uniform"){
		param.prior2[["log10bmax"]] <- ((param.prior2[["Q975log10b"]] +
						param.prior2[["Q25log10b"]]) +
					(param.prior2[["Q975log10b"]] -
						param.prior2[["Q25log10b"]])/0.95)/2
		
		param.prior2[["log10bmin"]] <- (param.prior2[["Q975log10b"]] +
					param.prior2[["Q25log10b"]]) -
				param.prior2[["log10bmax"]]
	}
	
	#	rate r
	#	omega
	if(!is.null(distr.prior[["log10omegadistr"]])){
		#	normal
		if(distr.prior[["log10omegadistr"]] == "normal"){
			param.prior2[["meanlog10omega"]] <- (param.prior2[["Q25log10omega"]] +
						param.prior2[["Q975log10omega"]])/2
			
			param.prior2[["sdlog10omega"]] <- (param.prior2[["Q975log10omega"]] -
						param.prior2[["Q25log10omega"]])/
					(2 * qnorm(0.975))
			
			#	calcul of the precision log10r
			param.prior2[["taulog10omega"]] <- 1/(param.prior2[["sdlog10omega"]])^2
			param.prior2[["sdlog10omega"]] <- NULL
		}
		#	uniform
		if(distr.prior[["log10omegadistr"]] == "uniform"){
			param.prior2[["log10omegamax"]] <- ((param.prior2[["Q975log10omega"]] +
							param.prior2[["Q25log10omega"]]) +
						(param.prior2[["Q975log10omega"]] -
							param.prior2[["Q25log10omega"]])/0.95)/2
			
			param.prior2[["log10rmin"]] <- (param.prior2[["Q975log10omega"]] +
						param.prior2[["Q25log10omega"]]) -
					param.prior2[["log10omegamax"]]
		}
	}
	#	clean list
	param.prior2[c("Q25log10e", "Q975log10e", "Q25d", "Q975d",
					"Q25log10b", "Q975log10b", 
					"Q25log10omega", "Q975log10omega")] <- NULL
	return(param.prior2)
}
.repro.load.model <-
		function(model.program, lr.bound.keep, data, n.chains, Nadapt = 3000, quiet = quiet)
#	internal function	#
#	create the JAGS model object and called by .repro.load.poisson.model
#	and .repro.load.gammapoisson.model
{
	# INPUTS:
	# - model.program: character string with a jags model description
	# - data: list with data to run jags model with 10 items:
	#   - n: number of lines of the tranformed dataframe
	#   - x: concentration vector. size concentration*replicat
	#   - Nindtime: vector of the number of individual-time
	#   - meanlog10EC50: log mean of data to define prior of log10EC50
	#   - taulog10EC50: log accuracy of data to define prior of log10EC50
	#   - meand: mean of Nreprocumul/Nindtime for control to define prior of d
	#   - taud: accuracy of Nreprocumul/Nindtime for control to define prior of d
	#   - y: vector of the number of collected data size concentration*replicate
	#   - lr.bound: bound to define prior of log10rate
	#   - lb.bound: bound to define prior of log10b
	# - Nchains: Nomber of chains desired
	# - Nadapt: length of the adaptation phase
	# - quiet: silent option
	# OUTPUT:
	# - JAGS model
	
	# delisting of lr.bound because not used in the function
	if(!lr.bound.keep) data[c("meanlog10omega", "taulog10omega", "log10omegamin", "log10omegamax")] <- NULL
	model.file <- tempfile()
	fileC <- file(model.file)
	writeLines(model.program, fileC)
	close(fileC)
	#	creation of the jags model
	model <- jags.model(file = model.file, data = data, n.chains = n.chains,
			n.adapt = Nadapt, quiet = quiet)
	unlink(model.file)
	return(model)
}
.repro.load.model.par <-
		function(cl, model.program, lr.bound.keep, data, n.chains, Nadapt = 3000, quiet = quiet)
#	internal function	#
#	create the JAGS model object and called by .repro.load.poisson.model
#	and .repro.load.gammapoisson.model
{
	
	# delisting of lr.bound because not used in the function
	if(!lr.bound.keep) data[c("meanlog10omega", "taulog10omega", "log10omegamin", "log10omegamax")] <- NULL
	model.file <- tempfile()
	fileC <- file(model.file)
	writeLines(model.program, fileC)
	close(fileC)
	#	creation of the jags model
	model <- parJagsModel(cl, name = "res", file = model.file,
			data = data, n.chains = n.chains,
			n.adapt = Nadapt, quiet = quiet)
	unlink(model.file)
	return(model)
}
.repro.load.poisson.model <-
		# sub function to load jags poisson model
		function(model.program, data, n.chains, quiet = quiet) .repro.load.model(model.program, F,
			data, n.chains, quiet = quiet)
.repro.load.poisson.model.par <-
		# sub function to load jags poisson model
		function(cl, model.program, data, n.chains, quiet = quiet) .repro.load.model.par(cl, model.program, F,
			data, n.chains, quiet = quiet)
.repro.load.gammapoisson.model <-
		# sub function to load jags gamma poisson model
		function(model.program, data, n.chains, quiet = quiet) .repro.load.model(model.program, T, 
			data, n.chains, quiet = quiet)
.repro.load.gammapoisson.model.par <-
		# sub function to load jags gamma poisson model
		function(cl, model.program, data, n.chains, quiet = quiet) .repro.load.model.par(cl, model.program, T, 
			data, n.chains, quiet = quiet)
.repro.model.sampling.parameters <-
		function(model, variables, n.chains, quiet = quiet)
#	internal function	
{
	# INPUTS: 
	# - jags model and variables from loading function
	# OUTPUTS:
	# - niter: number of iteration (mcmc)
	# - thin: thining rate parameter
	# - burnin: number of iteration burned
	
	# number of iteration for the pilote run required by raftery.diag
	# default value: 3746
	niter.init <- 5000
	prog.b <- ifelse(quiet == TRUE, "none", "text")
	mcmc <- coda.samples(model, variables, n.iter=niter.init, thin=1, progress.bar = prog.b)
	RL <- raftery.diag(mcmc)
	
	# check raftery.diag result (number of sample for the diagnostic procedure)
	
	if(n.chains < 2)
		stop('2 or more parallel chains required')
	else{
		resmatrix <- RL[[1]]$resmatrix
		for(i in 2: length(RL)){
			resmatrix <- rbind(resmatrix, RL[[i]]$resmatrix)
		}
	}
	#	creation of sampling parameters
	thin <- round(max(resmatrix[,"I"])+0.5)
	niter <- max(resmatrix[,"Nmin"]) * thin
	burnin <- max(resmatrix[,"M"])
	
	return(list(niter = niter, thin = thin, burnin = burnin))
}
.repro.model.sampling.parameters.par <-
		function(cl, model, variables, n.chains, quiet = quiet)
#	internal function	
{
	# INPUTS: 
	# - jags model and variables from loading function
	# OUTPUTS:
	# - niter: number of iteration (mcmc)
	# - thin: thining rate parameter
	# - burnin: number of iteration burned
	
	# number of iteration for the pilote run required by raftery.diag
	# default value: 3746
	niter.init <- 5000
	prog.b <- ifelse(quiet == TRUE, "none", "text")
	mcmc <- parCodaSamples(cl,"res", variables, n.iter=niter.init, thin=1, progress.bar = prog.b)
	RL <- raftery.diag(mcmc)
	
	# check raftery.diag result (number of sample for the diagnostic procedure)
	
	if(n.chains < 2)
		stop('2 or more parallel chains required')
	else{
		resmatrix <- RL[[1]]$resmatrix
		for(i in 2: length(RL)){
			resmatrix <- rbind(resmatrix, RL[[i]]$resmatrix)
		}
	}
	#	creation of sampling parameters
	thin <- round(max(resmatrix[,"I"])+0.5)
	niter <- max(resmatrix[,"Nmin"]) * thin
	burnin <- max(resmatrix[,"M"])
	
	return(list(niter = niter, thin = thin, burnin = burnin))
}
.repro.listllparameters <- 
		function(model, model.label, sampling.parameters, variables, DIC)
{
#	internal	#
#	for rpero.fit log logistic
#	create a list with
#	- model
#	- model.label
#	- niter
#	- thin
#	- nburnin
#	- variables
#	- DIC
	return(list(model = model,
					model.label = model.label,
					niter = sampling.parameters$niter,
					thin = sampling.parameters$thin,
					nburnin = sampling.parameters$burnin,
					variables = variables,
					DIC = DIC))
}	
.repro.PARAMS <-
		function(mcmc, MODEL="P") #, det.part)
{
	# INPUT:
	# - mcmc:  list of estimated parameters for the model with each item representing a chain
	# - MODEL: a position flag model with P: poisson model and GP: gammapoisson model
	# OUTPUT:
	# - data frame with 3 columns (values, CIinf, CIsup) and 3-4rows (the estimated parameters)
	
	## Retrieving parameters of the model
	res.M <- summary(mcmc)
	
#	if(det.part == "loglogistic"){
	d <- res.M$quantiles["d", "50%"]
	b <- 10^res.M$quantiles["log10b","50%"]
	e <- 10^res.M$quantiles["log10e","50%"]	
	dinf <- res.M$quantiles["d", "2.5%"]
	binf <- 10^res.M$quantiles["log10b","2.5%"]
	einf <- 10^res.M$quantiles["log10e","2.5%"]	
	dsup <- res.M$quantiles["d", "97.5%"]
	bsup <- 10^res.M$quantiles["log10b","97.5%"]
	esup <- 10^res.M$quantiles["log10e","97.5%"]	
#	}
	## Definition of the variable storage and storage data
	## If Poisson Model
#	if(det.part == "loglogistic" && MODEL=="P")
	if(MODEL == "P")
	{
		rownames = c("d", "b", "e")
		params = c(d, b, e)		
		CIinf = c(dinf, binf, einf)		
		CIsup = c(dsup, bsup, esup)		
	}
	## If Gamma Poisson Model
#	if(det.part == "loglogistic" && MODEL=="GP")
	if(MODEL == "GP")			
	{
#		## Calculation of the parameter rate
		## Calculation of the parameter omega
		omega <- 10^res.M$quantiles["log10omega","50%"]
		omegainf <- 10^res.M$quantiles["log10omega","2.5%"]
		omegasup <- 10^res.M$quantiles["log10omega","97.5%"]
		## Definition of the variable storage and storage data
		rownames = c("d", "b", "e", "omega")		
		params = c(d, b, e, omega)
		CIinf = c(dinf, binf, einf, omegainf)		
		CIsup = c(dsup, bsup, esup, omegasup)
	}
	res <- data.frame(median = params, Q2.5 = CIinf, Q97.5 = CIsup, row.names = rownames)
	return(res)
}
.repro.ECX <-
		function(mcmc, x)
{
	# INPUT:
	# - mcmc:  list of estimated parameters for the model with each item representing a chains
	# - x: vector of values of ECx
	# OUTPUT:
	# - data frame with the estimated ECx and their CIs 95% (3 columns (values, CIinf, CIsup) and 
	#   length(x) rows)
	
	# Retrieving estimated parameters of the model
	
	mctot <- mcmc[[1]]
	for(i in 2:length(mcmc)){
		mctot <- rbind(mctot,mcmc[[i]])
	}
	
	b <- 10^mctot[,"log10b"]
	e <- 10^mctot[,"log10e"]
	
	# Calculation ECx
	ECx <-sapply(x, function(x){e*((100/(100-x))-1)^(1/b)})
	
	q50 <- apply(ECx, 2, function(ECx){quantile(ECx, probs = 0.5)})
	qinf95 <- apply(ECx, 2, function(ECx){quantile(ECx, probs = 0.025)})
	qsup95 <- apply(ECx, 2, function(ECx){quantile(ECx, probs = 0.975)})
	
	# defining names
	ECname <- sapply(x, function(x){paste("EC", x, sep='')})
	
	res <- data.frame(median = q50, Q2.5 = qinf95, Q97.5 = qsup95, row.names = ECname)
	
	return(res)
}
.repro.DIC <-
		function(m.M, sampling.parameters, quiet = quiet)
#	internal function	#
####### DIc calcualtion
{
	# INPUTS
	# - m.M:  jags model object
	# - niter: number of iterations for the sampling
	# - thin
	# OUTPUT:
	# - numeric value of the DIC
	
	prog.b <- ifelse(quiet == TRUE, "none", "text")
	
	dic <- dic.samples(m.M, n.iter=sampling.parameters$niter, thin=sampling.parameters$thin, progress.bar = prog.b)
	return(round(sum(sapply(dic$deviance, mean) + sapply(dic$penalty, mean))))
}
.repro.fullsurvplot.generic <-
		function(data, 
				xlab,
				ylab,
				addlegend)

# internal function #
### plot of survival data: one subplot for each concentration, and one color for each replicate
{
	# INPUTS
	# - data: raw dataframe with 5 columns with:
	#   - replicate: replicate indentification
	#   - conc: tested concentrations
	#   - time: time of the observation
	#   - Nsurv: number of alive individuals at time "time" and at concentration "conc"
	#   - Nrepro: number of collected offspring at concentration "conc" and at time "time"
	# OUTPUT:
	# - Plot
	
	.convert <- function(x)
	## conversion of a replicate name in a number coding for color
	{
		# INPUT
		# - x: name of replicate
		# OUTPUT
		#  !!!!!! data are supposed to be sorted by replicate, conc and time !!!!!!
		# - a position of replicate in vector to display color
		replicate <- unique(data$replicate)
		mat <- matrix(ncol=2, nrow=length(replicate))
		mat[,1]=as.character(replicate)
		mat[,2]=as.character(seq(1, length(replicate)))
		return(as.integer(mat[mat[,1]==as.character(x)][2]))
	}
	
	.NbPlot <- function(conc)
	## definition of the number of subplots
	{
		# INPUT
		# - conc: vector of tested concentrations
		# OUTPUT
		# - vector defining the number of columns and rows in par(mfrow)
		
		nbconc <- length(conc)
		PlotPar <- c(c(2,2), c(2,3), c(2,4), c(3,3), c(2,5), c(3,4), c(3,5), c(4,4))
		NbPlotTheo <- matrix(ncol=2, nrow=8)
		NbPlotTheo[,1] <- c(1, 3, 5, 7, 9, 11, 13, 15)
		NbPlotTheo[,2] <- c(4, 6, 8, 9, 10, 12, 15, 16)
		if(nbconc<15){
			PlotTrue <- NbPlotTheo[NbPlotTheo[,2]-nbconc>0,][1,1]
		}
		else{PlotTrue=15}
		return(c(PlotPar[PlotTrue], PlotPar[PlotTrue+1]))
	}
	
	
	## creation of a vector of colors
	colors=rainbow(length(unique(data$replicate)))
	
	## split of the graphical window in subplots
	par(mfrow=.NbPlot(unique(data$conc)))
	
	by(data, data$conc, function(x){
				plot(x$time, rep(0, length(x$time)), 
						xlab = xlab,
						ylab = ylab, 
						ylim = c(0,max(x$Nsurv)), 
						type = "p", 
						col = 'white',
						yaxt = 'n')
				axis(side = 2, at = .repro.keep.only.ints(pretty(c(0,max(x$Nsurv)))))
				by(x, x$replicate, function(y){
							lines(y$time, y$Nsurv, type="b", 
									col = colors[.convert(unique(y$replicate))])
						})
				title(paste("Conc ", unique(x$conc), sep=""))
			})
	
	if(addlegend){		
		### creation of an empty plot to display legend
		plot(0, 0, 
				xlab = "", 
				ylab = "", 
				xlim = c(0,5),
				ylim = c(0,5),
				type = "n",
				xaxt = "n",
				yaxt = "n",
				bty = "n")
		
		### Display legend
		title.legend <- "Replicate"
		mat <- matrix(nrow = length(unique(data$replicate)), ncol = 2)
		mat[,1] <- rep(title.legend, length(unique(data$replicate)))
		mat[,2] <- unique(as.character(data$replicate))
		name <- apply(mat, 1, function(x){paste(x[1], x[2], sep = " ")})
		
		legend(0,5, name, 
				lty = rep(1,length(unique(data$replicate))), 
				pch = 21,
				col = colors, 
				bty = "n", 
				cex=1)
	}
	par(mfrow = c(1,1))
}
.repro.fullsurvplot.l <-
		function(data, 
				xlab,
				ylab,
				addlegend)

{
# INPUTS
# - data: raw dataframe with 5 columns with:
#   - replicate: replicate indentification
#   - conc: tested concentrations
#   - time: time of the observation
#   - Nsurv: number of alive individuals at time "time" and at concentration "conc"
#   - Nrepro: number of collected offspring at concentration "conc" and at time "time"
# OUTPUT:
# - Plot
	if(addlegend){
		xyplot(Nsurv ~ time | factor(conc), 
				data = data,
				group = replicate,
				type = "b",
				pch = 16,
				xlab = xlab,
				ylab = ylab,
				auto.key = list(space = "right",
						title = "Replicate",
						cex.title = 1,
						type = "b",
						pch = 16,
						lines = TRUE,
						points = FALSE))
	}
	else{
		xyplot(Nsurv ~ time | factor(conc), 
				data = data,
				group = replicate,
				type = "b",
				pch = 16,
				xlab = xlab,
				ylab = ylab,
				auto.key = list(draw = FALSE))
	}
}

.repro.fullsurvplot.gg <-
		function(data, 
				xlab,
				ylab,
				addlegend)
{
# INPUTS
# - data: raw dataframe with 5 columns with:
#   - replicate: replicate indentification
#   - conc: tested concentrations
#   - time: time of the observation
#   - Nsurv: number of alive individuals at time "time" and at concentration "conc"
#   - Nrepro: number of collected offspring at concentration "conc" and at time "time"
# OUTPUT:
# - Plot
	time = NULL
	Nsurv = NULL
	title.legend <- "Replicate"
	fg <- ggplot(data,
					aes(time,Nsurv,colour = factor(replicate))) +
			geom_point() + 
			geom_line() +
			labs(x = xlab, y = ylab) +
			facet_wrap(~conc,nrow = 2) +
			scale_x_continuous(breaks = unique(data$time))
	
	# legend option
	if(addlegend)
		fd <- fg + scale_colour_hue(title.legend)
	else{
		fd <- fg + theme(legend.position = "none")
	}
	return(fd)	
}
.repro.ggplotfit.legend <-
		function(a.gplot){
	#	create an independent legend for each ggplot object
	# INPUT
	# - a.gplot: ggplot object
	# OUPUT
	# - legend
	tmp <- ggplot_gtable(ggplot_build(a.gplot))
	leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
	legend <- tmp$grobs[[leg]]
	return(legend)
}
.poisson.model.text.p1 <-
		"model # Poisson model\n\t\t{\n\t\t#\t\n\t\tfor (j in 1:n) # loop on replicates\n\t\t{\n\t\t\t# Explicit writting of a Poisson law for each replicate\n\t\t\t# mean is given by the theoretical curve\n\t\t\tytheo[j] <- d / (1 + pow(xconc[j]/e, b))\n\t\t\tnbtheo[j] <- ytheo[j]*Nindtime[j]\n\t\t\tyNcumul[j] ~ dpois(nbtheo[j])\n\t\t}\n\t\t# Prior distributions\n\t\td ~ "
.poisson.model.text.p2 <- 
		"\n\t\tlog10b ~ "
.poisson.model.text.p3 <-
		"\n\t\tlog10e ~ "
.poisson.model.text.p4 <-
		"\n\t\t\n\t\tb <- pow(10,log10b)\n\t\te <- pow(10,log10e)\n\t\t}"

.gammapoisson.model.text.p1 <-
		"\n\t\tmodel #Gamma poisson model\n\t\t{\n\t\t#\t\n\t\tfor (j in 1:n) # loop on replicates\n\t\t{\n\t\t\t# Explicit writting of a gamma-Poisson law for each replicate\n\t\t\t# the mean is given by a gamma law centered on the theoretical curve\n\t\t\trate[j] <- d / (1 + pow(xconc[j]/e, b)) / omega\n\t\t\tp[j] <- 1 / (Nindtime[j] * omega + 1)\n\t\t\tyNcumul[j] ~ dnegbin(p[j], rate[j])\n\t\t}\n\t\t# Prior distributions\n\t\td ~ "
.gammapoisson.model.text.p2 <-
		"\n\t\tlog10b ~ "
.gammapoisson.model.text.p3 <-
		"\n\t\tlog10e ~ "
.gammapoisson.model.text.p4 <-
		"\n\t\tlog10omega ~ "
.gammapoisson.model.text.p5 <-
		"\n\t\t\t\t\n\t\tomega <- pow(10,log10omega)\n\t\tb <- pow(10,log10b)\n\t\te <- pow(10,log10e)\n\t\t}"
