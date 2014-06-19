plot.repro.fit <-
		function(x,
				xlab,
				ylab,
				fitcol,
				fitlty,
				fitlwd,
				ci = FALSE,
				cicol,
				cilty,
				cilwd,
				addlegend = TRUE,
				log.scale = FALSE,
				type = "generic",
				...)	

#	user function	#		
{
	# INPUTS
	# - x:  repro.fit object
	# - xlab : label x
	# - ylab : label y
	# - fitcol : color fitted curve
	# - fitlty : type line fitted curve
	# - fitlwd : width line fitted curve
	# - ci : credible interval
	# - cicol : color ci
	# - cilty : type line ci
	# - cilwd : width line ci
	# - addlegend
	# - log.scale : x log option
	# - type : generic ou ggplot
	# OUTPUT:
	# - plot of fitted regression
	
	# library
	if(type == "ggplot")
		require(ggplot2, quietly = TRUE)
	
	## Define data
	concentrations <- x$transformed.data$conc
	response <- x$transformed.data$Nreprocumul / x$transformed.data$Nindtime
	
	
	## Fitted curve parameters
	X <- seq(min(concentrations), max(concentrations), length=100)
	res.M <- summary(x$mcmc)
	
	## Choose median of posteriors as parameter values
	
	d <- res.M$quantiles["d","50%"]
	b <- 10^res.M$quantiles["log10b","50%"]
	e <- 10^res.M$quantiles["log10e","50%"]
	fNcumulpidtheo <- d/(1+(X/e)^b)
	
	## IC 95
	if(ci){
		mcmctot <- do.call("rbind", x$mcmc)
		k <- nrow(mcmctot)
		#	parameters
		d2 <- mcmctot[,"d"]
		log10b2 <- mcmctot[,"log10b"]
		b2 <- 10^log10b2
		log10e2 <- mcmctot[,"log10e"]
		e2 <- 10^log10e2
		# quantile
		qinf95 = NULL
		qsup95 = NULL
		
		# poisson
		if(x$model.label == "P"){		
			for (i in 1:length(X))
			{
				theomean <- d2/(1 + (X[i] / e2)^(b2))
				
				# IC 95%
				qinf95[i] <- quantile(theomean, probs = 0.025, na.rm = TRUE)
				qsup95[i] <- quantile(theomean, probs = 0.975, na.rm = TRUE)				
			}		
		}		
		
		# gamma poisson
		if(x$model.label == "GP"){		
			# parameter		
			log10omega2 <- mcmctot[,"log10omega"]	
			omega2 <- 10^(log10omega2)	
			for (i in 1:length(X))
			{
				theomean <- d2 / (1 + (X[i] / e2)^(b2))		
				theo <- rgamma(n = k, shape = theomean / omega2, rate = 1 / omega2)
				
				# IC 95%
				qinf95[i] <- quantile(theo, probs = 0.025, na.rm = TRUE)
				qsup95[i] <- quantile(theo, probs = 0.975, na.rm = TRUE)
				
			}	
		}
	}
	
	## Define visual parameters
	mortality <- c(0,1)
	
	nomortality <- match(x$transformed.data$Nsurv == x$transformed.data$Ninit,
			c(TRUE, FALSE))	# valid if at least one replicat
# without mortality
	mortality <- mortality[nomortality] #vector of 0 and 1
	
# encodes mortality empty dots (1) and not mortality solid dots (19)
	if(type == "generic")
		mortality[which(mortality == 0)] <- 19
	if(type == "ggplot"){
		mortality[which(mortality == 0)] <- "No"
		mortality[which(mortality == 1)] <- "Yes"
	}
	
	# default axis parameters
	if(missing(xlab))
		xlab <- "Concentrations"
	if(missing(ylab))
		ylab <- "Response"
	
	# default legend parameters
	
	if(missing(fitcol))
		fitcol <- "red"
	if(missing(fitlty))
		fitlty <- 1
	if(missing(fitlwd))
		fitlwd <- 1
	if(type == "generic"){
		legend.position <- "bottomleft"
		legend.position.ci <- "left"
	}	
	legend.title <- "Mortality"
	legend.name.no <- "No"
	legend.name.yes <- "Yes"
	
	# IC parameters
	if(missing(cicol))
		cicol <- "red"
	if(missing(cilty))
		cilty <- 2
	if(missing(cilwd))
		cilwd <- 1
	
	## Plotting data
	if(type == "generic"){
		if(!ci){
			plot(concentrations, response,
					xlab = xlab,
					ylab = ylab,
					pch = mortality,
					xaxt = "n",
					yaxt = "n",
					log = if(log.scale) "x" else "",					
					...)
			axis(side = 2, at = .repro.keep.only.ints(pretty(c(0,max(response)))))
			axis(side = 1, at = unique(concentrations), labels = unique(concentrations))
			
			## Plotting the theoretical curve
			
			#fitted curve
			
			lines(X, fNcumulpidtheo, col = fitcol, 
					lty = fitlty, lwd = fitlwd, 
					type = "l")
		}
		if(ci){			
			# plotting data
			
			plot(concentrations, response,
					xlab = xlab,
					ylab = ylab,
					pch = mortality,
					xaxt = "n",
					yaxt = "n",
					ylim = c(0,max(qsup95)),
					log = if(log.scale) "x" else "",
					...)
			axis(side = 2, at = .repro.keep.only.ints(pretty(c(0,max(qsup95)))))
			axis(side = 1, at = unique(concentrations), labels = unique(concentrations))
			
			## Plotting the theoretical curve
			
			#fitted curve
			
			lines(X, fNcumulpidtheo, col = fitcol, 
					lty = fitlty, lwd = fitlwd, 
					type = "l")
			
			# C.I.
			
			lines(X, qsup95, type = "l", col = cicol, lty = cilty, lwd = cilwd)
			lines(X, qinf95, type = "l", col = cicol, lty = cilty, lwd = cilwd)	
		}
		
		## legend
		if(addlegend  && !ci){
			legend(legend.position, title = legend.title, pch = c(19, 1, NA),
					lty = c(0, 0, fitlty),
					lwd = c(1, 1, fitlwd),
					col = c(1, 1, fitcol),
					legend = c(legend.name.no, legend.name.yes, "Fitted curve"),
					bty = "n")
		}
		
		if(addlegend  && ci){
			legend(legend.position, title = legend.title, pch = c(19, 1, NA, NA),
					lty = c(0, 0, fitlty, cilty),
					lwd = c(1, 1, fitlwd, cilwd),
					col = c(1, 1, fitcol, cicol),
					legend = c(legend.name.no, legend.name.yes,
							"Fitted curve", "Credible limits"),
					bty = "n")
		}
	}
	if(type == "ggplot"){
		#	dataframes points and curve
		Line = NULL
		
		data.one <- data.frame(concentrations, 
				response, 
				mortality)
		
		data.two <- data.frame(X, 
				fNcumulpidtheo, 
				Line = "Fitted curve")
		
		# colors
		# points vector
		n <- 2
		cols <- hcl(h=seq(15, 375-360/n, length=n)%%360, c=100, l=65)
		cols1 <- cols[1:2]
		names(cols1) <- c("No", "Yes")
		# fitted curve
		cols2 <- fitcol
		names(cols2) <- c("Fitted curve")
		
		# points (to create the legend)
		plt_1 <- ggplot(data.one) +
				geom_point(data = data.one, aes(concentrations, response,
								color = mortality)) + scale_color_manual(values = cols1)
		
		# curve (to create the legend)
		plt_2 <- ggplot(data.one) +
				geom_line(data = data.two, aes(X, fNcumulpidtheo, color = Line),
						linetype = fitlty, size = fitlwd) +
				scale_color_manual(values = cols2)
		
		
		if(ci){
			
			# IC
			Ci = NULL
			
			data.three <- data.frame(X, qinf95, qsup95, 
					Ci = "Credible limits")  		
			# colors 
			cols3 <- cicol
			names(cols3) <- c("Credible limits")
			
			plt_3 <- ggplot(data.one) +
					geom_line(data = data.three, aes(X, qinf95, color = Ci),
							linetype = cilty, size = cilwd) + 
					geom_line(data = data.three, aes(X, qsup95, color = Ci),
							linetype = cilty, size = cilwd) +
					scale_color_manual(values = cols3)
			
			# final plot
			plt_4 <- ggplot(data.one) +
					geom_point(data=data.one, aes(concentrations, response, color = mortality)) +
					geom_line(aes(X, fNcumulpidtheo), data.two, linetype = fitlty, size = fitlwd,
							color = cols2) +
					geom_line(aes(X, qinf95), data.three, linetype = cilty, size = cilwd, color = cols3) +
					geom_line(aes(X, qsup95), data.three, linetype = cilty, size = cilwd, color = cols3) +
					scale_color_discrete(guide = "none") +
					labs(x = xlab,
							y = ylab)
		}
		
		if(!ci){
			plt_4 <- ggplot(data.one) +
					geom_point(data=data.one, aes(concentrations, response, color = mortality)) +
					geom_line(aes(X, fNcumulpidtheo), data.two, linetype = fitlty, size = fitlwd,
							color = cols2) +
					scale_color_discrete(guide = "none") +
					labs(x = xlab,
							y = ylab)					
		}
		
		
		
		if(addlegend){
			#	create legends
			mylegend_1 <- .repro.ggplotfit.legend(plt_1)
			mylegend_2 <- .repro.ggplotfit.legend(plt_2)
			if(ci) mylegend_3 <- .repro.ggplotfit.legend(plt_3)
			
			if (log.scale) 
				plt_5 <- plt_4 +
						scale_x_log10(breaks = unique(data.one$concentrations)[2:length(unique(data.one$concentrations))],
								labels = unique(data.one$concentrations)[2:length(unique(data.one$concentrations))])				
			else
				plt_5 <- plt_4 + scale_x_continuous(breaks = unique(data.one$concentrations))	
			
			
			library(gridExtra, quietly = TRUE)
			
			if(!ci)
				grid.arrange(plt_5, arrangeGrob(mylegend_1, mylegend_2, nrow = 6),
						ncol = 2, widths = c(7,1))
			if(ci)
				grid.arrange(plt_5, arrangeGrob(mylegend_1, mylegend_2, mylegend_3, nrow = 6),
						ncol = 2, widths = c(7,1))		
		}
		else{
			if(log.scale) 
				plt_5 <- plt_4 +
						scale_x_log10(breaks = unique(data.one$concentrations)[2:length(unique(data.one$concentrations))],
								labels = unique(data.one$concentrations)[2:length(unique(data.one$concentrations))])	
			else 
				plt_5 <- plt_4 + scale_x_continuous(breaks = unique(data.one$concentrations))	
			
			
			return(plt_5)
		}
		
	}
}
