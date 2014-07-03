repro.survplot <-
		function(data, 
				xlab,
				ylab,
				pch,
				type = "generic",				
				log.scale = FALSE,
				addlegend = TRUE,
				...)


# user function #
### plot of survival data
{
	# INPUTS
	# - data: raw dataframe with 5 columns with:
	#   - replicate: replicate indentification
	#   - conc: tested concentrations
	#   - time: time of the observation
	#   - Nsurv: number of alive individuals at time "time" and at concentration "conc"
	#   - Nrepro: number of collected offspring at concentration "conc" and at time "time"
	# - xlab, ylab
	# - type: generic or ggplot
	# - [log.scale: bool]: set log scale 
	# - addlegend: bool legend option
	# OUTPUT:
	# - Plot
	
	#	check data
	if (!is.data.frame(data)) 
		stop("data.frame expected")
	
	ref.names <- c("replicate","conc","time","Nsurv","Nrepro")
	test.names <- ref.names[which(is.na(match(ref.names,names(data))))]
	if (length(test.names) != 0)
		stop(paste("\nThe column",test.names,					
						"is missing or wrong name in the data.frame 'data' !\n",sep=" "))	
	
	rm(ref.names)
	rm(test.names)
	
	#	default legend
	if(missing(xlab))
		xlab <- "Concentrations"
	if(missing(ylab))
		ylab <- "Number of survivors"
	legend.position <- "right"
	
	#######
	
	timeend <- max(data$time)
	data <- data[data$time==timeend,]
	# select only no null concentration datapoint for log reprsentation 
	sel <- if(log.scale) data$conc > 0  else rep(TRUE, length(data$conc))
	if(log.scale){
		data$conc2[sel] <- log(data$conc[sel])
	}
	else{
		data$conc2[sel] <- data$conc[sel]
	}
	#	generic	
	if(type == "generic"){	

		#	calculate number of point by coordinate
		tt <- xyTable(na.omit(cbind(data$conc2, data$Nsurv)))
		
		# plot
		if(missing(pch))
			pch <- 16
		
		plot(tt$x,tt$y, 
				cex = (tt$number)/3, 
				xlab = xlab,
				ylab = ylab,
				pch = pch,
				xaxt = "n",
				yaxt = "n",
				...)
		
		# axis option
		axis(side = 2, at = .repro.keep.only.ints(pretty(c(0,max(data$Nsurv[sel])))))
		axis(side = 1, at = unique(data$conc2[sel]), labels = unique(data$conc[sel]))
		
		# legend option
		if(addlegend){
		legend(legend.position,
				pt.cex = sort(unique((tt$number)/3)), 
				pch = rep(pch, length(unique(tt$number))), 
				legend = paste(sort(unique(tt$number)), "replic.", sep = " "),
				bty = "n")
		}
	}
	
	
	#	ggplot
	if(type == "ggplot"){
		
		require(ggplot2, quietly = TRUE)
		
		conc = NULL
		Nsurv = NULL
		..n.. = NULL
		
		data <- data[sel,]	
		sp <- ggplot(data,aes(conc,Nsurv)) +
				stat_sum(aes(size = factor(..n..))) +
				labs(x = xlab,
						y = ylab)
		
		# legend option
		if(addlegend){
			if (log.scale) 
				fd <- sp + scale_x_log10(breaks = unique(data$conc)[2:length(unique(data$conc))],
								labels = unique(data$conc)[2:length(unique(data$conc))]) + 
						theme(legend.position = legend.position) +
						theme(legend.title = element_blank())
			else
				fd <- sp + scale_x_continuous(breaks = unique(data$conc)) + 
						theme(legend.position = legend.position) +
						theme(legend.title=element_blank())
		}
		else{
			if (log.scale) 
				fd <- sp + scale_x_log10(breaks = unique(data$conc)[2:length(unique(data$conc))],
								labels = unique(data$conc)[2:length(unique(data$conc))])  + 
						theme(legend.position = "none")
			else
				fd <- sp + scale_x_continuous(breaks = unique(data$conc)) + 
						theme(legend.position = "none")
		}	
		return(fd)
	}
}
