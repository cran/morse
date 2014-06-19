repro.cumulplot <-
		function(data,
				xlab,
				ylab,
				type = "generic",
				log.scale = FALSE, 		 
				addlegend = TRUE,
				...)

# user function  #
### plot of cumulative reproduction data without taking into account mortality
{
	# INPUTS
	# - data: raw dataframe with 5 columns with:
	#   - replicate: replicate indentification
	#   - conc: tested concentrations
	#   - time: time of the observation
	#   - Nsurv: number of alive individuals at time "time" and at concentration "conc"
	#   - Nrepro: number of collected offspring at concentration "conc" and at time "time"
	#  !!!!!! data are supposed to be sorted by replicate, conc and time !!!!!!
	# - [log.scale: bool]: set log scale 
	# - type: generic or ggplot2
	# - addlegend: boolean plot legend option
	# - ... generic option
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
	
	conc = NULL
	
	# calculating the sum of Nrepro
	Nreprocumul <- by(data, data$conc, function(c)
			{
				by(c, c$replicate, 
						function(x){sum = sum(x["Nrepro"])})
			})
	# transformation vector
	Nreprocumul <- c(do.call("cbind", Nreprocumul))
	
	# calculating the mortality
	mortality <- by(data, data$conc, function(c)
			{
				by(c, c$replicate,
						function(x) x[x$time==0,]$Nsurv - x[x$time==max(x$time),]$Nsurv)
			})
	# transformation vector
	mortality <- c(do.call("cbind", mortality))
	# encodes mortality empty dots (1) and not mortality solid dots (19)
	mortality[mortality>0] <- 1 # mortality
	mortality[mortality==0] <- 19 #no mortality
	
	# plotting data
	
	if(missing(xlab))
		xlab <- "Concentrations"
	if(missing(ylab))
		ylab <- "Nreprocumul"
	
	#	default legend argument
	legend.position <- "right"
	legend.title <- "Mortality"
	legend.name.no <- "No"
	legend.name.yes <- "Yes"
	
	#generic
	
	if(type == "generic"){
		plot(sort(rep(unique(data$conc), length(unique(data$replicate)))), 
				Nreprocumul, 
				xlab = xlab,
				ylab = ylab,
				pch = mortality, 
				yaxt = "n",
				xaxt = "n",
				log = if(log.scale) "x" else "",
				...)
		axis(side = 2, at = .repro.keep.only.ints(pretty(c(0,max(Nreprocumul)))))
		axis(side = 1, at = unique(data$conc), labels = unique(data$conc))
		
		#	legend
		
		if(addlegend)
			legend(legend.position,title = legend.title,pch = c(19,1),
					bty="n", 
					legend=c(legend.name.no,legend.name.yes))
	}
	
	#ggplot2
	if(type == "ggplot"){
		
		require(ggplot2, quietly = TRUE)
		
		df <- data.frame(cbind(conc = sort(rep(unique(data$conc), 
										length(unique(data$replicate)))),
						Nreprocumul,mortality))
		
		# plot
		gp <- ggplot(df, aes(conc,Nreprocumul,colour = factor(mortality))) + 
				geom_point(size=2.5) +
				labs(x = xlab,
						y = ylab)
		
		# legend yes
		if(addlegend){
			if (log.scale) 
				fd <- gp +
						scale_x_log10(breaks = unique(data$conc)[2:length(unique(data$conc))],
								labels = unique(data$conc)[2:length(unique(data$conc))]) +
						scale_colour_hue(legend.title,
								breaks = c("19","1"),
								labels = c(legend.name.no, legend.name.yes)) +					
						theme(legend.position = legend.position)
			else
				fd <- gp + scale_x_continuous(breaks = unique(data$conc)) + 
						scale_colour_hue(legend.title,
								breaks = c("19","1"),
								labels = c(legend.name.no, legend.name.yes)) +
						theme(legend.position = legend.position)		
		}
		# legend no
		else{
			if (log.scale) 
				fd <- gp +
						scale_x_log10(breaks = unique(data$conc)[2:length(unique(data$conc))],
								labels = unique(data$conc)[2:length(unique(data$conc))]) + 
						theme(legend.position = "none")
			else
				fd <- gp + scale_x_continuous(breaks = unique(data$conc),
								labels = unique(data$conc)) + 
						theme(legend.position = "none")				
		}
		#return a ggplot object
		return (fd)	
	}
}
