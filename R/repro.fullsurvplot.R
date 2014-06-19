repro.fullsurvplot <-
		function(data,
				xlab,
				ylab,
				type = "generic",				 
				addlegend = TRUE)
# user function #
### plot of survival data: one subplot for each concentration, and one color
#	for each replicate
# INPUT: - type: type of plot generic, lattice or ggplot2
#		  - xlab
#		  - ylab
#		  - addlegend
{
	
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
	
	# default argument
	
	if(missing(xlab))
		xlab <- "Time"
	if(missing(ylab))
		ylab <- "Number of survivors"
	
#	Plot
	
#	generic option	
	if (type == "generic"){
		.repro.fullsurvplot.generic(data, 
				xlab,
				ylab,
				addlegend)
	}else{
		
		# lattice option
		if (type == "lattice"){
			
			require(lattice, quietly = TRUE)		
			
			.repro.fullsurvplot.l(data, 
					xlab,
					ylab,
					addlegend)
		}else{
			
			# ggplot2 option
			if(type == "ggplot"){
				
				require(ggplot2, quietly = TRUE)
				
				.repro.fullsurvplot.gg(data, 
						xlab,
						ylab,
						addlegend)
			}
		}	
	}
}
