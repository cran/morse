repro.data <- 
		function(data, target.time)
# user function #
####### Function to define data4JAGS which will be the argument data in call to jags.model
{
	# INPUTS:
	#	-data: raw data frame
	#	-target.time: expected period of a part
	#			of the experience 
	#			(default duration of experience) 	
	
	# OUTPUT:
	# - rdata: list with data to run jags model:
	#	- raw.data: raw dataframe (Input)
	#	- transform.data: a dataframe with 6 columns corresponding to survival, Nindtime and repro data
	#					   at the time given in the input of the function
	#	- target.time
	
	#	function of data transformation
	
	#	check data
	if (!is.data.frame(data)) 
		stop("data.frame expected")
	
	#	dataset
	raw.data <- data
	
	#target.time default
	if(missing(target.time))
		target.time <- data$time[length(data$time)]
	
	# correct target time
	
	if(!any(data$time == target.time))
		stop("target.time is not one of the possible time !")
	
	#	transformed data
	transformed.data <- .repro.TransformData(data, target.time)
	
	#	rdata object
	rdata <- list(raw.data = raw.data,	
			transformed.data = na.omit(transformed.data),
			target.time = target.time)
	
	class(rdata) <- "repro.data"
	
	return(rdata)
	
}
