summary.repro.data <-
		function(object,...)
#	INPUT: object repro.doata object
{
	#	raw.data
	raw <- object$raw.data
	#	number of rows
	ans1 <- dim(raw)[1]
	#	table of datapoints per replicate
	ans2 <- table(raw$replicate)
	# table of datapoints per concentration
	ans3 <- table(raw$conc)
	# table of datapoints per time
	ans4 <- table(raw$time)
	
	#	transformed data
	transform <- object$transformed.data
	#	number of rows
	ans5 <- dim(transform)[1]
	#	table of datapoints per replicate 
	ans6 <- table(transform$replicate)
	# table of datapoints per concentration
	ans7 <- table(transform$conc)
	
	#	target.time
	target.time <- object$target.time
	
	cat("Summary: \n\n")
	cat("Raw data: \n\n")
	cat("Total number of datapoints: ", ans1, "\n")
	cat("\nNumber of datapoints per replicate: \n")
	print(ans2)	
	cat("\nNumber of datapoints per concentration: \n")
	print(ans3)
	cat("\nNumber of datapoints per time: \n")
	print(ans4)
	cat("\n\nTransformed data: \n")
	cat("\nNumber of the datapoints at target time",
			target.time,":", ans5, "\n\n")
	cat("\nNumber of datapoints per replicate: \n")
	print(ans6)	
	cat("\nNumber of datapoints per concentration: \n")
	print(ans7)
}
