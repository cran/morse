print.repro.data <-
		function(x,...)
#	print list of raw data, transformed data and control data
# INPUT: -x: repro.data object
{
	transformed.data <- x$transformed.data
	target.time <- x$target.time
	
	#output
	cat("\n\nTransformed data:\n")
	print(transformed.data)
	cat("\n\nTarget time:\n")
	print(target.time)	
}
