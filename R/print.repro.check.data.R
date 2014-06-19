print.repro.check.data <-
		function(x,...)
#	print list of raw data, transformed data and control data
# INPUT: -x: repro.check.data object
# OUTPUT: - message of errors dataframe
{
	if(is.null(x$id))
		cat("The dataframe seems to be correct !\n")
	else{
		cat("Warning(s):\n")
		print(x$msg)
	}
}