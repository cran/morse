print.repro.fit <-
		function(x,...)

{
	model <- x$model
	#	M.C.M.C. informations
	ans1 <- t(data.frame("n.chain" = x$n.chain,
					"n.iter" = x$n.iter,
					"n.burnin" = x$n.burnin,
					"n.thin" = x$n.thin,
					"DIC" = x$DIC))
	colnames(ans1) <- "values"
	#	estimated parameters
	estim.par <- x$estim.par
	#	estimated ECX
	estim.ECx <- x$estim.ECx
	
	cat("Model:\n")
	print(model)
	cat("\nComputing information:\n\n")
	print(ans1)
}
