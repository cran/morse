repro.check.data <-
		function(data, diagnos.plot = TRUE)	
# user function #
##### test data consistency
## Performs several tests on the data and returns at the first it
## encounters
{
	# INPUTS
	# - data: raw dataframe with 5 columns with:
	#   - replicate: replicate indentification
	#   - conc: tested concentrations
	#   - time: time of the observation
	#   - Nsurv: number of alive individuals at time "time" and at concentration "conc"
	#   - Nrepro: number of collected offspring at concentration "conc" and at time "time"
	#  !!!!!! data are supposed to be sorted by replicate, conc and time !!!!!!
	# - diagnos.plot: option to call the repro.fullsurvplot
	# OUTPUT
	# - dataframe with 2 columns:
	#   - id: the ID of the test
	# 	- msg: the message error with warning print
	
	#	function to generate rows of dataframe
	error <- function(id, msg) {
		data.frame(id = id, msg = msg, stringsAsFactors = FALSE)
	}
	#	errors dataframe
	errors <- data.frame()	
	
	### 1 ### test if the column names are correct
	
	ref.names <- c("replicate","conc","time","Nsurv","Nrepro")
	test.names <- ref.names[which(is.na(match(ref.names,names(data))))]
	
	if (length(test.names) != 0){
		err <- error("missingcolumn",
				paste("The column",test.names,"is missing or have a wrong name.",sep=" "))
		errors <- rbind(errors,err)
		
		class(errors) <- "repro.check.data"
		return(errors)
	}
	
	### storage the numbers of time
	lengthTime <- length(unique(data$time))
	
	### 2 ### test if the first time is zero
	
	subdata <- split(data, list(data$replicate, data$conc))
	subdata2 <- subdata[sapply(subdata,dim)[1,] != 0]
	if(any(lapply(subdata2,function(x) x$time[1] != 0) == TRUE)){
		err <- error("firstTime0",
				"Data are required at time 0 for each concentration and each replicate.")
		errors <- rbind(errors, err)   
	}
	
	### 3 ### test if concentrations are numeric
	
	if(!is.numeric(data$conc)){
		err <- error("concNumeric",
				"Column 'conc' must contain only numerical values.")
		errors <- rbind(errors, err)
	}
	
	### 4 ### test if Nsurv are integer							
	
	if(!is.integer(data$Nsurv)){
		err <- error("NsurvInteger",
				"Column 'Nsurv' must contain only integer values.")
		errors <- rbind(errors, err)
	}
	
	### 5 ### test if Nrepro are integer	
	
	if(!is.integer(data$Nrepro)){		
		err <- error("NreproInteger",
				"Column 'Nrepro' must contain only integer values.")
		errors <- rbind(errors, err)
	}
	
	###	6 ### positivity test table
	
	#	setting table
	table <- subset(data, select = -c(replicate))
	
	#	test
	if(any(table[table<0.0])){
		err <- error("tablePositive",
				"Data must contain only positive values.")
		errors <- rbind(errors, err)			 
	}
	
	### 7 ### test Nrepro = 0 at time 0
	
	# setting vector
	datatime0 <- data[data$time == 0,]
	# test if Nrepro > 0 at time 0
	if(any(datatime0$Nrepro > 0)){ 
		
		err <- error("Nrepro0T0",
				"Nrepro should be 0 at time 0 for each concentration and each replicate.")
		errors <- rbind(errors, err)
	}
	
	.consistency <- function(subdata)	
	### Function to be used on a subdataset corresponding to one replicate at one concentration.
	### This function checks:
	###   - if each replicate appears once and only once at each time
	###   - if Nsurv is never increasing with time
	###   - if at each time T for which Nsurv = 0, Nrepro = 0 at time T+1
	{
		# errors consitency dataframe 
		errors2 <- data.frame()	
		
		#### 8 ### test if each replicate appears once and only once at each conc and time		
		
		if(nrow(subdata) != length(unique(subdata$time)))
		{			
			err2 <- error("onlyReplicate",
					paste("Replicate ", 
							unique(subdata$replicate),
							" appears on different lines for the same time point and the same
							concentration ", 
							unique(subdata$conc), ".", sep=""))  				
			errors2 <- rbind(errors2, err2)
		}
		
		### 9 ### test if there is no lack of replicate at each conc and time
		
		if(length(subdata$replicate)!=lengthTime)
		{
			err2 <- error("missingReplicate",paste("Replicate ", 
							unique(subdata$replicate),
							" is missing for at least one time points at concentration ", 
							unique(subdata$conc), ".", sep=""))
			errors2 <- rbind(errors2, err2)
		}
		
		### 10 ### test Nsurv monotony
		
		nonmonotonous <- subdata$Nsurv[-length(subdata$Nsurv)] < subdata$Nsurv[-1]
		if(any(nonmonotonous))
		{
			err2 <- error("NsurvMonotone",
					paste("For replicate ", unique(subdata$replicate),
							" and concentration ", unique(subdata$conc),
							", Nsurv increases at some time points compared to the previous one.", sep=""))
			errors2 <- rbind(errors2, err2)
		}
		
		### 11 ### test Nsurv = 0 at time t => Nrepro > 0 at time t-1
		
		NsurvT <- subdata$Nsurv[-length(subdata$Nsurv)]
		NreproTplus1 <- subdata$Nrepro[-1]
		if(any(NreproTplus1[NsurvT==0]>0))
		{
			err2 <- error("Nsurvt0Nreprotp1P",
					paste("For replicate ", 
							unique(subdata$replicate),
							" and concentration ", unique(subdata$conc), 
							", there are some Nsurv = 0 followed by Nrepro > 0 at the next time point.", sep=""))
			errors2 <- rbind(errors2, err2)
		}
		return(errors2)			
	}
	
	res <- by(data, list(data$replicate, data$conc), .consistency)
	err <- c(do.call("rbind", res))
	
	if(length(err)!=0){
		errors <- rbind(errors,err)		
	}
	
	#print diagnos plot option
	if(length(err)!=0){
		if(diagnos.plot == TRUE && "NsurvMonotone" %in% err == TRUE){
			repro.fullsurvplot(data)
		}}	
		
	class(errors) <- "repro.check.data"
	return(errors)	
}
