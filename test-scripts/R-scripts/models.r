buildDNAmodel <- function(inp, avoid) {

	modelName <- inp[[1]]
	free_params <- inp[[2]]
	num_ti_rates <- inp[[3]]
	num_tv_rates <- inp[[4]]
	equal_base_frequencies <- inp[[5]]
	rate_variation <- inp[[6]]
	substitution_scheme <- inp[[7]]
	model_titv <- inp[[8]]
	model_ratematrix <- inp[[9]]
	model_partition  <- inp[[10]]
	model_num_partitions <- inp[[11]]

	model <- c("",0,0,0,0,0,0,0,0,0,0,0,0,0,0)

	repeat {
		modelid <- sample(1:22,1,replace=T)
		if (!modelid %in% avoid) break
	}

	f_model <- equal_base_frequencies[modelid]
	g_model <- (rate_variation[modelid]==2 || rate_variation[modelid]==3)
	i_model <- (rate_variation[modelid]==1 || rate_variation[modelid]==3)

	model[1] <- modelName[modelid]
	if (f_model == 1) {
		model[2:5] <- 0.25
	} else {
	  base_frequencies <- numeric(4)
	  while (min(base_frequencies) < 0.1) {
	    base_frequencies <- rdirichlet(1, c(1,1,1,1) )
	  }
	  model[2:5] <- formatC(base_frequencies,digits=4)
    	}
    	
    	if (model_titv[modelid] == 1) {
	    if (modelid > 8) {
		#Transition/Transversion rate from a Gamma (2,1) truncated between 2 and 10
		titv <- rtrunc(1,"gamma", 2, 1, a=2, b=10)
		if (f_model == 1) {
		    kappa <- titv*2
		} else {
		    kappa <- titv*(base_frequencies[1]+base_frequencies[3])*(base_frequencies[2]+base_frequencies[4])/(base_frequencies[1]*base_frequencies[3] + base_frequencies[2]*base_frequencies[4])
		}
	    } else {
	  	titv <- 0.5
		kappa <- 1
	    }
	    model[6] <- formatC(kappa,digits=4)
	}
		
	if (model_ratematrix[modelid] == 1) {
	    partition <- substring(model_partition[modelid], seq(1,6,1), seq(1,6,1))

	    #R-matrix parameters from a Dirichlet(1.0,6.0,1.0,1.0,6.0,1.0) scaled with the last rate. The expected ti/tv is 12/4 = 3
	    validRates=0
	    while (!validRates) {
	    	switch(model_num_partitions[modelid],
		    {rates <- rdirichlet(1, c(1))},
		    {rates <- rdirichlet(1, c(1,6))},
		    {
		  	if (as.numeric(partition[5]) == 2) {
		   	    ti2=20
		  	} else {
		    	    ti2=2
		  	}
			rates <- rdirichlet(1, c(6,16,ti2))
		    },
		    {rates <- rdirichlet(1, c(6,16,2,8))},
		    {rates <- rdirichlet(1, c(6,16,2,8,4))},
		    {rates <- rdirichlet(1, c(6,16,2,8,20,4))}
	    	)

	    	rmatrix <- 0

		# build rate matrix
		for (m in 1:6) {
		    rmatrix[m] = rates[as.numeric(partition[m])+1]
		}
		
		#scaled with the last rate
		rmatrix_scaled <- rmatrix/rmatrix[6]
		validRates=(sum(rmatrix_scaled<0.1) == 0) && (sum(rmatrix_scaled>100) == 0)
		if (validRates) {
		    for (ind1 in 2:6) {
		    	for (ind2 in 1:(ind1-1)) {
			    if (rmatrix_scaled[ind1] != rmatrix_scaled[ind2]) {
			        if (abs(rmatrix_scaled[ind1]-rmatrix_scaled[ind2]) < 1) {
				    validRates=0
				    break
				}
			    }
			}
		    }
		}
	    }
   	    model[7:12] <- formatC(rmatrix_scaled,digits=4)
    	}
    
       #Gamma shape from an Exponential (1,1) => mean=1
	gamma_shape <- rtrunc(1,"exp", 1, 2, a=0.5,b=5)
	model[13] <- formatC(gamma_shape,digits=4)
	       	
    model[14] <- model_partition[modelid]
    model[15] <- modelid
   
    model

}
