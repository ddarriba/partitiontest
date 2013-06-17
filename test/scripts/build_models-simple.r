# Sampling model parameters for models of substitution
# January 2012
# D. Posada
# D. Darriba

library (MCMCpack)
library (ape)
library (phylosim)

#### CONSTANTS #### (Actually they are not, but I encourage you to
                  # not change their value outside this declaration)
SAMPLES <- 10  # Total number of samples of each model group. i.e, this 
                 # script will generate a result set with 4 times SAMPLES
                 # elements
IN_MODELS_FILE  <- "data.in"         # Data input file
OUT_TREE_FILE   <- "treefile.out"    # Output file for models
OUT_MODELS_FILE <- "modelsfile.out"  # Output file for trees
TAXA_COUNT <- 20  # Number of taxa in the output trees

source("truncated.r")

# delete output files
unlink(OUT_TREE_FILE)
unlink(OUT_MODELS_FILE)

# load models and properties from data.in
inp <- scan(IN_MODELS_FILE,list("",0,0,0,0,0,0,0,0,"",0))
model <- inp[[1]]
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

# allocate results matrix
results <- matrix(nrow=(SAMPLES*4),ncol=16,byrow=TRUE)
colnames(results) <- c("ModelName", "Params", "fA","fC","fG","fT","titv", "Ra", "Rb", "Rc", "Rd", "Re", "Rf", "pInv", "shape", "partition")
results <- as.table(results)

# main loop (SAMPLES times)
for(i in 0:(SAMPLES-1)) {
# nested loop (iterates for each group of models: M, M+I, M+G, M+I+G)
for(j in 1:4) {
    modelid <- sample(0:21,1,replace=T)*4+j
    results[i*4+j,1] <- model[modelid]
    results[i*4+j,2] <- free_params[modelid]

    if (equal_base_frequencies[modelid] == 1) {
	results[i*4+j,3:6] <- 0.25
    } else {
	#Base frequencies from a Dirichlet(1.0,1.0,1.0,1.0)
	base_frequencies <- rdirichlet(1, c(1,1,1,1) )
	results[i*4+j,3:6] <- formatC(base_frequencies,digits=4)
    }

    if (model_titv[modelid] == 1) {
	if (modelid > 8) {
	  #titv <- rtrunc(1,"norm", 2, 2, a=0.5, b=10)
	  #Transition/Transversion rate from a Gamma (2,1) truncated between 2 and 10
	  titv <- rtrunc(1,"gamma", 2, 1, a=2, b=10)
	} else {
	  titv <- 0.5
	}
	results[i*4+j,7] <- formatC(titv,digits=4)
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
   	results[i*4+j,8:13] <- formatC(rmatrix_scaled,digits=4)
    }

    if (rate_variation[modelid]==1 || rate_variation[modelid]==3) {
	#Proportion of invariable sites from uniform (0,1)
	pinv <- rtrunc(1,"beta", 1, 3, a=0.2,b=0.8)
	results[i*4+j,14] <- formatC(pinv,digits=4)
    }

    if (rate_variation[modelid]==2 || rate_variation[modelid]==3) {
	#Gamma shape from an Exponential (1,1) => mean=1
	gamma_shape <- rtrunc(1,"exp", 1, 2, a=0.5,b=5)
	results[i*4+j,15] <- formatC(gamma_shape,digits=4)
    }

    results[i*4+j,16] <- model_partition[modelid]

    ###SIMULATE TREE
    # Generate a random non-ultrametric tree with branches according to a exponential with mean 1/10=0.1. Note that a rooted tree has 2n-2 branches (40 taxa => 78 branches => expected treeLength = 7.8) although this does not matter because we are going to scale the treelength next) 
    TAXA_COUNT <- sample(10:100,1)
    tree <- rtree(TAXA_COUNT,TRUE,NULL,rexp,rate=10)
    PhyloSim(tree)$treeLength
    write.tree(tree)

    # Scale total tree length so the tree length uniformly distributed in the [0.5, 10] range
    runiform <- runif(1,0,1)
    scale_factor <- (2.0 + 10.0*runiform)
    scaledPhyloTree<-PhyloSim(tree)
    scaleTree(scaledPhyloTree,scale_factor/scaledPhyloTree$treeLength)
    scaledPhyloTree$treeLength
    scaledTree <- getPhylo(scaledPhyloTree)
    write.tree(scaledTree, file=OUT_TREE_FILE, append=TRUE)
}
}

# Write the output models
write.table(results,file=OUT_MODELS_FILE,append=FALSE,quote=FALSE,col.names=FALSE,row.names=FALSE)









