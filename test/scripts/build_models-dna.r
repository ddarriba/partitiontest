# Sampling model parameters for models of substitution
# May 2013
# D. Posada
# D. Darriba

# This script generates 4 files:
#   modelsfile.out         list of amino-acid replacement models
#   genestopartitions.out  genes mapping to partitions (thus, models)
#   partitionsfile.out     definition of partitions
#   treeefile.out          one tree per partitioning scheme (i.e., simulation)

library (MCMCpack)
library (ape)
library (phylosim)

SAMPLES <- 1000   # Total number of samples
TAXA_COUNT <- 20  # Number of taxa in the output trees
GENE_LEN <- 500   # Length of each gene

MIN_GENES <- 10
MAX_GENES <- 100
MIN_PARTITIONS <- 1

IN_MODELS_FILE  <- "data.11.in"      # Data input file

OUT_TREE_FILE   <- "sims/treefile.out"    # Output file for trees
OUT_MODELS_FILE <- "sims/modelsfile.out"  # Output file for models
OUT_GENTOPART_FILE <- "sims/genestopartitions.out"  # Output file for gene mapping to partitions
OUT_PARTITIONS_FILE <- "sims/partitionsfile.out"  # Output file for partitions
# GENTOPART file has the following format
# N(NUM_GENES) K(NUM_PARTITIONS)
# 0	MODEL_TO_GENE0
# 1	MODEL_TO_GENE1
# ...
# N	MODEL_TO_GENEN
# N(NUM_GENES) K(NUM_PARTITIONS)
# 0	MODEL_TO_GENE0
# 1	MODEL_TO_GENE1
# ...
# N	MODEL_TO_GENEN

source("truncated.r")
source("clustering.r")

# delete output files
unlink(OUT_TREE_FILE)
unlink(OUT_MODELS_FILE)
unlink(OUT_PARTITIONS_FILE)
unlink(OUT_GENTOPART_FILE)

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

indelible_model <- inp[[1]]
indelible_index <- inp[[2]]
phyml_model <- inp[[3]]


current_index <- 0


for(sample_index in 0:(SAMPLES-1)) {

	num_genes <- sample(MIN_GENES:MAX_GENES,1,replace=T)

	# allocate results matrix
	models_mat <- matrix(nrow=(num_genes),ncol=16,byrow=TRUE)
	colnames(models_mat) <- c("ModelName", "Params", "fA","fC","fG","fT","kappa", "Ra", "Rb", "Rc", "Rd", "Re", "Rf", "pInv", "shape", "partition")
	models_mat <- as.table(models_mat)

	# construct one model per partition
	for(i in 1:(num_genes)) {

	    modelid <- sample(1:88,1,replace=T)

	    f_model <- equal_base_frequencies[modelid]
	    g_model <- (rate_variation[modelid]==2 || rate_variation[modelid]==3)
	    i_model <- (rate_variation[modelid]==1 || rate_variation[modelid]==3)


	    models_mat[i,1] <- model[modelid]
	    models_mat[i,2] <- free_params[modelid]
	    
	    completename <- phyml_model[modelid]
	    if (f_model == 1) {
	    	models_mat[i,3:6] <- 0.25
	    } else {
	    	#Base frequencies from a Dirichlet(1.0,1.0,1.0,1.0)
			base_frequencies <- rdirichlet(1, c(1,1,1,1) )
			models_mat[i,3:6] <- formatC(base_frequencies,digits=4)
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
	    models_mat[i,7] <- formatC(kappa,digits=4)
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
   	    models_mat[i,8:13] <- formatC(rmatrix_scaled,digits=4)
    	}
    
    	if (rate_variation[modelid]==1 || rate_variation[modelid]==3) {
	    #Proportion of invariable sites from uniform (0,1)
	    pinv <- rtrunc(1,"beta", 1, 3, a=0.2,b=0.8)
	    models_mat[i,14] <- formatC(pinv,digits=4)
    	} else {
    	    models_mat[i,14] <- 0
    	}
    	
	if (rate_variation[modelid]==2 || rate_variation[modelid]==3) {
	       #Gamma shape from an Exponential (1,1) => mean=1
			gamma_shape <- rtrunc(1,"exp", 1, 2, a=0.5,b=5)
			models_mat[i,15] <- formatC(gamma_shape,digits=4)
    	} else {
    		models_mat[i,15] <- 100
    	}
	       	
    models_mat[i,16] <- model_partition[modelid]

    }

    # Write the output models
    write.table(models_mat,file=OUT_MODELS_FILE,append=TRUE,quote=FALSE,col.names=FALSE,row.names=FALSE)

	# Assign models to genes
	boxes <- rchinese(num_genes,sample(5:100,1))
	# boxes <- cluster(max_partitions,num_genes)
	genes_mat <- matrix(nrow=num_genes,ncol=2,byrow=TRUE)
	colnames(genes_mat) <- c("GeneNumber", "ModelNumber")
	genes_mat <- as.table(genes_mat)
	hash <- matrix(0, 1, length(unique(boxes)))
	hashsize <- 1
	hash[1] <- boxes[1]

	partitionstring <- matrix(1:num_genes)
	for(i in 1:(num_genes)) {
	  genes_mat[i,1] <- i-1
	  genes_mat[i,2] <- boxes[i] + current_index
	  indexed <- F
	  cur_index <- -1
	  for(j in 1:hashsize) {
	    if (hash[j] == boxes[i]) {
	      cur_index <- j
	      indexed <- T
	    }
	  }
	  if (indexed == F) {
	    hashsize <- hashsize + 1
	    cur_index <- hashsize
	    hash[cur_index] <- boxes[i]
	  }
	  partitionstring[i] <- cur_index-1
	}
	partitionstringDEC <- matrix(1:num_genes)
	for(i in 1:(num_genes)) {
	  partitionstringDEC[i] <- partitionstring[i] %/% 10
	  partitionstring[i] <- partitionstring[i] %% 10
	}
	partitionstring <- toString(partitionstring)
	partitionstring <- gsub(", ","",partitionstring)
	partitionstringDEC <- toString(partitionstringDEC)
	partitionstringDEC <- gsub(", ","",partitionstringDEC)

	header <- array(1:3)
	header[1] <- sample_index
	header[2] <- num_genes
	header[3] <- length(unique(boxes))
	write(header, file=OUT_GENTOPART_FILE, append=TRUE)
	write.table(genes_mat,file=OUT_GENTOPART_FILE,append=TRUE,quote=FALSE,col.names=FALSE,row.names=FALSE)

	write(header, file=OUT_PARTITIONS_FILE, append=TRUE)
	write(partitionstringDEC, file=OUT_PARTITIONS_FILE, append=TRUE)
	write(partitionstring, file=OUT_PARTITIONS_FILE, append=TRUE)

	current_index <- current_index + num_genes 


	   ###SIMULATE TREE
    # Generate a random non-ultrametric tree with branches according to a exponential with mean 1/10=0.1. Note that a rooted tree has 2n-2 branches (40 taxa => 78 branches => expected treeLength = 7.8) although this does not matter because we are going to scale the treelength next) 
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

} # end SAMPLES






