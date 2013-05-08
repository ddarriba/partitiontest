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

SAMPLES <- 10    # Total number of samples
TAXA_COUNT <- 30  # Number of taxa in the output trees
GENE_LEN <- 1000  # Length of each gene

MIN_GENES <- 2
MAX_GENES <- 100
MIN_PARTITIONS <- 1

IN_MODELS_FILE  <- "data-aa.in"      # Data input file

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
inp <- scan(IN_MODELS_FILE,list("",0,""))
indelible_model <- inp[[1]]
indelible_index <- inp[[2]]
phyml_model <- inp[[3]]


current_index <- 0


for(sample_index in 0:(SAMPLES-1)) {

	num_genes <- sample(MIN_GENES:MAX_GENES,1,replace=T)
	max_partitions <- sample(MIN_PARTITIONS:num_genes,1,replace=T)

	# allocate results matrix
	models_mat <- matrix(nrow=(max_partitions),ncol=9,byrow=TRUE)
	colnames(models_mat) <- c("ModelName", "ModelIndex", "PhymlName","CompleteName", "+F", "+I", "+G", "pInv", "shape")
	models_mat <- as.table(models_mat)

	# construct one model per partition
	for(i in 1:(max_partitions)) {
	    f_model <- rbinom(1, 1, 0.5)
	    g_model <- rbinom(1, 1, 0.5)
	    i_model <- rbinom(1, 1, 0.5)

	    modelid <- sample(1:14,1,replace=T)

	    models_mat[i,1] <- indelible_model[modelid]
	    models_mat[i,2] <- indelible_index[modelid]
	    models_mat[i,3] <- phyml_model[modelid]
	    completename <- phyml_model[modelid]
	    if (i_model == 1) {
	       completename <- paste(completename, "+I", sep="")
	    }
	    if (g_model == 1) {
	       completename <- paste(completename, "+G", sep="")
	    }
	    if (f_model == 1) {
	       completename <- paste(completename, "+F", sep="")
	    }
	    models_mat[i,4] <- completename
	    models_mat[i,5] <- f_model
	    models_mat[i,6] <- i_model
	    models_mat[i,7] <- g_model

	    if (i_model==1) {
		#Proportion of invariable sites from uniform (0,1)
		pinv <- rtrunc(1,"beta", 1, 3, a=0.2,b=0.8)
		models_mat[i,8] <- formatC(pinv,digits=4)
	    } else {
		models_mat[i,8] <- 0.0
	    }
	
	    if (g_model == 1) {
		#Gamma shape from an Exponential (1,1) => mean=1
		gamma_shape <- rtrunc(1,"exp", 1, 2, a=0.5,b=5)
		models_mat[i,9] <- formatC(gamma_shape,digits=4)
	    } else {
		models_mat[i,9] <- 100.0
	    }

	}

	# Write the output models
	write.table(models_mat,file=OUT_MODELS_FILE,append=TRUE,quote=FALSE,col.names=FALSE,row.names=FALSE)

	# Assign models to genes
	boxes <- cluster(max_partitions,num_genes)
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

	current_index <- current_index + max_partitions


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






