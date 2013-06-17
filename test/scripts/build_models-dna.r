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

SAMPLES <- 20  # Total number of samples
TAXA_COUNT <- 20  # Number of taxa in the output trees
GENE_LEN <- 500   # Length of each gene

MIN_GENES <- 1000
MAX_GENES <- 1000
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
source("models.r")

# delete output files
unlink(OUT_TREE_FILE)
unlink(OUT_MODELS_FILE)
unlink(OUT_PARTITIONS_FILE)
unlink(OUT_GENTOPART_FILE)

# load models and properties from data.in
inp <- scan(IN_MODELS_FILE,list("",0,0,0,0,0,0,0,0,"",0))

current_index <- 0


for(sample_index in 0:(SAMPLES-1)) {

	num_genes <- 1000 #sample(MIN_GENES:MAX_GENES,1,replace=T)

	# Assign models to genes
	boxes <- rchinese(num_genes,sample(1:20,1))
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

	num_partitions <- length(unique(boxes))

	# allocate results matrix
	models_mat <- matrix(nrow=(num_partitions),ncol=17,byrow=TRUE)
	colnames(models_mat) <- c("ModelName", "Params", "fA","fC","fG","fT","kappa", "Ra", "Rb", "Rc", "Rd", "Re", "Rf", "pInv", "shape", "partition", "id")
	models_mat <- as.table(models_mat)

	# construct one model per partition
	avoid <- c()
	for(i in 1:num_partitions) {

		model <- buildDNAmodel(inp, avoid)
		models_mat[i,1:17] <- model
		id <- as.numeric(model[17])

		if (id <= 80) {
			avoid <- c(id, avoid)
			if (id %% 2 > 0) {
				avoid <- c( id + 1, avoid)
			} else {
				avoid <- c( id - 1, avoid)
			}
		}
	}

	# Write the output models
	write.table(models_mat,file=OUT_MODELS_FILE,append=TRUE,quote=FALSE,col.names=FALSE,row.names=FALSE)

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
	header[3] <- num_partitions
	write(header, file=OUT_GENTOPART_FILE, append=TRUE)
	write.table(genes_mat,file=OUT_GENTOPART_FILE,append=TRUE,quote=FALSE,col.names=FALSE,row.names=FALSE)

	write(header, file=OUT_PARTITIONS_FILE, append=TRUE)
	write(partitionstringDEC, file=OUT_PARTITIONS_FILE, append=TRUE)
	write(partitionstring, file=OUT_PARTITIONS_FILE, append=TRUE)

	current_index <- current_index + num_partitions


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






