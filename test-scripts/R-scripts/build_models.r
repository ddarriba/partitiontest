# Sampling model parameters for models of substitution
# May 2013
# D. Posada
# D. Darriba

# This script generates 4 files:
#   modelsfile.out          list of amino-acid replacement models
#   blockstopartitions.out  blocks mapping to partitions (thus, models)
#   partitionsfile.out      definition of partitions
#   treeefile.out           one tree per partitioning scheme (i.e., simulation)

library (MCMCpack)
library (ape)
library (phylosim)

options("scipen"=100, "digits"=6)

args = commandArgs(trailingOnly = T)
if (length(args) < 9) {
  cat("Number of arguments is ", length(args), "\n")
  stop("Required arguments: Prefix, DataType={aa, dna}, Num_Samples, Min_Taxa, Max_Taxa, Min_DataBlocks, Max_DataBlocks, Min_BlockLen, Max_BlockLen")
}
PREFIX      = args[1]
DATA_TYPE   = args[2]
SAMPLES     = as.numeric(args[3])   # Total number of samples
MIN_TAXA    = as.numeric(args[4])   # Number of taxa in the output trees
MAX_TAXA    = as.numeric(args[5])   
MIN_DATA_BLOCKS   = as.numeric(args[6])	# Number of data blocks
MAX_DATA_BLOCKS   = as.numeric(args[7])
MIN_BLOCKLEN = as.numeric(args[8])	# Length of each data block
MAX_BLOCKLEN = as.numeric(args[9])

cat("\n------ Sim parameters ------\n")
cat("Exec Prefix:       ",PREFIX,"\n")
cat("Number of samples: ",SAMPLES,"\n")
cat("Number of taxa:    ",MIN_TAXA,"-",MAX_TAXA,"\n")
cat("Number of blocks:  ",MIN_DATA_BLOCKS,"-",MAX_DATA_BLOCKS,"\n")
cat("Block length:      ",MIN_BLOCKLEN,"-",MAX_BLOCKLEN,"\n")
cat("------ -------------- ------\n\n")

MIN_PARTITIONS = 1
MAX_PARTITIONS = MAX_DATA_BLOCKS

BASE_DIR = paste("../sims.",PREFIX, sep="")
OUT_TREE_FILE   = paste(BASE_DIR,"/treefile.out", sep="")               # Output file for trees
OUT_SUMMARY_FILE   = paste(BASE_DIR,"/true.summary", sep="")
OUT_MODELS_FILE = paste(BASE_DIR,"/modelsfile.out", sep="")             # Output file for models
OUT_GENTOPART_FILE = paste(BASE_DIR,"/blockstopartitions.out", sep="")  # Output file for data block mapping to partitions
OUT_PARTITIONS_FILE = paste(BASE_DIR,"/partitionsfile.out", sep="")     # Output file for partitions
OUT_SIMS_DIR = paste(BASE_DIR,"/simulations/", sep="")

source("truncated.r")
source("clustering.r")
source("models.r")

# delete output files
unlink(BASE_DIR, recursive=TRUE)
dir.create(BASE_DIR)
dir.create(OUT_SIMS_DIR)

# load models and properties from data.in
if (DATA_TYPE == "aa") {
	inp = scan("data-aa.in",list(indelible_model="",indelible_index=0,""), quiet=TRUE)
} else {
	inp = scan("data-dna.in",list("",0,0,0,0,0,0,0,0,"",0), quiet=TRUE)
}

current_index = 0

for(sample_index in 1:(SAMPLES)) {

	if (MIN_DATA_BLOCKS < MAX_DATA_BLOCKS) {
		num_data_blocks = sample(MIN_DATA_BLOCKS:MAX_DATA_BLOCKS,1,replace=T)
	} else {
		num_data_blocks = MIN_DATA_BLOCKS
	}
	
	# Assign models to data blocks
  n_partitions = sample(MIN_PARTITIONS:MAX_PARTITIONS,1)
	boxes = cluster(n_partitions,num_data_blocks)
	data_blocks_mat = matrix(nrow=num_data_blocks,ncol=4,byrow=TRUE)
	colnames(data_blocks_mat) = c("BlockNumber", "ModelNumber", "LocalNumber", "NumSites")
	data_blocks_mat = as.table(data_blocks_mat)
	hash = matrix(0, 1, length(unique(boxes)))
	hashsize = 1
	hash[1] = boxes[1]

	converttable = numeric(length(boxes))
	cur_index=1
	for (i in 1:length(boxes)) {
  	  if (converttable[boxes[i]] == 0) {
	    converttable[boxes[i]] = cur_index
	    cur_index = cur_index+1
	  }
	  boxes[i] = converttable[boxes[i]]
	}

	OUT_CURSIM_FILE = paste(OUT_SIMS_DIR, "/alignment", sample_index, ".desc", sep="")
	sim_desc = matrix(nrow=(num_data_blocks),ncol=6,byrow=TRUE)
	colnames(sim_desc) = c("block", "partid", "model", "length", "start", "end")

	num_partitions = length(unique(boxes))

	# allocate results matrix
	if (DATA_TYPE == "aa") {
		models_mat = matrix(nrow=(num_partitions),ncol=26,byrow=TRUE)
		colnames(models_mat) = c("ModelName", "+F",
                             "f1", "f2", "f3", "f4", "f5",
                             "f6", "f7", "f8", "f9", "f10",
                             "f11","f12","f13","f14","f15",
                             "f16","f17","f18","f19","f20",
                             "pinv", "shape", "partition", "id")
	} else {
		models_mat = matrix(nrow=(num_partitions),ncol=16,byrow=TRUE)
		colnames(models_mat) = c("ModelName", "fA","fC","fG","fT","kappa",
                             "Ra", "Rb", "Rc", "Rd", "Re", "Rf", 
                             "pinv", "shape", "partition", "id")
	}
	models_mat = as.table(models_mat)

	# construct one model per partition

	avoid = c()
	for(i in 1:(num_partitions)) {
	  if (DATA_TYPE == "aa") {
		  f_model = rbinom(1, 1, 0.5)
		  g_model = rbinom(1, 1, 0.5)
		  i_model = rbinom(1, 1, 0.5)

		  modelid = sample(1:14,1,replace=T)

		  models_mat[i,25] = inp$indelible_model[modelid]
		  models_mat[i,26] = inp$indelible_index[modelid]
		  completename = inp$indelible_model[modelid]
		  completename = paste(completename, "+G", sep="")
		  if (f_model == 1) {
		   completename = paste(completename, "+F", sep="")
		   base_frequencies <- numeric(20)
		   while (min(base_frequencies) < 0.005) {
		       base_frequencies <- rdirichlet(1, c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1) )
		   }
		   models_mat[i,3:22] <- formatC(base_frequencies,digits=6)
		  }
		  models_mat[i,1] = completename
		  models_mat[i,2] = f_model

      if (i_model)
  		  p_inv = rtrunc(1,"exp", 1, 2, a=0.2,b=0.8)
      else
        p_inv = 0
		  models_mat[i,23] = formatC(p_inv,digits=4)

      if (g_model)
  		  gamma_shape = rtrunc(1,"exp", 1, 2, a=0.5,b=5)
      else
        gamma_shape = 100
		  models_mat[i,24] = formatC(gamma_shape,digits=4)

    } else {

		  model = buildDNAmodel(inp, avoid)
		  models_mat[i,1:16] = model
		  id = as.numeric(model[16])
		  if ((id %% 22) < 13) {
			  avoid = c(id, avoid)
#			  if (id %% 2 > 0) {
#				  avoid = c( id + 1, avoid)
#			  } else {
#				  avoid = c( id - 1, avoid)
#			  }

		  }
	  }
	}

	# Write the output models
	write.table(models_mat,file=OUT_MODELS_FILE,append=TRUE,quote=FALSE,col.names=FALSE,row.names=FALSE)

	partitionstring = matrix(1:num_data_blocks)
	startposition = 1
	for(i in 1:(num_data_blocks)) {

	  data_blocks_mat[i,1] = i-1
	  data_blocks_mat[i,2] = boxes[i] + current_index
	  data_blocks_mat[i,3] = boxes[i]
	  if (MIN_BLOCKLEN < MAX_BLOCKLEN) {
	  	data_blocks_mat[i,4] = sample(MIN_BLOCKLEN:MAX_BLOCKLEN,1,replace=T)
	  } else {
		data_blocks_mat[i,4] = MIN_BLOCKLEN
	  }

	  sim_desc[i,1] = i
	  sim_desc[i,2] = boxes[i] + current_index
	  sim_desc[i,3] = models_mat[boxes[i],1]
	  sim_desc[i,4] = data_blocks_mat[i,4]
	  sim_desc[i,5] = startposition
	  sim_desc[i,6] = startposition + data_blocks_mat[i,4] - 1
	  startposition = startposition + data_blocks_mat[i,4]

	  partitionstring[i] = boxes[i]-1
	}
	write.table(sim_desc,file=OUT_CURSIM_FILE,append=FALSE,quote=FALSE,col.names=FALSE,row.names=FALSE)

	partitionstring100 = matrix(1:num_data_blocks)
	partitionstring10  = matrix(1:num_data_blocks)
	for(i in 1:(num_data_blocks)) {
	  partitionstring100[i] = partitionstring[i] %/% 100
	  partitionstring10[i] = (partitionstring[i] %% 100) %/% 10
	  partitionstring[i] = partitionstring[i] %% 10
	}
	partitionstring = toString(partitionstring)
	partitionstring = gsub(", ","",partitionstring)
	partitionstring10 = toString(partitionstring10)
	partitionstring10 = gsub(", ","",partitionstring10)
	partitionstring100 = toString(partitionstring100)
	partitionstring100 = gsub(", ","",partitionstring100)

	header = array(1:3)
	header[1] = sample_index
	header[2] = num_data_blocks
	header[3] = num_partitions
	write(header, file=OUT_GENTOPART_FILE, append=TRUE)
	write.table(data_blocks_mat,file=OUT_GENTOPART_FILE,append=TRUE,quote=FALSE,col.names=FALSE,row.names=FALSE)

	write(header, file=OUT_PARTITIONS_FILE, append=TRUE)
	write(partitionstring100, file=OUT_PARTITIONS_FILE, append=TRUE)
	write(partitionstring10, file=OUT_PARTITIONS_FILE, append=TRUE)
	write(partitionstring, file=OUT_PARTITIONS_FILE, append=TRUE)

	current_index = current_index + num_partitions

	   ###SIMULATE TREE
    # Generate a random non-ultrametric tree with branches according to a exponential with mean 1/10=0.1. Note that a rooted tree has 2n-2 branches (40 taxa => 78 branches => expected treeLength = 7.8) although this does not matter because we are going to scale the treelength next) 
    if (MIN_TAXA < MAX_TAXA) {
	TAXA_COUNT = sample(MIN_TAXA:MAX_TAXA,1,replace=T)
    } else {
	TAXA_COUNT = MIN_TAXA
    }
    goodTree = FALSE
    tree_it = 400
    while (!goodTree) {
      tree = rtree(TAXA_COUNT,TRUE,NULL,rexp,rate=10)
      PhyloSim(tree)$treeLength
    
      write.tree(tree)

      # Scale total tree length so the tree length uniformly distributed in the [2, 12] range
      min_treelen = 2
      max_treelen = 12
      runiform = runif(1,min_treelen,max_treelen)
      nbranches = 2*TAXA_COUNT - 3
      scale_factor = runiform # ((.1 + runiform)*nbranches)
      scaledPhyloTree=PhyloSim(tree)
      scaleTree(scaledPhyloTree,scale_factor/scaledPhyloTree$treeLength)
      scaledPhyloTree$treeLength
      scaledTree = getPhylo(scaledPhyloTree)
      goodTree = (max(tree$edge.length)/min(tree$edge.length) < 1000 && max(scaledTree$edge.length) < 2 && min(scaledTree$edge.length)>0.001)

      # for security (avoid infinite loops)
      tree_it = tree_it - 1
      stopifnot(tree_it > 0)
    }
    write.tree(scaledTree, file=OUT_TREE_FILE, append=TRUE)

    treeStr = write.tree(scaledTree)
    treeLen = sum(scaledTree$edge.length)

    tList = list(id=sample_index, ndatablocks=num_data_blocks, nparts=num_partitions,  
                 part2=partitionstring100, part1=partitionstring10, part0=partitionstring, 
                 ntaxa=TAXA_COUNT, seqlen=(startposition-1),treelen=treeLen,tree=treeStr)
    write.table(tList, file=OUT_SUMMARY_FILE, append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE)
    cat(sample_index,"/",SAMPLES, " : ntaxa=",TAXA_COUNT, " ndatablocks=", num_data_blocks, " nparts=", num_partitions,"\n")
} # end SAMPLES
cat("\nDone R script\n")






