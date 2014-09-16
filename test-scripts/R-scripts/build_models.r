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

options("scipen"=100, "digits"=6)

args = commandArgs(trailingOnly = T)
if (length(args) < 9) {
  cat("Number of arguments is ", length(args), "\n")
  stop("Required arguments: Prefix, DataType={aa, dna}, Num_Samples, Min_Taxa, Max_Taxa, Min_Genes, Max_Genes, Min_GeneLen, Max_GeneLen")
}
PREFIX      = args[1]
DATA_TYPE   = args[2]
SAMPLES     = as.numeric(args[3])   # Total number of samples
MIN_TAXA    = as.numeric(args[4])   # Number of taxa in the output trees
MAX_TAXA    = as.numeric(args[5])   # Length of each gene
MIN_GENES   = as.numeric(args[6])
MAX_GENES   = as.numeric(args[7])
MIN_GENELEN = as.numeric(args[8])
MAX_GENELEN = as.numeric(args[9])

cat("\n------ Sim parameters ------\n")
cat("Exec Prefix:       ",PREFIX,"\n")
cat("Number of samples: ",SAMPLES,"\n")
cat("Number of taxa:    ",MIN_TAXA,"-",MAX_TAXA,"\n")
cat("Number of genes:   ",MIN_GENES,"-",MAX_GENES,"\n")
cat("Gene length:       ",MIN_GENELEN,"-",MAX_GENELEN,"\n")
cat("------ -------------- ------\n\n")

MIN_PARTITIONS = 1
MAX_PARTITIONS = 100

BASE_DIR = paste("../sims.",PREFIX, sep="")
OUT_TREE_FILE   = paste(BASE_DIR,"/treefile.out", sep="")              # Output file for trees
OUT_SUMMARY_FILE   = paste(BASE_DIR,"/true.summary", sep="")
OUT_MODELS_FILE = paste(BASE_DIR,"/modelsfile.out", sep="")            # Output file for models
OUT_GENTOPART_FILE = paste(BASE_DIR,"/genestopartitions.out", sep="")  # Output file for gene mapping to partitions
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

	if (MIN_GENES < MAX_GENES) {
		num_genes = sample(MIN_GENES:MAX_GENES,1,replace=T)
	} else {
		num_genes = MIN_GENES
	}
	max_partitions = sample(MIN_PARTITIONS:min(num_genes,MAX_PARTITIONS),1,replace=T)

	

	# Assign models to genes
	boxes = cluster(max_partitions,num_genes)
	genes_mat = matrix(nrow=num_genes,ncol=4,byrow=TRUE)
	colnames(genes_mat) = c("GeneNumber", "ModelNumber", "LocalNumber", "NumSites")
	genes_mat = as.table(genes_mat)
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
	sim_desc = matrix(nrow=(num_genes),ncol=6,byrow=TRUE)
	colnames(sim_desc) = c("gene", "partid", "model", "length", "start", "end")

	num_partitions = length(unique(boxes))

	# allocate results matrix
	models_mat = matrix(nrow=(num_partitions),ncol=15,byrow=TRUE)
	if (DATA_TYPE == "aa") {
		colnames(models_mat) = c("ModelName", "+F","","","", "", "","","","","","", "shape", "partition", "id")
	} else {
		colnames(models_mat) = c("ModelName", "fA","fC","fG","fT","kappa", "Ra", "Rb", "Rc", "Rd", "Re", "Rf", "shape", "partition", "id")
	}
	models_mat = as.table(models_mat)

	# construct one model per partition

	avoid = c()
	for(i in 1:(num_partitions)) {
	    if (DATA_TYPE == "aa") {
		f_model = rbinom(1, 1, 0.5)

		modelid = sample(1:14,1,replace=T)

		models_mat[i,14] = inp$indelible_model[modelid]
		models_mat[i,15] = inp$indelible_index[modelid]
		completename = inp$indelible_model[modelid]
		completename = paste(completename, "+G", sep="")
		if (f_model == 1) {
		   completename = paste(completename, "+F", sep="")
		}
		models_mat[i,1] = completename
		models_mat[i,2] = f_model

		#Gamma shape from an Exponential (1,1) => mean=1
		gamma_shape = rtrunc(1,"exp", 1, 2, a=0.5,b=5)
		models_mat[i,13] = formatC(gamma_shape,digits=4)
	    } else {
		model = buildDNAmodel(inp, avoid)
		models_mat[i,1:15] = model
		id = as.numeric(model[15])
		if (id <= 20) {
			avoid = c(id, avoid)
			if (id %% 2 > 0) {
				avoid = c( id + 1, avoid)
			} else {
				avoid = c( id - 1, avoid)
			}
		}
	    }
	}

	# Write the output models
	write.table(models_mat,file=OUT_MODELS_FILE,append=TRUE,quote=FALSE,col.names=FALSE,row.names=FALSE)

	partitionstring = matrix(1:num_genes)
	startposition = 1
	for(i in 1:(num_genes)) {

	  genes_mat[i,1] = i-1
	  genes_mat[i,2] = boxes[i] + current_index
	  genes_mat[i,3] = boxes[i]
	  if (MIN_GENELEN < MAX_GENELEN) {
	  	genes_mat[i,4] = sample(MIN_GENELEN:MAX_GENELEN,1,replace=T)
	  } else {
		genes_mat[i,4] = MIN_GENELEN
	  }

	  sim_desc[i,1] = i
	  sim_desc[i,2] = boxes[i] + current_index
	  sim_desc[i,3] = models_mat[boxes[i],1]
	  sim_desc[i,4] = genes_mat[i,4]
	  sim_desc[i,5] = startposition
	  sim_desc[i,6] = startposition + genes_mat[i,4] - 1
	  startposition = startposition + genes_mat[i,4]

	  partitionstring[i] = boxes[i]-1
	}
	write.table(sim_desc,file=OUT_CURSIM_FILE,append=FALSE,quote=FALSE,col.names=FALSE,row.names=FALSE)

	partitionstringDEC = matrix(1:num_genes)
	for(i in 1:(num_genes)) {
	  partitionstringDEC[i] = partitionstring[i] %/% 10
	  partitionstring[i] = partitionstring[i] %% 10
	}
	partitionstring = toString(partitionstring)
	partitionstring = gsub(", ","",partitionstring)
	partitionstringDEC = toString(partitionstringDEC)
	partitionstringDEC = gsub(", ","",partitionstringDEC)

	header = array(1:3)
	header[1] = sample_index
	header[2] = num_genes
	header[3] = num_partitions
	write(header, file=OUT_GENTOPART_FILE, append=TRUE)
	write.table(genes_mat,file=OUT_GENTOPART_FILE,append=TRUE,quote=FALSE,col.names=FALSE,row.names=FALSE)

	write(header, file=OUT_PARTITIONS_FILE, append=TRUE)
	write(partitionstringDEC, file=OUT_PARTITIONS_FILE, append=TRUE)
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
    while (!goodTree) {
      tree = rtree(TAXA_COUNT,TRUE,NULL,rexp,rate=10)
      PhyloSim(tree)$treeLength
    
      write.tree(tree)

      # Scale total tree length so the tree length uniformly distributed in the [0.5, 10] range
      runiform = runif(1,0,1)
      nbranches = 2*TAXA_COUNT - 3
      scale_factor = ((.1 + runiform)*nbranches)
      scaledPhyloTree=PhyloSim(tree)
      scaleTree(scaledPhyloTree,scale_factor/scaledPhyloTree$treeLength)
      scaledPhyloTree$treeLength
      scaledTree = getPhylo(scaledPhyloTree)
      goodTree = (max(tree$edge.length)/min(tree$edge.length) < 1000 && max(scaledTree$edge.length) < 2 && min(scaledTree$edge.length)>0.001)
# if (!goodTree)
# cat("Reply ", TAXA_COUNT, " " , max(tree$edge.length), " ", min(tree$edge.length), " ", max(tree$edge.length)/min(tree$edge.length), " ", max(scaledTree$edge.length), " ", max(scaledTree$edge.length < 2) , " " , min(scaledTree$edge.length), " ", min(scaledTree$edge.length>0.001), "\n") 
    }
    write.tree(scaledTree, file=OUT_TREE_FILE, append=TRUE)

    treeStr = write.tree(scaledTree)
    treeLen = sum(scaledTree$edge.length)

    tList = list(id=sample_index, ngenes=num_genes, nparts=num_partitions, part0=partitionstring, 
      part1=partitionstringDEC, ntaxa=TAXA_COUNT, seqlen=(startposition-1),treelen=treeLen,tree=treeStr)
    write.table(tList, file=OUT_SUMMARY_FILE, append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE)
    cat(sample_index,"/",SAMPLES, " : ntaxa=",TAXA_COUNT, " ngenes=", num_genes, " nparts=", num_partitions,"\n")
} # end SAMPLES
cat("\nDone R script\n")






