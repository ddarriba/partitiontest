R_LOG_FILE="rscript.log"
IND_LOG_FILE="indelible.log"
LOG_DIR="log"

prefix=testbed
datatype=aa
NUM_SIMS=100
mintaxa=10
maxtaxa=50
mingenes=5
maxgenes=50
minlen=500
maxlen=1500

simsdir=sims.$prefix

INPUT_TREEFILE="${simsdir}/treefile.out"
INPUT_MODELFILE="${simsdir}/modelsfile.out"
INPUT_GEN2PARTFILE="${simsdir}/genestopartitions.out"
INPUT_PARTITIONSFILE="${simsdir}/partitionsfile.out"

OUTPUT_DIR="indelible"
OUTPUT_ALIGNS_DIR="alignments"
OUTPUT_PARTEST_DIR="partitiontest"
OUTPUT_TRUE_DIR="truepart"
OUTPUT_PART_SUMMARY="scripts/sims/partitionssummary.out"

# Generate models, partitions and trees
echo "[1] Running R scripts"
cd R-scripts
R --vanilla --slave --args $prefix $datatype $NUM_SIMS $mintaxa $maxtaxa $mingenes $maxgenes $minlen $maxlen < build_models.r 2>&1
cd -


	echo "[2] Building INDELible control files"

	rm -rf ${LOG_DIR}
	mkdir ${LOG_DIR}
	rm -rf ${OUTPUT_DIR}
	mkdir ${OUTPUT_DIR}
	cp bin/indelible ${OUTPUT_DIR}

	AUX_FILE=${OUTPUT_DIR}/aux.txt
	part_header_line=1
	# Loop over simulations
	for ((sim_index=1; sim_index<=${NUM_SIMS}; sim_index++)); do

	  # Indelible control file
	  IND_FILE=${OUTPUT_DIR}/control${sim_index}.txt
	  touch ${IND_FILE}

	  # Write header
	  echo [TYPE] NUCLEOTIDE 2 >> ${IND_FILE}
	  echo " " >> ${IND_FILE}

	  # Write tree
	  echo [TREE] tree `sed "${sim_index}q;d" ${INPUT_TREEFILE}` >> ${AUX_FILE}

	  part_header=`sed "${part_header_line}q;d" ${INPUT_GEN2PARTFILE}`
	  num_genes=`echo ${part_header} | cut -d' ' -f2`
	  num_partitions=`echo ${part_header} | cut -d' ' -f3`
	  echo " " >> ${AUX_FILE}

	  echo "    ${sim_index}/${NUM_SIMS} : ${num_partitions} partitions out of ${num_genes} genes - $((num_partitions*100/num_genes))"

	  # Write partitions
	  echo [PARTITIONS] partitions >> ${AUX_FILE}
	  for ((partition_index=1; partition_index<=${num_genes}; partition_index++)); do
	    partmodel=`sed "$((part_header_line + partition_index))q;d" ${INPUT_GEN2PARTFILE}`
	    model_index=`echo ${partmodel} | cut -d' ' -f2`
	    gene_length=`echo ${partmodel} | cut -d' ' -f4`
	    echo [tree model${model_index} ${gene_length}] >> ${AUX_FILE}
	  done

	  #Write models
	  models=`sed -n -e $((part_header_line+1)),$((part_header_line+num_genes))p ${INPUT_GEN2PARTFILE} | cut -d' ' -f2`
	  models=`for each in ${models}; do echo $each; done | sort | uniq`
	  for each in ${models}; do
	    model=`sed -n -e ${each}p ${INPUT_MODELFILE}`
	    modelname=`echo $model | cut -d' ' -f1`
	    fA=`echo $model | cut -d' ' -f3`
	    fC=`echo $model | cut -d' ' -f4`
	    fG=`echo $model | cut -d' ' -f5`
	    fT=`echo $model | cut -d' ' -f6`
	    
	    kappa=`echo $model | cut -d' ' -f7`
	    
	    rA=`echo $model | cut -d' ' -f8`
	    rB=`echo $model | cut -d' ' -f9`
	    rC=`echo $model | cut -d' ' -f10`
	    rD=`echo $model | cut -d' ' -f11`
	    rE=`echo $model | cut -d' ' -f12`
	    rF=`echo $model | cut -d' ' -f13`

	    pInv=`echo $model | cut -d' ' -f14`
	    shape=`echo $model | cut -d' ' -f15`
	   
	    if [ $shape != 0 ]; then
	      ncat=4
	    else
	      ncat=0
	    fi
	    
	    echo [MODEL] model${each} >> ${IND_FILE}

	    if [ $rA != NA -a $rA != 0 ]; then
	        echo [submodel] GTR ${rE} ${rC} ${rF} ${rA} ${rD} ${rB} >> ${IND_FILE}
	    else
	        echo [submodel] HKY ${kappa} >> ${IND_FILE}
	    fi

	    if [ $shape == 100 ]; then
		shape=0
	    fi

	    echo [rates] ${pInv} ${shape} ${ncat} >> ${IND_FILE}
	    echo [statefreq] $fT $fC $fA $fG >> ${IND_FILE}

	    echo " " >> ${IND_FILE}
	  done

	  cat ${AUX_FILE} >> ${IND_FILE}

	  echo " " >> ${IND_FILE}
	  echo [EVOLVE] partitions 1 alignment${sim_index} >> ${IND_FILE}

	  rm ${AUX_FILE}

	  part_header_line=$((part_header_line + num_genes + 1))
	done
