PERGENE_LENGTH=1000

# JOB CONTROL
RUN_R_SCRIPT=0 # [1] Execute R scripts for defining the simulations
BUILD_INDEL_FILES=0 # [2] Build INDELible control files from R results
RUN_INDELIBLE=0 # [3] Build alignments
BUILD_CONTROL_FILES=1 # [4] Build control files for Partitiontest

R_LOG_FILE="rscript.log"
IND_LOG_FILE="indelible.log"
LOG_DIR="log"

INPUT_TREEFILE="scripts/sims/treefile.out"
INPUT_MODELFILE="scripts/sims/modelsfile.out"
INPUT_GEN2PARTFILE="scripts/sims/genestopartitions.out"
OUTPUT_DIR="indelible"
OUTPUT_ALIGNS_DIR="alignments"
OUTPUT_PARTEST_DIR="partitiontest"

if [ $RUN_R_SCRIPT -eq 1 ]; then
	# Generate models, partitions and trees
  echo "[1] Running R scripts"
  cd scripts
  R --vanilla -q < build_models.r 2>&1 > ../${LOG_DIR}/${R_LOG_FILE}
  cd -
else
  echo "[1] Running R scripts (AVOID)"
fi

NUM_SIMS=`wc -l ${INPUT_TREEFILE} | cut -d' ' -f1`
echo "    There are $((NUM_SIMS)) simulations"



if [ $BUILD_INDEL_FILES -eq 1 ]; then
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
	  echo [TYPE] AMINOACID 1 >> ${IND_FILE}
	  echo " " >> ${IND_FILE}

	  # Write tree
	  echo [TREE] tree `sed "${sim_index}q;d" ${INPUT_TREEFILE}` >> ${AUX_FILE}

	  part_header=`sed "${part_header_line}q;d" ${INPUT_GEN2PARTFILE}`
	  num_genes=`echo ${part_header} | cut -d' ' -f2`
	  num_partitions=`echo ${part_header} | cut -d' ' -f3`
	  echo " " >> ${AUX_FILE}

	  echo "    ${sim_index}/${NUM_SIMS} : ${num_partitions} partitions out of ${num_genes} genes"

	  # Write partitions
	  echo [PARTITIONS] partitions >> ${AUX_FILE}
	  for ((partition_index=1; partition_index<=${num_genes}; partition_index++)); do
	    partmodel=`sed "$((part_header_line + partition_index))q;d" ${INPUT_GEN2PARTFILE}`
	    model_index=`echo ${partmodel} | cut -d' ' -f2`
	    echo [tree model${model_index} ${PERGENE_LENGTH}] >> ${AUX_FILE}
	  done

	  #Write models
	  models=`sed -n -e $((part_header_line+1)),$((part_header_line+num_genes))p ${INPUT_GEN2PARTFILE} | cut -d' ' -f2`
	  models=`for each in ${models}; do echo $each; done | sort | uniq`
	  for each in ${models}; do
	    model=`sed -n -e ${each}p ${INPUT_MODELFILE}`
	    indeliblemodel=`echo $model | cut -d' ' -f1`
	    modelind=`echo $model | cut -d' ' -f2`
	    phymlmodel=`echo $model | cut -d' ' -f3`
	    partestmodel=`echo $model | cut -d' ' -f4`
	    isF=`echo $model | cut -d' ' -f5`
	    isI=`echo $model | cut -d' ' -f6`
	    isG=`echo $model | cut -d' ' -f7`
	    pInv=`echo $model | cut -d' ' -f8`
	    if [ $isG -eq 1 ]; then
	      ncat=4
	      shape=`echo $model | cut -d' ' -f9`
	    else
	      ncat=0
	      shape=0
	    fi
	    echo [MODEL] model${each} >> ${IND_FILE}
	    echo [submodel] ${indeliblemodel} >> ${IND_FILE}
	    echo [rates] ${pInv} ${shape} ${ncat} >> ${IND_FILE}
	    if [ $isF -eq 1 ]; then
	      echo [statefreq] 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 >> ${IND_FILE}
	    fi
	    echo " " >> ${IND_FILE}
	  done

	  cat ${AUX_FILE} >> ${IND_FILE}

	  echo " " >> ${IND_FILE}
	  echo [EVOLVE] partitions 1 alignment${sim_index} >> ${IND_FILE}

	  rm ${AUX_FILE}

	  part_header_line=$((part_header_line + num_genes + 1))
	done
else
  echo "[2] Building INDELible control files (AVOID)"
fi

if [ $RUN_INDELIBLE -eq 1 ]; then
  echo "[3] Building alignments with INDELible"
  cd ${OUTPUT_DIR}
  rm ${IND_LOG_FILE}

  # Loop over simulations
  for ((sim_index=1; sim_index<=${NUM_SIMS}; sim_index++)); do
    echo "    ${sim_index}/${NUM_SIMS} `date`"
    cp control${sim_index}.txt control.txt
    ./indelible >> ${IND_LOG_FILE}
    cat LOG.txt >> ${IND_LOG_FILE}
    mv alignment${sim_index}_TRUE.phy alignment${sim_index}.phy
  done

  echo "    Cleanup"
  rm *.fas LOG.txt trees.txt control.txt
  rm -rf ../${OUTPUT_ALIGNS_DIR}
  mkdir ../${OUTPUT_ALIGNS_DIR}
  mv alignment* ../${OUTPUT_ALIGNS_DIR}
  cd ..
else
  echo "[3] Building alignments with INDELible (AVOID)"
fi

if [ $BUILD_CONTROL_FILES -eq 1 ]; then
  echo "[4] Building PartitionTest control files"
  rm -rf ${OUTPUT_PARTEST_DIR}
  mkdir ${OUTPUT_PARTEST_DIR}

# Loop over simulations
  for ((sim_index=1; sim_index<=${NUM_SIMS}; sim_index++)); do
    align_filename=${OUTPUT_ALIGNS_DIR}/alignment${sim_index}.phy
    partest_filename=${OUTPUT_PARTEST_DIR}/alignment${sim_index}.conf
    num_genes=`head -n 1 ${align_filename} | tr -s ' ' | cut -d' ' -f2`
    num_genes=$((num_genes / 1000))
    echo "    ${sim_index}/${NUM_SIMS} : ${num_genes} genes"
    echo [PARTITIONS] > ${partest_filename}
    next_start=1
    for ((gene_index=1; gene_index<=${num_genes}; gene_index++)); do
      echo GENE${gene_index}=${next_start}-$((next_start+999)) >> ${partest_filename}
      next_start=$((next_start+1000))
    done
    echo [OUTPUT] >> ${partest_filename}
    echo models=models${sim_index}.out >> ${partest_filename}
    echo partitions=partitions${sim_index}.out >> ${partest_filename}
    echo schemes=schemes${sim_index}.out >> ${partest_filename}
  done

else
  echo "[4] Building PartitionTest control files (AVOID)"
fi

echo "[5] Done!"