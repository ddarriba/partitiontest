PERGENE_LENGTH=1000

R_LOG_FILE="rscript.log"
LOG_DIR="log"

INPUT_TREEFILE="scripts/sims/treefile.out"
INPUT_MODELFILE="scripts/sims/modelsfile.out"
INPUT_GEN2PARTFILE="scripts/sims/genestopartitions.out"
OUTPUT_DIR="indelible"

rm -rf ${LOG_DIR}
mkdir ${LOG_DIR}

rm -rf ${OUTPUT_DIR}
mkdir ${OUTPUT_DIR}

# Generate models, partitions and trees
cd scripts
R --vanilla -q < build_models.r 2>&1 > ../${LOG_DIR}/${R_LOG_FILE}
cd -

NUM_SIMS=`wc -l ${INPUT_TREEFILE} | cut -d' ' -f1`
echo There are $((NUM_SIMS)) simulations

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
  echo [TREE] tree `sed "${sim_index}q;d" ${INPUT_TREEFILE}` >> ${IND_FILE}

  part_header=`sed "${part_header_line}q;d" ${INPUT_GEN2PARTFILE}`
  num_genes=`echo ${part_header} | cut -d' ' -f2`
  num_partitions=`echo ${part_header} | cut -d' ' -f3`
  echo " " >> ${IND_FILE}

  echo Parsed ${num_genes} genes in ${num_partitions} partitions

  # Write partitions
  for ((partition_index=1; partition_index<=${num_genes}; partition_index++)); do
    partmodel=`sed "$((part_header_line + partition_index))q;d" ${INPUT_GEN2PARTFILE}`
    model_index=`echo ${partmodel} | cut -d' ' -f2`
    echo [PARTITION] partition${partition_index} >> ${IND_FILE}
    echo [tree model${model_index} ${PERGENE_LENGTH}] >> ${IND_FILE}
    echo " " >> ${IND_FILE}
  done

  #Write models
  echo BEGIN MODELS
  models=`sed -n -e $((part_header_line+1)),$((part_header_line+num_genes))p ${INPUT_GEN2PARTFILE} | cut -d' ' -f2`
  models=`for each in ${models}; do echo $each; done | sort | uniq`
  for each in ${models}; do
    model=`sed -n -e ${each}p ${INPUT_MODELFILE}`
    echo [MODEL] model${each} >> ${IND_FILE}
    echo [submodel] whatever >> ${IND_FILE}
    echo [pinv] whatever >> ${IND_FILE}
    echo " " >> ${IND_FILE}
  done

  part_header_line=$((part_header_line + num_genes + 1))
done
