R_LOG_FILE="rscript.log"
IND_LOG_FILE="indelible.log"


prefix=testbed
datatype=dna
NUM_SIMS=10
mintaxa=10
maxtaxa=50
mingenes=5
maxgenes=50
minlen=500
maxlen=1500

simsdir=sims.$prefix

LOG_DIR="${simsdir}/log"
INPUT_TREEFILE="${simsdir}/treefile.out"
INPUT_MODELFILE="${simsdir}/modelsfile.out"
INPUT_GEN2PARTFILE="${simsdir}/genestopartitions.out"
INPUT_PARTITIONSFILE="${simsdir}/partitionsfile.out"

OUTPUT_INDELIBLE_DIR="${simsdir}/indelible"
OUTPUT_ALIGNS_DIR="${simsdir}/alignments"
OUTPUT_PARTEST_DIR="${simsdir}/partitiontest"
OUTPUT_RKN_DIR="${simsdir}/rkn"
OUTPUT_RKT_DIR="${simsdir}/rkt"
OUTPUT_PART_SUMMARY="${simsdir}/partitionssummary.out"

# Generate models, partitions and trees
echo "[1] Running R scripts"
cd R-scripts
R --vanilla --slave --args $prefix $datatype $NUM_SIMS $mintaxa $maxtaxa $mingenes $maxgenes $minlen $maxlen < build_models.r 2>&1
cd -

rm -rf ${LOG_DIR} ${OUTPUT_INDELIBLE_DIR} ${OUTPUT_ALIGNS_DIR} ${OUTPUT_PARTEST_DIR} ${OUTPUT_RKN_DIR} ${OUTPUT_RKT_DIR}
mkdir  ${LOG_DIR} ${OUTPUT_INDELIBLE_DIR} ${OUTPUT_ALIGNS_DIR} ${OUTPUT_PARTEST_DIR} ${OUTPUT_RKN_DIR} ${OUTPUT_RKT_DIR}
cp bin/indelible ${OUTPUT_INDELIBLE_DIR}/

AUX_FILE=${OUTPUT_INDELIBLE_DIR}/aux.txt
part_header_line=1
# Loop over simulations
for ((sim_index=1; sim_index<=${NUM_SIMS}; sim_index++)); do

	echo "TRACE ${sim_index}/${NUM_SIMS}"
	# Indelible control file
	IND_FILE=${OUTPUT_INDELIBLE_DIR}/control${sim_index}.txt
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
		fA=`echo $model | cut -d' ' -f2`
		fC=`echo $model | cut -d' ' -f3`
		fG=`echo $model | cut -d' ' -f4`
		fT=`echo $model | cut -d' ' -f5`

		kappa=`echo $model | cut -d' ' -f6`

		rA=`echo $model | cut -d' ' -f7`
		rB=`echo $model | cut -d' ' -f8`
		rC=`echo $model | cut -d' ' -f9`
		rD=`echo $model | cut -d' ' -f10`
		rE=`echo $model | cut -d' ' -f11`
		rF=`echo $model | cut -d' ' -f12`

		pInv=0
		shape=`echo $model | cut -d' ' -f13`

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


	echo "TRACE    build raxml k1 file"
	unset trueparts
	declare -a trueparts
	next_pos=1
	for ((partition_index=1; partition_index<=${num_genes}; partition_index++)); do
		partline=`sed "$((part_header_line + partition_index))q;d" ${INPUT_GEN2PARTFILE}`
		cur_part=`echo ${partline} | cut -d' ' -f3`
		cur_len=`echo ${partline} | cut -d' ' -f4`
		if [ "" == "${trueparts[$cur_part]}" ]; then
			if [ "$datatype" == "dna" ]; then
				trueparts[$cur_part]="DNA, PART$cur_part=${next_pos}-$((next_pos + cur_len))"
			else
				trueparts[$cur_part]="AUTO, PART$cur_part=${next_pos}-$((next_pos + cur_len))"
			fi
		else
			trueparts[$cur_part]="${trueparts[$cur_part]},${next_pos}-$((next_pos + cur_len))"
		fi
		next_pos=$((next_pos + cur_len + 1))
	done

	for ((partition_index=1; partition_index<=${num_partitions}; partition_index++)); do
		echo ${trueparts[$partition_index]} >> $OUTPUT_RKT_DIR/control${sim_index}
	done

	part_header_line=$((part_header_line + num_genes + 1))

	echo "TRACE    build alignment"
	cd ${OUTPUT_INDELIBLE_DIR}
	cp control${sim_index}.txt control.txt
	./indelible >> ${IND_LOG_FILE}
	cat LOG.txt >> ${IND_LOG_FILE}
	mv alignment${sim_index}_TRUE.phy ../../${OUTPUT_ALIGNS_DIR}/alignment${sim_index}.phy
	rm *.fas LOG.txt trees.txt control.txt
	cd - > /dev/null

	align_filename=${OUTPUT_ALIGNS_DIR}/alignment${sim_index}.phy
	rkn_filename=${OUTPUT_RKN_DIR}/control${sim_index}
	partest_filename=${OUTPUT_PARTEST_DIR}/alignment${sim_index}.conf
	partition_line=`sed "$((sim_index * 3 - 2))q;d" $INPUT_PARTITIONSFILE`
	num_genes=`echo $partition_line | cut -d' ' -f 2`
	echo "TRACE    ${num_genes} genes"
	echo [PARTITIONS] > ${partest_filename}

	echo "TRACE    build raxml kn file and partest"
	next_start=1
	for ((gene_index=1; gene_index<=${num_genes}; gene_index++)); do
		partmodel=`sed "$((part_header_line + gene_index))q;d" ${INPUT_GEN2PARTFILE}`
		gene_length=`echo ${partmodel} | cut -d' ' -f4`
		echo GENE${gene_index}=${next_start}-$((next_start+gene_length-1)) >> ${partest_filename}
		if [ "$datatype" == "dna" ]; then
		echo "DNA, GENE${gene_index}=${next_start}-$((next_start+gene_length-1))" >> ${rkn_filename}
		else
		echo "AUTO, GENE${gene_index}=${next_start}-$((next_start+gene_length-1))" >> ${rkn_filename}
		fi
		next_start=$((next_start+gene_length))
	done
	echo [OUTPUT] >> ${partest_filename}
	echo models=models${sim_index}.out >> ${partest_filename}
	echo partitions=partitions${sim_index}.out >> ${partest_filename}
	echo schemes=schemes${sim_index}.out >> ${partest_filename}
done
