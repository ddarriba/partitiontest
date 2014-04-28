/*
 * ConfigParser.cpp
 *
 *  Created on: May 1, 2013
 *      Author: diego
 */

#include "ConfigParser.h"
#include "util/Utilities.h"
#include "IniParser.h"
#include <string.h>
#include <assert.h>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <parsePartition.h>

namespace partest {

struct comparePartitionInfos {
	inline bool operator()(partitionInfo p1, partitionInfo p2) {

		if (p1.numberOfSections < 1)
			return true;
		if (p2.numberOfSections < 1)
			return false;
		if (p1.start[0] != p2.start[0])
			return p1.start[0] < p2.start[0];
		else
			return p1.stride[0] < p2.stride[0];
	}
};

ConfigParser::ConfigParser(const char * configFile) {
	maxSamples = 1;

	dataType = DT_DEFAULT;
	optimizeMode = OPT_SEARCH;

	vector<const char *> * keys;

	if (configFile != 0 && strcmp(configFile, "")) {

		int partitionId = 0;
		const char * value;

#ifdef SIMPLEINI
		CSimpleIniA ini;
		ini.SetUnicode();
		SI_Error rc = ini.LoadFile(configFile);
		if (rc < 0)
		exit_partest(EX_IOERR);

		CSimpleIniA::TNamesDepend keys;
#else
		IniParser ini(configFile);
#endif
		/** INPUT **/
		keys = ini.GetAllKeys(INPUT_TAG);
		if (keys->size() > 0) {
			value = ini.GetValue(INPUT_TAG, INPUT_MSA_TAG, "");
			if (strcmp(value, "")) {
				strcpy(input_file, value);
			}
			value = ini.GetValue(INPUT_TAG, INPUT_TREE_TAG, "");
			if (strcmp(value, "")) {
				strcpy(user_tree, value);
			}
			value = ini.GetValue(INPUT_TAG, INPUT_DATATYPE_TAG, "nt");
			if (!strcmp(value, "aa")) {
				dataType = DT_PROTEIC;
			} else if (!strcmp(value, "nt")) {
				dataType = DT_NUCLEIC;
			} else {
				dataType = DT_NUCLEIC;
			}
		}

		/** CANDIDATE MODELS **/
#ifdef SIMPLEINI
		ini.GetAllKeys(MODELS_TAG, keys);
		if (keys.size() > 0) {
#else
#endif
		value = ini.GetValue(MODELS_TAG, MODELS_INCLUDE_TAG, "all");
		if (!strcmp(value, "all")) {
			optimizeMode = OPT_SEARCH;
		} else if (!strcmp(value, "gtr")) {
			optimizeMode = OPT_GTR;
		} else {
			// parse protein matrices
			protModels = 0;
			optimizeMode = OPT_CUSTOM;
			istringstream iss(value);
			string curMatrix;
			while (iss) {
				unsigned int curValue = 0;
				iss >> curMatrix;
				if (!strcmp(curMatrix.c_str(), "dayhoff")) {
					curValue = PROT_MATRIX_DAYHOFF;
				} else if (!strcmp(curMatrix.c_str(), "dcmut")) {
					curValue = PROT_MATRIX_DCMUT;
				} else if (!strcmp(curMatrix.c_str(), "jtt")) {
					curValue = PROT_MATRIX_JTT;
				} else if (!strcmp(curMatrix.c_str(), "mtrev")) {
					curValue = PROT_MATRIX_MTREV;
				} else if (!strcmp(curMatrix.c_str(), "wag")) {
					curValue = PROT_MATRIX_WAG;
				} else if (!strcmp(curMatrix.c_str(), "cprev")) {
					curValue = PROT_MATRIX_CPREV;
				} else if (!strcmp(curMatrix.c_str(), "rtrev")) {
					curValue = PROT_MATRIX_RTREV;
				} else if (!strcmp(curMatrix.c_str(), "vt")) {
					curValue = PROT_MATRIX_VT;
				} else if (!strcmp(curMatrix.c_str(), "blosum62")) {
					curValue = PROT_MATRIX_BLOSUM62;
				} else if (!strcmp(curMatrix.c_str(), "mtmam")) {
					curValue = PROT_MATRIX_MTMAM;
				} else if (!strcmp(curMatrix.c_str(), "mtart")) {
					curValue = PROT_MATRIX_MTART;
				} else if (!strcmp(curMatrix.c_str(), "hivb")) {
					curValue = PROT_MATRIX_HIVB;
				} else if (!strcmp(curMatrix.c_str(), "hivw")) {
					curValue = PROT_MATRIX_HIVW;
				} else if (!strcmp(curMatrix.c_str(), "mtzoa")) {
					curValue = PROT_MATRIX_MTZOA;
				} else if (!strcmp(curMatrix.c_str(), "pmb")) {
					curValue = PROT_MATRIX_PMB;
				} else if (!strcmp(curMatrix.c_str(), "flu")) {
					curValue = PROT_MATRIX_FLU;
				} else if (!strcmp(curMatrix.c_str(), "auto")) {
					curValue = PROT_MATRIX_AUTO;
				}
				protModels |= Utilities::binaryPow(curValue);
			}
		}

		/** SEARCH ALGORITHM **/
		keys = ini.GetAllKeys(SEARCH_TAG);
		searchData = keys->size() > 0;
		if (searchData) {
			value = ini.GetValue(SEARCH_TAG, SEARCH_ALGORITHM_TAG, "none");
			if (!strcmp(value, "greedy")) {
				searchAlgorithm = SearchGreedy;
			} else if (!strcmp(value, "exhaustive")) {
				searchAlgorithm = SearchExhaustive;
			} else if (!strcmp(value, "random")) {
				searchAlgorithm = SearchRandom;
			} else if (!strcmp(value, "hcluster")) {
				searchAlgorithm = SearchHCluster;
			} else {
				cerr << "Invalid search algorithm : " << value << endl;
				exit_partest(EX_SOFTWARE);
			}

			value = ini.GetValue(SEARCH_TAG, SEARCH_ALGORITHM_REPS, "1");\
			for (int i = 0; value[i] != 0; i++)
				if (!isdigit(value[i])) {
					cerr << "[ERROR] Number of replicates \"" << value
							<< "\" must be an integer." << endl;
					exit_partest(EX_USAGE);
				}
			maxSamples = atoi(value);
		}

		/** PARTITIONS **/
		keys = ini.GetAllKeys(PARTITIONS_TAG);
		number_of_genes = keys->size();

		partitions = new vector<partitionInfo>(number_of_genes);
		singleGeneNames = (string **) malloc(number_of_genes * sizeof(string*));

		char * lineBuffer = (char *) malloc(150);
		for (const char * key : *keys) {
			partitions->at(partitionId).name = key;
			singleGeneNames[partitionId] = new string(key);

			strcpy(lineBuffer,
					ini.GetValue(PARTITIONS_TAG, key, "default"));
			parsePartitionDetails(lineBuffer, &partitions->at(partitionId));
			partitionId++;
		}
		free(lineBuffer);

		std::sort(partitions->begin(), partitions->end(),
				comparePartitionInfos());

		pllPartitionRegion * pregion;
		pllPartitionInfo * pinfo;

		pllQueueInit(&pllPartsQueue);

		for (unsigned int i = 0; i < number_of_genes; i++) {
			partitions->at(i).partitionId.push_back(i);
			pinfo = (pllPartitionInfo *) malloc(sizeof(pllPartitionInfo));
			pllQueueInit(&(pinfo->regionList));
			pllQueueAppend(pllPartsQueue, (void *) pinfo);

			pinfo->partitionName = (char *) malloc(
					(partitions->at(i).name.size() + 1) * sizeof(char));
			strcpy(pinfo->partitionName, partitions->at(i).name.c_str());
			pinfo->partitionModel = (char *) malloc(1);

			pinfo->protModels = -1;
			pinfo->protFreqs = -1;
			pinfo->dataType = PLL_DNA_DATA;
			pinfo->optimizeBaseFrequencies = PLL_TRUE;
			for (int j = 0; j < partitions->at(i).numberOfSections; j++) {
				pregion = (pllPartitionRegion *) malloc(
						sizeof(pllPartitionRegion));
				pregion->start = partitions->at(i).start[j];
				pregion->end = partitions->at(i).end[j];
				pregion->stride = partitions->at(i).stride[j];
				pllQueueAppend(pinfo->regionList, (void *) pregion);
			}
		}

		/** OUTPUT **/
		value = ini.GetValue(OUTPUT_TAG, OUTPUT_BASE_PATH, 0);
		if (value) {
			if (value[strlen(value) - 1] != char_separator)
				outputBasePath = string(value) + os_separator;
			else
				outputBasePath = string(value);
		}
		value = ini.GetValue(OUTPUT_TAG, OUTPUT_TMP_PATH, 0);
		if (value) {
			if (value[strlen(value) - 1] != char_separator)
				outputTmpPath = string(value) + os_separator;
			else
				outputTmpPath = string(value);
		} else {
			outputTmpPath = outputBasePath;
		}
		value = ini.GetValue(OUTPUT_TAG, OUTPUT_RESULTS_TAG, 0);
		if (value) {
			outputFileResults = outputBasePath + string(value);
		}
		value = ini.GetValue(OUTPUT_TAG, OUTPUT_MODELS_TAG, 0);
		if (value) {
			outputFileModels = outputBasePath + string(value);
		}
		value = ini.GetValue(OUTPUT_TAG, OUTPUT_PARTS_TAG, 0);
		if (value) {
			outputFilePartitions = outputBasePath + string(value);
		}
		value = ini.GetValue(OUTPUT_TAG, OUTPUT_SCHEMES_TAG, 0);
		if (value) {
			outputFileSchemes = outputBasePath + string(value);
		}

		if ((optimizeMode == OPT_CUSTOM) & (dataType == DT_NUCLEIC)) {
			cerr
					<< "[ERROR] Custom model set is not available for nucleic data."
					<< endl;
			cerr << "        Please use \"gtr\" or \"all\" option." << endl;
			exit_partest(EX_USAGE);
		} else if ((optimizeMode == OPT_GTR) & (dataType == DT_PROTEIC)) {
			cerr << "[ERROR] GTR option is not available for proteic data."
					<< endl;
			cerr
					<< "        Please use \"all\" option or define the set of models."
					<< endl;
			exit_partest(EX_USAGE);
		}
	}

}

int ConfigParser::parsePartitionDetails(char * line,
		struct partitionInfo * pInfo) {
	int numberOfSections = 0;
	char * parsed = strtok(line, "-");

	while (parsed != NULL) {
		int start = atoi(parsed);
		parsed = strtok(NULL, ",");
		int end = atoi(parsed);
		pInfo->start[numberOfSections] = start;
		pInfo->end[numberOfSections] = end;
		pInfo->stride[numberOfSections] = 1;

		numberOfSections++;
		parsed = strtok(NULL, "-");
	}
	pInfo->numberOfSections = numberOfSections;

	return 0;
}

int ConfigParser::parsePartitionLine(char * line,
		struct partitionInfo * pInfo) {
	int numberOfSections = 0;
	char * parsed = strtok(line, "-");

	while (parsed != NULL) {
		int start = atoi(parsed);
		parsed = strtok(NULL, ",");
		int end = atoi(parsed);
		pInfo->start[numberOfSections] = start;
		pInfo->end[numberOfSections] = end;
		pInfo->stride[numberOfSections] = 1;

		numberOfSections++;
		parsed = strtok(NULL, ",");
	}
	pInfo->numberOfSections = numberOfSections;

	return 0;
}

ConfigParser::~ConfigParser() {
	if (partitions) {
		delete partitions;
	}
}

struct partitionInfo ConfigParser::getPartition(int index) {
	if (!partitions) {
		cerr << "ERROR: No partitions were defined" << endl;
		exit_partest(EX_SOFTWARE);
	}
	return partitions->at(index);

}

void ConfigParser::printFormat() {
	cout << "Config file format:" << endl << endl;
	cout << "   ;THIS IS A COMMENT" << endl;
	cout << endl << "   ; Start of input data" << endl;
	cout << "   [" << INPUT_TAG << "]" << endl;
	cout << "   " << INPUT_MSA_TAG << "=INPUT_ALIGNMENT_FILE" << endl;
	cout << "   " << INPUT_TREE_TAG << "=INPUT_TREE_FILE" << endl;
	cout << "   " << INPUT_DATATYPE_TAG << "={nt|aa} (default: nt)" << endl;
	cout << endl << "   ; Start of searching options" << endl;
	cout << "   [" << SEARCH_TAG << "]" << endl;
	cout << "   " << SEARCH_ALGORITHM_TAG
			<< "={greedy|hcluster} (default: greedy)";
	cout << "   " << SEARCH_ALGORITHM_REPS << "=# (default: 1)" << endl;
	cout << endl << "   ; Start of candidate models description" << endl;
	cout << "   [" << MODELS_TAG << "]" << endl;
	cout << "   " << MODELS_INCLUDE_TAG << "={all | gtr | [LIST]} (default:all)"
			<< endl;
	cout << endl << "   ; Start of partitions" << endl;
	cout << "   [" << PARTITIONS_TAG << "]" << endl;
	cout << "   PART1=INI1-END1" << endl;
	cout << "   PART1=INI2-END2" << endl;
	cout << "   ..." << endl;
	cout << "   PART1=INI3-END3" << endl;
	cout << endl << "   ; Start of output section" << endl;
	cout << "   [" << OUTPUT_TAG << "]" << endl;
	cout << "   " << OUTPUT_BASE_PATH
			<< "=OUTPUT_BASE_URL (default:current directory)" << endl;
	cout << "   " << OUTPUT_MODELS_TAG << "=OUTPUT_MODELS_FILE" << endl;
	cout << "   " << OUTPUT_PARTS_TAG << "=OUTPUT_PARTITIONS_FILE" << endl;
	cout << "   " << OUTPUT_SCHEMES_TAG << "=OUTPUT_SCHEMES_FILE" << endl;
	cout << "   " << OUTPUT_RESULTS_TAG << "=OUTPUT_RESULTS_FILE" << endl;
	cout << endl << "Example:" << endl << endl;
	cout << "   [" << PARTITIONS_TAG << "]" << endl;
	cout << "   DNA1=1-976" << endl;
	cout << "   DNA2=976-1803" << endl;
	cout << "   DNA3=1804-2700" << endl;
	cout << "   DNA4=2701-3815" << endl;
	cout << "   ; End of partitions for file.phy" << endl << endl;
	cout << INPUT_TAG << "/" << INPUT_DATATYPE_TAG << endl;
	cout << "      nt - Nucleic (DNA) data" << endl;
	cout << "      aa - Amino Acid (protein) data" << endl;
	cout << SEARCH_TAG << "/" << SEARCH_ALGORITHM_TAG << endl;
	cout << "      greedy   - Greedy search" << endl;
	cout << "      hcluster - Hierarchical Clustering search" << endl;
	cout << SEARCH_TAG << "/" << SEARCH_ALGORITHM_REPS << endl;
	cout << "      (int) # - Maximum number of replicates on each HCluster step"
			<< endl;
	cout << MODELS_TAG << "/" << MODELS_INCLUDE_TAG << endl;
	cout << "      all    - Evaluate the whole set of models" << endl;
	cout << "      gtr    - Evaluate only gtr models (Only for DNA data)"
			<< endl;
	cout
			<< "      [LIST] - List of matrices to evaluate (Only for protein data)"
			<< endl;
	cout
			<< "               { dayhoff, dcmut, jtt, mtrev, cprev, rtrev, wag, vt,"
			<< endl;
	cout
			<< "                 blossum62, mtmam, mtart, hivb, hivw, mtzoa, pmb, flu }"
			<< endl;
	cout << "               e.g.,  " << MODELS_INCLUDE_TAG
			<< "=dayhoff dcmut jtt" << endl;
	cout << OUTPUT_TAG << endl;
	cout
			<< "      Define urls for the output files. set to N/A for avoid output"
			<< endl;
	cout << "               e.g., " << OUTPUT_MODELS_TAG << "=N/A" << endl;
}

void ConfigParser::createTemplate() {
	cout << ";THIS IS A COMMENT" << endl;
	cout << "; Start of input data" << endl;
	cout << "[" << INPUT_TAG << "]" << endl;
	cout << INPUT_MSA_TAG << "=input.phy" << endl;
	cout << INPUT_TREE_TAG << "=input.tree" << endl;
	cout << INPUT_DATATYPE_TAG << "=nt" << endl;
	cout << "; Start of candidate models description" << endl;
	cout << "[" << MODELS_TAG << "]" << endl;
	cout << MODELS_INCLUDE_TAG << "=all" << endl;
	cout << "; Start of searching options" << endl;
	cout << "[" << SEARCH_TAG << "]" << endl;
	cout << SEARCH_ALGORITHM_TAG << "=hcluster" << endl;
	cout << SEARCH_ALGORITHM_REPS << "=1" << endl;
	cout << "[" << PARTITIONS_TAG << "]" << endl;
	cout << "PART1=INI1-END1" << endl;
	cout << "PART2=INI2-END2" << endl;
	cout << "..." << endl;
	cout << "PARTn=INIn-ENDn" << endl;
}

vector<partitionInfo> * ConfigParser::getPartitions() {
	return partitions;
}

} /* namespace partest */
