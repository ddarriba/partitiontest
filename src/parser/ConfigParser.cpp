/*  PartitionTest, fast selection of the best fit partitioning scheme for
 *  multi-gene data sets.
 *  Copyright May 2013 by Diego Darriba
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 3 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 *
 *  For any other inquiries send an Email to Diego Darriba
 *  ddarriba@udc.es
 */

/**
 * @file ConfigParser.cpp
 * @author Diego Darriba
 */

#include "ConfigParser.h"
#include "util/Utilities.h"
#include "INIReader.h"
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

void ConfigParser::createSinglePartition() {
	if (!input_file) {
		cerr << "[ERROR] Input file has not been set" << endl;
		exit_partest(EX_CONFIG);
	}

	number_of_genes = 1;
	partitions = new vector<partitionInfo>(number_of_genes);
	singleGeneNames = (string **) malloc(sizeof(string*));
	singleGeneNames[0] = new string("SinglePartition");

	phylip = pllParseAlignmentFile(PLL_FORMAT_PHYLIP, input_file->c_str());
	if (!phylip) {
		cerr << "[ERROR] There was an error parsing input data." << endl;
		exit_partest(EX_IOERR);
	}
	num_taxa = phylip->sequenceCount;
	seq_len = phylip->sequenceLength;

	pllPartitionRegion * pregion;
	pllPartitionInfo * pinfo;
	pllQueueInit(&pllPartsQueue);

	partitions->at(0).partitionId.push_back(0);
	pinfo = (pllPartitionInfo *) malloc(sizeof(pllPartitionInfo));
	pinfo->ascBias = PLL_FALSE;
	pllQueueInit(&(pinfo->regionList));
	pllQueueAppend(pllPartsQueue, (void *) pinfo);
	pinfo->partitionName = (char *) malloc(2);
	strcpy(pinfo->partitionName, "S");
	pinfo->partitionModel = (char *) malloc(1);
	if (starting_topology == StartTopoDEFAULT) {
		starting_topology = DEFAULT_STARTING_TOPOLOGY;
	}
	if (data_type == DT_DEFAULT) {
		data_type = dataType == DT_DEFAULT ? DEFAULT_DATA_TYPE : dataType;
	}
	switch (data_type) {
	case DT_NUCLEIC:
		pinfo->protModels = -1;
		pinfo->protUseEmpiricalFreqs = -1;
		pinfo->dataType = PLL_DNA_DATA;
		pinfo->optimizeBaseFrequencies = PLL_TRUE;
		break;
	case DT_PROTEIC:
		pinfo->protModels = PLL_AUTO;
		pinfo->protUseEmpiricalFreqs = PLL_FALSE;
		pinfo->dataType = PLL_AA_DATA;
		pinfo->optimizeBaseFrequencies = PLL_TRUE;
		break;
	case DT_DEFAULT:
		exit_partest(EX_SOFTWARE);
		break;
	}
	pregion = (pllPartitionRegion *) malloc(sizeof(pllPartitionRegion));
	pregion->start = 1;
	pregion->end = seq_len;
	pregion->stride = 1;
	pllQueueAppend(pinfo->regionList, (void *) pregion);
}

ConfigParser::ConfigParser(const char * configFile) {
	maxSamples = 1;

	dataType = DT_DEFAULT;
	optimizeMode = OPT_SEARCH;

	if (configFile != 0 && strcmp(configFile, "")) {

		const char * value;

		//IniParser ini(configFile);
		INIReader ini(configFile);
		/** INPUT **/
		value = ini.Get(INPUT_TAG, INPUT_MSA_TAG, "").c_str();
		if (strcmp(value, "")) {
			int prefixLen = 0;
			char basedir[strlen(configFile)];
			strcpy(basedir, configFile);
			char * lastSlash = strrchr(basedir, char_separator);
			if (lastSlash) {
				*(lastSlash + 1) = 0;
				strcpy(inputFile, basedir);
				prefixLen = strlen(basedir);
			}
			strcpy(inputFile + prefixLen, value);
		}

		if (starting_topology == StartTopoDEFAULT)
		value = ini.Get(INPUT_TAG, INPUT_STARTTOPO_TAG, "").c_str();
			if (!strcmp(value, "fixed")) {
				starting_topology = StartTopoFIXED;
			} else if (!strcmp(value, "ml")) {
				starting_topology = StartTopoML;
			} else if (!strcmp(value, "mp")) {
				starting_topology = StartTopoMP;
			} else if (!strcmp(value, "user")) {
				starting_topology = StartTopoUSER;
			} else {
				starting_topology = DEFAULT_STARTING_TOPOLOGY;
			}

		value = ini.Get(INPUT_TAG, INPUT_TREE_TAG, "").c_str();
		if (strcmp(value, "")) {
			strcpy(userTree, value);
		}

		value = ini.Get(INPUT_TAG, INPUT_DATATYPE_TAG, "nt").c_str();
		if (!strcmp(value, "aa")) {
			dataType = DT_PROTEIC;
		} else if (!strcmp(value, "nt")) {
			dataType = DT_NUCLEIC;
		} else {
			dataType = DT_NUCLEIC;
		}

		value = ini.Get(INPUT_TAG, INPUT_KEEPBRANCHES_TAG, "false").c_str();
		if (!strcmp(value, "true") && reoptimize_branch_lengths == true) {
			reoptimize_branch_lengths = false;
		}

		/** CANDIDATE MODELS **/
		value = ini.Get(MODELS_TAG, MODELS_INCLUDE_TAG, "all").c_str();
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
				size_t curValue = 0;
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
		if (search_algo == SearchDefault) {
			value = ini.Get(SEARCH_TAG, SEARCH_ALGORITHM_TAG, "auto").c_str();
			if (!strcmp(value, "greedy")) {
				searchAlgorithm = SearchGreedy;
			} else if (!strcmp(value, "exhaustive")) {
				searchAlgorithm = SearchExhaustive;
			} else if (!strcmp(value, "random")) {
				searchAlgorithm = SearchRandom;
			} else if (!strcmp(value, "hcluster")) {
				searchAlgorithm = SearchHCluster;
			} else if (!strcmp(value, "auto")) {
				searchAlgorithm = SearchAuto;
			} else {
				cerr << "Invalid search algorithm : " << value << endl;
				exit_partest(EX_CONFIG);
			}
		}

		maxSamples = ini.GetInteger(SEARCH_TAG, SEARCH_ALGORITHM_REPS, 1);

		/** PARTITIONS **/
		std::vector<std::string> * keys = ini.getGenes(
		PARTITIONS_TAG);
		number_of_genes = keys->size()/2;

		if (data_type == DT_DEFAULT) {
			data_type = (dataType == DT_DEFAULT) ? DEFAULT_DATA_TYPE : dataType;
		}

		if (number_of_genes) {
			partitions = new vector<partitionInfo>(number_of_genes);
			singleGeneNames = (string **) malloc(
					number_of_genes * sizeof(string*));

			for (size_t partitionId = 0; partitionId < number_of_genes; partitionId++) {
				partitions->at(partitionId).name = keys->at(partitionId*2);
				singleGeneNames[partitionId] = new string(keys->at(partitionId*2));

				char lineBuffer[keys->at(partitionId * 2 + 1).length() + 1];
				strcpy(lineBuffer, keys->at(partitionId * 2 + 1).c_str());

				parsePartitionDetails(lineBuffer, &partitions->at(partitionId));
			}
			delete keys;

			std::sort(partitions->begin(), partitions->end(),
					comparePartitionInfos());

			pllPartitionRegion * pregion;
			pllPartitionInfo * pinfo;

			pllQueueInit(&pllPartsQueue);

			for (size_t i = 0; i < number_of_genes; i++) {
				partitions->at(i).partitionId.push_back(i);
				pinfo = (pllPartitionInfo *) malloc(sizeof(pllPartitionInfo));
				pinfo->ascBias = PLL_FALSE;
				pllQueueInit(&(pinfo->regionList));
				pllQueueAppend(pllPartsQueue, (void *) pinfo);

				pinfo->partitionName = (char *) malloc(
						(partitions->at(i).name.size() + 1) * sizeof(char));
				strcpy(pinfo->partitionName, partitions->at(i).name.c_str());
				pinfo->partitionModel = (char *) malloc(1);

				switch (data_type) {
				case DT_NUCLEIC:
					pinfo->protModels = -1;
					pinfo->protUseEmpiricalFreqs = -1;
					pinfo->dataType = PLL_DNA_DATA;
					pinfo->optimizeBaseFrequencies = PLL_TRUE;
					break;
				case DT_PROTEIC:
					pinfo->protModels = PLL_AUTO;
					pinfo->protUseEmpiricalFreqs = PLL_FALSE;
					pinfo->dataType = PLL_AA_DATA;
					pinfo->optimizeBaseFrequencies = PLL_TRUE;
					break;
				case DT_DEFAULT:
					exit_partest(EX_SOFTWARE);
					break;
				}

				for (int j = 0; j < partitions->at(i).numberOfSections; j++) {
					pregion = (pllPartitionRegion *) malloc(
							sizeof(pllPartitionRegion));
					pregion->start = partitions->at(i).start[j];
					pregion->end = partitions->at(i).end[j];
					pregion->stride = partitions->at(i).stride[j];
					pllQueueAppend(pinfo->regionList, (void *) pregion);
				}
			}
		} else {
			createSinglePartition();
		}

		/** SCHEMES **/
		vector<string> * defSchemes = ini.getSchemes(SCHEMES_TAG);
		number_of_schemes = defSchemes->size();
		if (number_of_schemes) {
			schemes = new vector<t_partitioningScheme>(number_of_schemes);

			int schemeId = 0;
			for (size_t i = 0; i < defSchemes->size(); i++) {
				string scheme = defSchemes->at(i);
				char lineBuffer[scheme.length() + 1];
				strcpy(lineBuffer, scheme.c_str());
				parseScheme(lineBuffer, &(schemes->at(schemeId)));
				schemeId++;
			}
		}
		delete defSchemes;

		/** OUTPUT **/
		value = ini.Get(OUTPUT_TAG, OUTPUT_BASE_PATH, "").c_str();
		if (strcmp(value, "")) {
			if (value[strlen(value) - 1] != char_separator)
				outputBasePath = string(value) + os_separator;
			else
				outputBasePath = string(value);
		}
		value = ini.Get(OUTPUT_TAG, OUTPUT_TMP_PATH, "").c_str();
		if (strcmp(value, "")) {
			if (value[strlen(value) - 1] != char_separator)
				outputTmpPath = string(value) + os_separator;
			else
				outputTmpPath = string(value);
		} else {
			outputTmpPath = outputBasePath;
		}
		value = ini.Get(OUTPUT_TAG, OUTPUT_RESULTS_TAG, "").c_str();
		if (strcmp(value, "")) {
			outputFileResults = outputBasePath + string(value);
		}
		value = ini.Get(OUTPUT_TAG, OUTPUT_MODELS_TAG, "").c_str();
		if (strcmp(value, "")) {
			outputFileModels = outputBasePath + string(value);
		}
		value = ini.Get(OUTPUT_TAG, OUTPUT_PARTS_TAG, "").c_str();
		if (strcmp(value, "")) {
			outputFilePartitions = outputBasePath + string(value);
		}
		value = ini.Get(OUTPUT_TAG, OUTPUT_SCHEMES_TAG, "").c_str();
		if (strcmp(value, "")) {
			outputFileSchemes = outputBasePath + string(value);
		}

		if ((optimizeMode == OPT_CUSTOM) & (dataType == DT_NUCLEIC)) {
			cerr
					<< "[ERROR] Custom model set is not available for nucleic data."
					<< endl;
			cerr << "        Please use \"gtr\" or \"all\" option." << endl;
			exit_partest(EX_CONFIG);
		} else if ((optimizeMode == OPT_GTR) & (dataType == DT_PROTEIC)) {
			cerr << "[ERROR] GTR option is not available for proteic data."
					<< endl;
			cerr
					<< "        Please use \"all\" option or define the set of models."
					<< endl;
			exit_partest(EX_CONFIG);
		}
	} else {
		// no configuration file provided. A single partition will be used
		if (search_algo == SearchDefault) {
			search_algo = SearchAuto;
		}
		createSinglePartition();
	}
	if (ic_type == IC_DEFAULT) {
		ic_type = DEFAULT_IC_TYPE;
	}
}

int ConfigParser::parsePartitionDetails(char * line,
		struct partitionInfo * pInfo) {

	pllQueue * partitionLine;
	pllPartitionInfo * pi;
	char the_line[strlen(line) + 10];
	if (data_type == DT_NUCLEIC) {
		strcpy(the_line, "DNA, P=");
	} else {
		strcpy(the_line, "AUTO, P=");
	}
	strcat(the_line, line);
	the_line[strlen(the_line)] = '\0';
	partitionLine = pllPartitionParseString(the_line);
	if (!partitionLine) {
		cerr << "[ERROR] Could not read partition:" << endl;
		cerr << "-> " << line << endl;
		exit_partest(EX_IOERR);
	}
	pi = (pllPartitionInfo *) partitionLine->head->item;
	int numberOfSections = 0;
	for (pllQueueItem * qItem = pi->regionList->head; qItem; qItem =
			qItem->next) {
		pllPartitionRegion * region = (pllPartitionRegion *) qItem->item;
		pInfo->start[numberOfSections] = region->start;
		pInfo->end[numberOfSections] = region->end;
		pInfo->stride[numberOfSections] = region->stride;
		numberOfSections++;
	}
	pInfo->numberOfSections = numberOfSections;

	return 0;
}

int ConfigParser::parseScheme(char * line, t_partitioningScheme * scheme) {

	char * parsed = strtok(line, "(");
	char * rest = strtok(NULL, "\0");
	parsed = strtok(parsed, ")");
	while (parsed != NULL) {
		char * parsedPart;
		parsedPart = strtok(parsed, ",");
		t_partitionElementId nextPart;
		while (parsedPart != NULL) {
			Utilities::toLower(parsedPart);
			size_t nextSingleElement = number_of_genes;
			for (size_t i = 0; i < partitions->size(); i++) {
				partitionInfo pInfo = partitions->at(i);
				if (!pInfo.name.compare(parsedPart)) {
					nextSingleElement = pInfo.partitionId.at(0);
					break;
				}
			}
			if (nextSingleElement >= number_of_genes) {
				cerr << "[ERROR] Partition " << parsedPart << " not found"
						<< endl;
				exit_partest(EX_IOERR);
			}
			nextPart.push_back(nextSingleElement);
			parsedPart = strtok(NULL, ",");
		}
		scheme->push_back(nextPart);

		parsed = strtok(rest, "(");
		rest = strtok(NULL, "\0");
		parsed = strtok(parsed, ")");
	}

	return 0;
}

ConfigParser::~ConfigParser() {
	if (partitions) {
		delete partitions;
	}
}

void ConfigParser::printFormat() {
	cout << "Config file format:" << endl << endl;
	cout << "   ;THIS IS A COMMENT" << endl;
	cout << endl << "   ; Start of input data" << endl;
	cout << "   [" << INPUT_TAG << "]" << endl;
	cout << "   " << INPUT_MSA_TAG << "=INPUT_ALIGNMENT_FILE" << endl;
	cout << "   " << INPUT_TREE_TAG << "=INPUT_TREE_FILE" << endl;
	cout << "   " << INPUT_DATATYPE_TAG << "={nt|aa} (default: nt)" << endl;
	cout << "   " << INPUT_KEEPBRANCHES_TAG << "={true|false} (default: false)" << endl;
	cout << endl << "   ; Start of searching options" << endl;
	cout << "   [" << SEARCH_TAG << "]" << endl;
	cout << "   " << SEARCH_ALGORITHM_TAG
			<< "={greedy|hcluster|random|auto} (default: auto)";
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
	cout << "   PARTn=INIn-ENDn" << endl;
	cout << endl << "   ; Start of schemes" << endl;
	cout << "   [" << SCHEMES_TAG << "]" << endl;
	cout << "   S1=(GENE1,GENE2)(GENE3)..." << endl;
	cout << "   S2=(GENE1,GENE2,GENE3)..." << endl;
	cout << "   ..." << endl;
	cout << "   Sn=(GENE1)(GENE2,GENE3,...)" << endl;
	cout << endl << "   ; Start of output section" << endl;
	cout << "   [" << OUTPUT_TAG << "]" << endl;
	cout << "   " << OUTPUT_BASE_PATH
			<< "=OUTPUT_BASE_URL (default:partitiontest_FILENAME" << endl;
	cout << "   " << OUTPUT_MODELS_TAG
			<< "=OUTPUT_MODELS_FILE (default: $path/models)" << endl;
	cout << "   " << OUTPUT_SCHEMES_TAG
			<< "=OUTPUT_SCHEMES_FILE (default: $path/schemes)" << endl;
	cout << "   " << OUTPUT_RESULTS_TAG
			<< "=OUTPUT_RESULTS_FILE (default: $path/results)" << endl;
	cout << endl << "Example:" << endl << endl;
	cout << "   [" << PARTITIONS_TAG << "]" << endl;
	cout << "   DNA1=1-976" << endl;
	cout << "   DNA2=976-1803" << endl;
	cout << "   DNA3=1804-2700" << endl;
	cout << "   DNA4=2701-3815" << endl << endl;
	cout << INPUT_TAG << "/" << INPUT_DATATYPE_TAG << endl;
	cout << "      nt - Nucleic (DNA) data" << endl;
	cout << "      aa - Amino Acid (protein) data" << endl;
	cout << SEARCH_TAG << "/" << SEARCH_ALGORITHM_TAG << endl;
	cout << "      greedy   - Greedy search" << endl;
	cout << "      hcluster - Hierarchical Clustering search" << endl;
	cout << "      random   - Random search" << endl;
	cout << "      auto     - Auto selected search algorithm" << endl;
	cout << SEARCH_TAG << "/" << SEARCH_ALGORITHM_REPS << endl;
	cout
			<< "      (int) # - Maximum number of replicates on each HCluster/Random step"
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
	cout << SCHEMES_TAG << endl;
	cout << "      Define a set of schemes to evaluate instead of searching."
			<< endl;
	cout << "               e.g., "
			<< "S1=(GENE1,GENE2)(GENE3,GENE5,GENE6)(GENE4)" << endl;
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
	cout << INPUT_MSA_TAG << "= input.phy" << endl;
	cout << INPUT_TREE_TAG << "= input.tree" << endl;
	cout << INPUT_DATATYPE_TAG << "= nt" << endl;
	cout << INPUT_KEEPBRANCHES_TAG << "= false" << endl << endl;
	cout << "; Start of candidate models description" << endl;
	cout << "[" << MODELS_TAG << "]" << endl;
	cout << MODELS_INCLUDE_TAG << "= all" << endl << endl;
	cout << "; Start of searching options" << endl;
	cout << "[" << SEARCH_TAG << "]" << endl;
	cout << SEARCH_ALGORITHM_TAG << "= auto" << endl;
	cout << SEARCH_ALGORITHM_REPS << "= 1" << endl;
	cout << "[" << PARTITIONS_TAG << "]" << endl;
	cout << "PART1 = INI1-END1" << endl;
	cout << "PART2 = INI2-END2" << endl;
	cout << "..." << endl;
	cout << "PARTn = INIn-ENDn" << endl;
}

vector<partitionInfo> * ConfigParser::getPartitions() {
	return partitions;
}

struct partitionInfo ConfigParser::getPartition(size_t index) {
	if (!partitions) {
		cerr << "[ERROR] No partitions were defined" << endl;
		exit_partest(EX_SOFTWARE);
	}
	if (index >= number_of_genes) {
		cerr << "[ERROR] Requested partition does not exist" << endl;
		exit_partest(EX_SOFTWARE);
	}
	return partitions->at(index);
}

} /* namespace partest */
