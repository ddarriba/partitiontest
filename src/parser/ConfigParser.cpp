/*
 * ConfigParser.cpp
 *
 *  Created on: May 1, 2013
 *      Author: diego
 */

#include "ConfigParser.h"
#include "util/Utilities.h"
#include <string.h>
#include <assert.h>
#include <algorithm>

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

ConfigParser::ConfigParser(const char * configFile) :
		configFile(configFile), partitions(0), numberOfPartitions(0), outputBasePath(
				DEFAULT_OUTPUT_BASE_PATH), outputFileModels(
				DEFAULT_OUTPUT_MODELS_TAG), outputFilePartitions(
				DEFAULT_OUTPUT_PARTS_TAG), outputFileSchemes(
				DEFAULT_OUTPUT_SCHEMES_TAG), outputFileResults(
				DEFAULT_OUTPUT_RESULTS_TAG) {

	if (configFile != 0 && strcmp(configFile, "")) {
		int partitionId = 0;
		const char * value;
		CSimpleIniA ini;
		ini.SetUnicode();
		SI_Error rc = ini.LoadFile(configFile);
		if (rc < 0)
			Utilities::exit_partest(EX_IOERR);

		CSimpleIniA::TNamesDepend keys;

		/** SEARCH ALGORITHM **/

		ini.GetAllKeys(SEARCH_TAG, keys);
		searchData = keys.size() > 0;
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
				Utilities::exit_partest(EX_SOFTWARE);
			}
		}

		/** PARTITIONS **/
		ini.GetAllKeys(PARTITIONS_TAG, keys);
		numberOfPartitions = keys.size();
//		assert( numberOfPartitions < sizeof(t_partitionElementId) * 8);
		partitions = new vector<partitionInfo>(numberOfPartitions);
		char * lineBuffer = (char *) malloc(150);
		for (CSimpleIniA::TNamesDepend::iterator it = keys.begin();
				it != keys.end(); it++) {
			CSimpleIniA::Entry entry = *it;
			partitions->at(partitionId).name = entry.pItem;
			strcpy(lineBuffer,
					ini.GetValue(PARTITIONS_TAG, entry.pItem, "default"));
			parsePartitionDetails(lineBuffer, &partitions->at(partitionId));
			partitionId++;
		}
		free(lineBuffer);

		std::sort(partitions->begin(), partitions->end(),
				comparePartitionInfos());

#ifdef _PLL
		struct pllPartitionRegion * pregion;
		struct pllPartitionInfo * pinfo;

		pllQueueInit(&parts);
#endif

		for (int i = 0; i < numberOfPartitions; i++) {
			partitions->at(i).partitionId.push_back(i);
#ifdef _PLL
			pinfo = (pllPartitionInfo *) malloc(
					sizeof(struct pllPartitionInfo));
			pllQueueInit(&(pinfo->regionList));
			pllQueueAppend(parts, (void *) pinfo);

			pinfo->partitionName= (char *) malloc ((partitions->at(i).name.size() + 1) * sizeof(char));
			strcpy(pinfo->partitionName, partitions->at(i).name.c_str());
			pinfo->partitionModel = (char *) malloc (1);

			pinfo->protModels = -1;
			pinfo->protFreqs = -1;
			pinfo->dataType = DNA_DATA;
			pinfo->optimizeBaseFrequencies = PLL_TRUE;
			for (int j = 0; j < partitions->at(i).numberOfSections; j++) {
				pregion = (struct pllPartitionRegion *) malloc(
						sizeof(struct pllPartitionRegion));
				pregion->start = partitions->at(i).start[j];
				pregion->end = partitions->at(i).end[j];
				pregion->stride = partitions->at(i).stride[j];
				pllQueueAppend(pinfo->regionList, (void *) pregion);
			}
#endif
		}

		/** OUTPUT **/
		value = ini.GetValue(OUTPUT_TAG, OUTPUT_BASE_PATH, 0);
		if (value) {
			outputBasePath = string(value);
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
#ifdef _PLL
	pllQueuePartitionsDestroy(&parts);
#endif
	if (partitions) {
		delete partitions;
	}
}

struct partitionInfo ConfigParser::getPartition(int index) {
	if (!partitions) {
		cerr << "ERROR: No partitions were defined" << endl;
		Utilities::exit_partest(EX_SOFTWARE);
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
	cout << endl << "   ; Start of searching options" << endl;
	cout << "   [" << SEARCH_TAG << "]" << endl;
	cout << "   " << SEARCH_ALGORITHM_TAG << "={greedy|random|exhaustive}"
			<< endl;
	cout << endl << "   ; Start of partitions for file.phy" << endl;
	cout << "   [" << PARTITIONS_TAG << "]" << endl;
	cout << "   PART1=INI1-END1[\\STRIDE1]" << endl;
	cout << "   PART1=INI2-END2[\\STRIDE2]" << endl;
	cout << "   ..." << endl;
	cout << "   PART1=INI3-END3[\\STRIDE3]" << endl;
	cout << endl << "Example:" << endl << endl;
	cout << "   [" << PARTITIONS_TAG << "]" << endl;
	cout << "   DNA1=1-976" << endl;
	cout << "   DNA2COD1=976-1803\\1" << endl;
	cout << "   DNA2COD2=976-1803\\2" << endl;
	cout << "   DNA2COD3=976-1803\\3" << endl;
	cout << "   ; End of partitions for file.phy" << endl;
}

void ConfigParser::createTemplate() {
	cout << ";THIS IS A COMMENT" << endl;
	cout << "; Start of input data" << endl;
	cout << "[" << INPUT_TAG << "]" << endl;
	cout << INPUT_MSA_TAG << "=INPUT_ALIGNMENT_FILE" << endl;
	cout << INPUT_TREE_TAG << "=INPUT_TREE_FILE" << endl;
	cout << "; Start of searching options" << endl;
	cout << "[" << SEARCH_TAG << "]" << endl;
	cout << SEARCH_ALGORITHM_TAG << "={greedy|random|exhaustive}" << endl;
	cout << "[" << PARTITIONS_TAG << "]" << endl;
	cout << "PART1=INI1-END1[\\STRIDE1]" << endl;
	cout << "PART1=INI2-END2[\\STRIDE2]" << endl;
	cout << "..." << endl;
	cout << "PART1=INI3-END3[\\STRIDE3]" << endl;
}

vector<partitionInfo> * ConfigParser::getPartitions() {
	return partitions;
}

} /* namespace partest */
