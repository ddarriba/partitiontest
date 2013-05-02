/*
 * ConfigParser.cpp
 *
 *  Created on: May 1, 2013
 *      Author: diego
 */

#include "ConfigParser.h"
#include "util/Utilities.h"
#include "SimpleIni.h"
#include <string.h>
#include <assert.h>

namespace partest {

ConfigParser::ConfigParser(const char * configFile) :
		configFile(configFile), numberOfPartitions(0) {

	int partitionId = 0;

	CSimpleIniA ini;
	ini.SetUnicode();
	SI_Error rc = ini.LoadFile(configFile);
	if (rc < 0)
		exit(-1);

	CSimpleIniA::TNamesDepend keys;

	/** SEARCH ALGORITHM **/

	ini.GetAllKeys(SEARCH_ALGO_TAG, keys);
	searchData = keys.size() > 0;
	if (searchData) {
		const char * value = ini.GetValue(SEARCH_ALGO_TAG, ALGORITHM_TAG,
				"none");
		if (!strcmp(value, "greedy")) {
			searchAlgorithm = SearchGreedy;
		} else if (!strcmp(value, "exhaustive")) {
			searchAlgorithm = SearchExhaustive;
		} else if (!strcmp(value, "random")) {
			searchAlgorithm = SearchRandom;
		} else {
			cerr << "Invalid search algorithm : " << value << endl;
			exit(-1);
		}
	}

	/** PARTITIONS **/

	ini.GetAllKeys(PARTITIONS_TAG, keys);
	numberOfPartitions = keys.size();

	assert(
			Utilities::binaryPow(numberOfPartitions) < sizeof(t_partitionElementId) * 8);
	partitions = new vector<partitionInfo>(numberOfPartitions);

	char * lineBuffer = (char *) malloc(30);
	for (CSimpleIniA::TNamesDepend::iterator it = keys.begin();
			it != keys.end(); it++) {
		CSimpleIniA::Entry entry = *it;
		partitions->at(partitionId).partitionId = Utilities::binaryPow(
				partitionId);
		strcpy(lineBuffer,
				ini.GetValue(PARTITIONS_TAG, entry.pItem, "default"));
		parsePartitionDetails(lineBuffer, &partitions->at(partitionId));
		partitionId++;
	}
	free(lineBuffer);
	exit(0);
}

int ConfigParser::parsePartitionDetails(char * line,
		struct partitionInfo * pInfo) {
	int start = atoi(strtok(line, "-"));
	int end = atoi(strtok(NULL, "\\"));
	char * strideStr = strtok(NULL, "\\");
	int stride = strideStr ? atoi(strideStr) : 0;
	pInfo->start = start;
	pInfo->end = end;
	pInfo->stride = stride;
	return 0;
}

int ConfigParser::parsePartitionLine(char * line,
		struct partitionInfo * pInfo) {
	char * name = strtok(line, "=");
	string nameStr(name);
	int start = atoi(strtok(NULL, "-"));
	int end = atoi(strtok(NULL, "\\"));
	char * strideStr = strtok(NULL, "\\");
	int stride = strideStr ? atoi(strideStr) : 0;
	pInfo->start = start;
	pInfo->end = end;
	pInfo->stride = stride;
	pInfo->name = nameStr;
//	if (line)
//		free (line);
//	line = (char *) NULL;
	return 0;
}

ConfigParser::~ConfigParser() {
	delete partitions;
}

struct partitionInfo ConfigParser::getPartition(int index) {
	return partitions->at(index);
}

void ConfigParser::printFormat() {
	cout << "Config file format:" << endl << endl;
	cout << "   [partitions]" << endl;
	cout << "   ;THIS IS A COMMENT" << endl;
	cout << "   PART1=INI1-END1[\\STRIDE1]" << endl;
	cout << "   PART1=INI2-END2[\\STRIDE2]" << endl;
	cout << "   ..." << endl;
	cout << "   PART1=INI3-END3[\\STRIDE3]" << endl;
	cout << endl << "Example:" << endl << endl;
	cout << "   ; Start of partitions for file.phy" << endl;
	cout << "   [partitions]" << endl;
	cout << "   DNA1=1-976" << endl;
	cout << "   DNA2COD1=976-1803\\1" << endl;
	cout << "   DNA2COD2=976-1803\\2" << endl;
	cout << "   DNA2COD3=976-1803\\3" << endl;
	cout << "   ; End of partitions for file.phy" << endl;
}

vector<partitionInfo> * ConfigParser::getPartitions() {
	return partitions;
}

} /* namespace partest */
