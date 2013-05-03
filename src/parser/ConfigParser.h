/*
 * ConfigParser.h
 *
 *  Created on: May 1, 2013
 *      Author: diego
 */

#ifndef CONFIGPARSER_H_
#define CONFIGPARSER_H_

#include "util/GlobalDefs.h"
#include "indata/PartitionElement.h"
#include "SimpleIni.h"
#include <vector>

#define PARTITIONS_TAG 		"partitions"

#define SEARCH_TAG 				"search"
#define SEARCH_ALGORITHM_TAG 	"algorithm"

#define OUTPUT_TAG 			"output"
#define OUTPUT_BASE_PATH 	"path"
#define OUTPUT_MODELS_TAG 	"models"
#define OUTPUT_PARTS_TAG 	"partitions"
#define OUTPUT_SCHEMES_TAG 	"schemes"
#define DEFAULT_OUTPUT_BASE_PATH 	""
#define DEFAULT_OUTPUT_MODELS_TAG 	"models.out"
#define DEFAULT_OUTPUT_PARTS_TAG 	"partitions.out"
#define DEFAULT_OUTPUT_SCHEMES_TAG 	"schemes.out"

using namespace std;
namespace partest {

struct partitionInfo {
	t_partitionElementId partitionId;
	int start;
	int end;
	int stride;
	string name;
	~partitionInfo(void) {
	}
};

class ConfigParser {
public:
	ConfigParser(const char * configFile);
	virtual ~ConfigParser();
	vector<partitionInfo> * getPartitions();
	struct partitionInfo getPartition(int index);
	int getNumberOfPartitions() {
		return numberOfPartitions;
	}
	static void printFormat();
private:
	int parsePartitionLine(char * line, struct partitionInfo * pInfo);
	int parsePartitionDetails(char * line, struct partitionInfo * pInfo);
	const char * configFile;
	int numberOfPartitions;
	vector<partitionInfo> * partitions;
	/** search algorithm **/
	bool searchData;
	SearchAlgo searchAlgorithm;
	string outputFileModels;
	string outputFilePartitions;
	string outputFileSchemes;
	string outputBasePath;
};

} /* namespace partest */
#endif /* CONFIGPARSER_H_ */
