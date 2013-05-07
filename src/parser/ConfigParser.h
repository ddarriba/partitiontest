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
 *  For any other enquiries send an Email to Diego Darriba
 *  ddarriba@udc.es
 */

/**
 * @file ConfigParser.h
 *
 * @brief Util for parsing configuration files.
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

#define INPUT_TAG 			"input"
#define INPUT_MSA_TAG 		"msa"
#define INPUT_TREE_TAG		"tree"
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
	static void createTemplate();

	const string& getOutputFileModels() const {
		return outputFileModels;
	}

	const string& getOutputFilePartitions() const {
		return outputFilePartitions;
	}

	const string& getOutputFileSchemes() const {
		return outputFileSchemes;
	}

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
