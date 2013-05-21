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
 * @file ConfigParser.h
 *
 * @brief Utility for parsing configuration files.
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

/**
 * @brief Structure with information about a single-gene partition
 */
struct partitionInfo {
	t_partitionElementId partitionId; /** Partition id */
	int start; /** Starting position */
	int end; /** Ending position */
	int stride; /** Stride for codon position (0 means no codon division) */
	string name; /** Name of the gene/partition */
	~partitionInfo(void) {
	}
};

/**
 * @brief Utility class for parsing configuration files.
 */
class ConfigParser {
public:
	/**
	 * @brief Creates a new ConfigParser
	 *
	 * @param configFile Configuration file name.
	 */
	ConfigParser(const char * configFile);
	virtual ~ConfigParser();

	/**
	 * @brief Gets the whole set of partitions.
	 *
	 * @return The set of partitions.
	 */
	vector<partitionInfo> * getPartitions();

	/**
	 * @brief Gets a single partition.
	 *
	 * @param index The index of the partition to be retrieved
	 *
	 * @return The partition at index
	 */
	struct partitionInfo getPartition(int index);

	/**
	 * @brief Gets the number of partitions.
	 *
	 * @return The number of partitions.
	 */
	int getNumberOfPartitions() {
		return numberOfPartitions;
	}

	/**
	 * @brief Prints format instructions for configuration files.
	 */
	static void printFormat();

	/**
	 * @brief Creates a template for configuration files.
	 */
	static void createTemplate();

	/**
	 * @brief Gets the file name for model selections output.
	 */
	const string& getOutputFileModels() const {
		return outputFileModels;
	}

	/**
	 * @brief Gets the file name for partition selections output.
	 */
	const string& getOutputFilePartitions() const {
		return outputFilePartitions;
	}

	/**
	 * @brief Gets the file name for scheme selections output.
	 */
	const string& getOutputFileSchemes() const {
		return outputFileSchemes;
	}

private:

	/**
	 * @brief Parses a single line from the configuration file.
	 *
	 * @param line The line to be parsed.
	 * @param[out] The partitionInfo structure.
	 *
	 * @return 0 if there was no error.
	 */
	int parsePartitionLine(char * line, struct partitionInfo * pInfo);

	/**
	 * @brief Parses the complete details of a partition.
	 *
	 * @param line The line to be parsed.
	 * @param[out] The partitionInfo structure.
	 *
	 * @return 0 if there was no error.
	 */
	int parsePartitionDetails(char * line, struct partitionInfo * pInfo);

	const char * configFile; /** Configuration file name */
	int numberOfPartitions; /** The numnber of partitions */
	vector<partitionInfo> * partitions; /** Vector of partitions */
	/** search algorithm **/
	bool searchData; /** TODO: Not sure what is this doing here */
	SearchAlgo searchAlgorithm; /** Search algorithm */
	string outputFileModels; /** File name for model selections output */
	string outputFilePartitions; /** File name for partition selections output */
	string outputFileSchemes; /** File name for scheme selections output */
	string outputBasePath; /** Base path for output files */
};

} /* namespace partest */
#endif /* CONFIGPARSER_H_ */
