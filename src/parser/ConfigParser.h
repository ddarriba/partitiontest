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

#define MAX_SECTIONS 20
#include "util/GlobalDefs.h"
#include <vector>

#define PARTITIONS_TAG 		"partitions"
#define SCHEMES_TAG	        "schemes"

#define SEARCH_TAG 				"search"
#define SEARCH_ALGORITHM_TAG 	"algorithm"
#define SEARCH_ALGORITHM_REPS 	"replicates"

#define MODELS_TAG			"models"
#define MODELS_INCLUDE_TAG	"include"

#define INPUT_TAG 			"input"
#define INPUT_MSA_TAG 		"msa"
#define INPUT_TREE_TAG		"tree"
#define INPUT_DATATYPE_TAG	"datatype"

#define OUTPUT_TAG 			"output"
#define OUTPUT_TMP_PATH 	"tmp"
#define OUTPUT_BASE_PATH 	"path"
#define OUTPUT_RESULTS_TAG 	"results"
#define OUTPUT_MODELS_TAG 	"models"
#define OUTPUT_PARTS_TAG 	"partitions"
#define OUTPUT_SCHEMES_TAG 	"schemes"
#ifdef _WIN32
#define DEFAULT_OUTPUT_TMP_PATH 	""
#else
#define DEFAULT_OUTPUT_TMP_PATH 	"/tmp"
#endif
#define DEFAULT_OUTPUT_BASE_PATH 	""
#define DEFAULT_OUTPUT_RESULTS_TAG 	"results.out"
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
	int start[MAX_SECTIONS]; /** Starting position */
	int end[MAX_SECTIONS]; /** Ending position */
	int stride[MAX_SECTIONS]; /** Stride for codon position (0 means no codon division) */
	int numberOfSections;
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
	struct partitionInfo getPartition(size_t index);

	/**
	 * @brief Gets the number of partitions.
	 *
	 * @return The number of partitions.
	 */
	int getNumberOfPartitions() {
		return number_of_genes;
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
	 * @brief Gets the file name for results output.
	 */
	const string& getOutputFileResults() const {
		return outputFileResults;
	}

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

	/**
	 * @brief Gets the base path for output files.
	 */
	const string& getOutputBasePath() const {
		return outputBasePath;
	}

	/**
	 * @brief Gets the path for temporary files.
	 */
	const string& getOutputTmpPath() const {
		return outputTmpPath;
	}

	/**
	 * @brief Gets the number of replicates for hcluster.
	 */
	const int getMaxSamples() const {
		return maxSamples;
	}

	DataType getDataType(void) {
		return dataType;
	}

	const char * getInputFile(void) {
		return inputFile;
	}

	const char * getUserTree(void) {
		return userTree;
	}

	OptimizeMode getOptimizeMode(void) {
		return optimizeMode;
	}

private:

	/**
	 * @brief Parses a partitioning scheme from a line definition.
	 *
	 * @param line The line to be parsed.
	 * @param[out] scheme The partitioning scheme.
	 *
	 * @return 0 if there was no error.
	 */
	int parseScheme(char * line, t_partitioningScheme * scheme);

	/**
	 * @brief Parses the complete details of a partition.
	 *
	 * @param line The line to be parsed.
	 * @param[out] pInfo The partitionInfo structure.
	 *
	 * @return 0 if there was no error.
	 */
	int parsePartitionDetails(char * line, struct partitionInfo * pInfo);

	void createSinglePartition();

	const char * configFile; /** Configuration file name */
	char inputFile[256]; /** User input alignment */
	char userTree[256]; /** User input tree */
	OptimizeMode optimizeMode; /** Mode of optimization */
	vector<partitionInfo> * partitions; /** Vector of partitions */

	/** search algorithm **/
	bool searchData; /** TODO: Not sure what is this doing here */
	SearchAlgo searchAlgorithm; /** Search algorithm */
	DataType dataType;
	int maxSamples; /** Maximum number of samples in HCluster */
	string outputFileResults; /** File name for results output */
	string outputFileModels; /** File name for model selections output */
	string outputFilePartitions; /** File name for partition selections output */
	string outputFileSchemes; /** File name for scheme selections output */
	string outputBasePath; /** Base path for output files */
	string outputTmpPath; /** Base path for output files */
};

} /* namespace partest */
#endif /* CONFIGPARSER_H_ */
