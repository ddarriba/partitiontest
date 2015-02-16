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
 * @file PartestParserUtils.h
 * @author Diego Darriba
 * @brief Set of utilities for the PartitionTest parser tool
 */

#ifndef PARTESTPARSERUTILS_H_
#define PARTESTPARSERUTILS_H_

#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#define MAX_SECTIONS 30

namespace partest_parser {

enum selection {
	icAIC, icAICc, icBIC, icUNDEFINED
};

enum pfSearch {
	searchAll, searchHCluster, searchRCluster, searchGreedy, searchUNDEFINED
};

/**
 * @brief Structure with information about a single-gene partition
 */
struct partitionInfo {
	int start[MAX_SECTIONS]; /** Starting position */
	int end[MAX_SECTIONS]; /** Ending position */
	int stride[MAX_SECTIONS]; /** Stride for codon position (0 means no codon division) */
	int numberOfSections;
	std::string name; /** Name of the gene/partition */
	~partitionInfo(void) {
	}
};

struct executionInfo {
	std::string alignment;
	std::string models;
	std::string modelSelection;
	pfSearch searchAlgorithm;
};

class PartestParserUtils {
public:
	PartestParserUtils(char * inputFile, char * outputFile);
	virtual ~PartestParserUtils();
	int parseRaxmlFile(std::vector<std::string> ** partitions);
	int parsePartitionFinderFile(std::vector<std::string> ** partitions,
			char ** alignment, char ** models, pfSearch * searchAlgo,
			selection * icSelection);
	static inline int isNum(char c) {
		return !(c < '0' || c > '9');
	}
	static inline bool existsFile(const std::string& name) {
		if (FILE *file = fopen(name.c_str(), "r")) {
			fclose(file);
			return true;
		} else {
			return false;
		}
	}
private:
	static void analyzeIdentifier(char **ch, char **modelName);
	static int myGetline(char **lineptr, int *n, FILE *stream);
	static void skipWhites(char **ch);
	static int lineContainsOnlyWhiteChars(char * str);

	char * inputFile;
	char * outputFile;
};

} /* namespace partest */

#endif /* PARTESTPARSERUTILS_H_ */
