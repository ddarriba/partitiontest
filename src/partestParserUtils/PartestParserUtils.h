/*
 * PartestParserUtils.h
 *
 *  Created on: Jul 1, 2014
 *      Author: diego
 */

#ifndef PARTESTPARSERUTILS_H_
#define PARTESTPARSERUTILS_H_

#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#define MAX_SECTIONS 30

namespace partest_parser {

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

class PartestParserUtils {
public:
	PartestParserUtils(char * inputFile, char * outputFile);
	virtual ~PartestParserUtils();
	int parseRaxmlFile(std::vector<partitionInfo> ** partitions);
	int parsePartitionFinderFile(std::vector<std::string> ** partitions);
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
	static void analyzeIdentifier(char **ch, int modelNumber, char **modelName);
	static int myGetline(char **lineptr, int *n, FILE *stream);
	static void skipWhites(char **ch);
	static int lineContainsOnlyWhiteChars(char * str);

	char * inputFile;
	char * outputFile;
};

} /* namespace partest */

#endif /* PARTESTPARSERUTILS_H_ */
