/*
 * PartitionTestParser.h
 *
 *  Created on: Jul 1, 2014
 *      Author: diego
 */

#ifndef PARTITIONTESTPARSER_H_
#define PARTITIONTESTPARSER_H_

#include "partestParserUtils/PartestParserUtils.h"

#include <string>
#include <vector>
#include <iostream>

namespace partest_parser {

enum ConfigFormat {
	cfPartitionFinder, /** PartitionFinder format */
	cfRAxML /** RAxML format */
};

class PartitionTestParser {
public:
	PartitionTestParser(int argc, char *argv[]);
	virtual ~PartitionTestParser();

	int parseConfigFile();
private:
	void printHelp(std::ostream & out);
	int getFormat(char * str, ConfigFormat * result);

	std::vector<std::string> * partitionStrings; /** Vector of partitions */
	char inputFile[256];
	char outputFile[256];
	char binName[50];
	ConfigFormat format;
};

} /* namespace partest */

#endif /* PARTITIONTESTPARSER_H_ */
