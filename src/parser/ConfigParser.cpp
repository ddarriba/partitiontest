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

namespace partest {

ConfigParser::ConfigParser(const char * configFile) :
		configFile(configFile), numberOfPartitions(0) {
	FILE *f;
	char *cc = (char *) NULL;
	int nbytes;
	int partitionId = 0;

	//int **partitions;

	f = Utilities::myfopen(configFile, "rb", true);

	while (Utilities::myGetline(&cc, &nbytes, f) > -1) {

		numberOfPartitions++;

		if (cc)
			free(cc);
		cc = (char *) NULL;
	}

	rewind(f);

	assert(
			Utilities::binaryPow(numberOfPartitions) < sizeof(t_partitionElementId) * 8);
	partitions = new vector<partitionInfo>(numberOfPartitions);

	while (Utilities::myGetline(&cc, &nbytes, f) > -1) {

		char * name = strtok(cc, "=");
		string nameStr(name);
		int start = atoi(strtok(NULL, "-"));
		int end = atoi(strtok(NULL, "\\"));
		char * strideStr = strtok(NULL, "\\");
		int stride = strideStr ? atoi(strideStr) : 0;

		/* partitionId is translated into partition mask */
		partitions->at(partitionId).partitionId = Utilities::binaryPow(
				partitionId);
		partitions->at(partitionId).start = start;
		partitions->at(partitionId).end = end;
		partitions->at(partitionId).stride = stride;
		partitions->at(partitionId).name = nameStr;

		partitionId++;

		if (cc)
			free(cc);

		cc = (char *) NULL;
	}
}

ConfigParser::~ConfigParser() {
	delete partitions;
}

struct partitionInfo ConfigParser::getPartition(int index) {
	return partitions->at(index);
}

void ConfigParser::printFormat() {
	cout << "Config file format:" << endl << endl;
	cout << "   #THIS IS A COMMENT" << endl;
	cout << "   PART1=INI1-END1[\\STRIDE1]" << endl;
	cout << "   PART1=INI2-END2[\\STRIDE2]" << endl;
	cout << "   ..." << endl;
	cout << "   PART1=INI3-END3[\\STRIDE3]" << endl;
	cout << endl << "Example:" << endl << endl;
	cout << "   # Start of partitions for file.phy" << endl;
	cout << "   DNA1=1-976" << endl;
	cout << "   DNA2COD1=976-1803\\1" << endl;
	cout << "   DNA2COD2=976-1803\\2" << endl;
	cout << "   DNA2COD3=976-1803\\3" << endl;
	cout << "   # End of partitions for file.phy" << endl;
}

vector<partitionInfo> * ConfigParser::getPartitions() {
	return partitions;
}

} /* namespace partest */
