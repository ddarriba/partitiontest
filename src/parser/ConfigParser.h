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
#include <vector>

using namespace std;
namespace partest {

struct partitionInfo {
	t_partitionElementId partitionId;
	int start;
	int end;
	int stride;
	string name;
	~partitionInfo(void){ }
};

class ConfigParser {
public:
	ConfigParser(const char * configFile);
	virtual ~ConfigParser();
	vector<partitionInfo> * getPartitions();
	struct partitionInfo getPartition(int index);
	int getNumberOfPartitions() { return numberOfPartitions; }
	static void printFormat();
private:
	int parsePartitionLine(char * line, struct partitionInfo * pInfo);
	const char * configFile;
	int numberOfPartitions;
	vector<partitionInfo> * partitions;
};

} /* namespace partest */
#endif /* CONFIGPARSER_H_ */
