/*
 * PartitionMap.h
 *
 *  Created on: Feb 13, 2013
 *      Author: diego
 */

#ifndef PARTITIONMAP_H_
#define PARTITIONMAP_H_

#include <stdlib.h>
#include <vector>
#include "PartitionElement.h"
#include "../util/GlobalDefs.h"

namespace partest {

struct partitionInfo {
	t_partitionElementId partitionId;
	PartitionElement * partitionElement;
	~partitionInfo(void){ }//delete partitionElement; }
};

/**
 * This class represent the mapping of the simplest and unbreakable partitions.
 * Each partition has a power of two Id, what allows to easy identify a complex partition from the id.
 * e.g., Partition 13 (1101) will be a merge of partitions 8 (1000),4 (0100) and 1 (0001).
 */
class PartitionMap {
public:
	PartitionMap(Alignment * alignment, unsigned int numberOfPartitions, bitMask rateVariation, DataType dataType);
	PartitionMap(const char * configFile, Alignment * alignment, bitMask rateVariation, DataType dataType);
	virtual ~PartitionMap();
	bool addPartitionElement(unsigned int partitionId, string name, unsigned int iniPosition,
			unsigned int endPosition, char stride);
	PartitionElement * getPartitionElement(t_partitionElementId partitionId);
	unsigned int getNumberOfElements() { return numberOfElements; }
	unsigned int getNumberOfPartitions() { return numberOfPartitions; }
private:
	Alignment * alignment;
	unsigned int numberOfElements;
	unsigned int numberOfPartitions;
	vector<partitionInfo> * partitions;
	bitMask rateVariation;
	DataType dataType;
};

} /* namespace partest */
#endif /* PARTITIONMAP_H_ */
