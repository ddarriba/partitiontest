/*
 * Partition.h
 *
 *  Created on: Mar 8, 2013
 *      Author: diego
 */

#ifndef PARTITION_H_
#define PARTITION_H_

#include "PartitionElement.h"
#include "PartitionMap.h"
#include "util/GlobalDefs.h"
#include <string>

namespace partest {

class Partition {
public:
	Partition(int numberOfElements);
	Partition(t_partition_elements * partition, PartitionMap * partitionMap);
	int getNumberOfElements() { return numberOfElements; }
	int getNumberOfBits() { return numberOfBits; }
	int addElement(PartitionElement * element);
	PartitionElement * getElement(int id);
	string toString();
	virtual ~Partition();
private:
	PartitionElement ** partitions;
	int currentElement;
	int numberOfElements;
	int numberOfBits;
	string * code;
};

} /* namespace partest */
#endif /* PARTITION_H_ */
