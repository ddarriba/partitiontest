/*
 * PartitionSelector.h
 *
 *  Created on: Mar 25, 2013
 *      Author: diego
 */

#ifndef PARTITIONSELECTOR_H_
#define PARTITIONSELECTOR_H_

#include "indata/Partition.h"
#include "util/GlobalDefs.h"

namespace partest {

typedef struct {
	Partition * partition;
	double lnL;
	double numParameters;
	double value;
	double delta;
	double weight;
} SelectionPartition;

class PartitionSelector {
public:
	PartitionSelector(Partition ** partition,
			int numberOfPartitions, InformationCriterion ic,
			SampleSize sampleSize, double sampleSizeValue = 0.0);
	virtual ~PartitionSelector();
	Partition * getBestPartition(void) { return bestPartition; }
	SelectionPartition * getBestSelectionPartition(void) { return partitionsVector->at(0); }
	void print(void);
private:
	Partition ** partitions;
	Partition * bestPartition;
	int numberOfPartitions;
	InformationCriterion ic;
	SampleSize sampleSize;
	double sampleSizeValue;
	vector<SelectionPartition *> * partitionsVector;
};

} /* namespace partest */
#endif /* PARTITIONSELECTOR_H_ */
