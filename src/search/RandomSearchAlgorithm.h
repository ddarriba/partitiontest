/*
 * RandomSearchAlgorithm.h
 *
 *  Created on: Apr 1, 2013
 *      Author: diego
 */

#ifndef RANDOMSEARCHALGORITHM_H_
#define RANDOMSEARCHALGORITHM_H_

#include "SearchAlgorithm.h"
#include "indata/Partition.h"

namespace partest {

class RandomSearchAlgorithm: public SearchAlgorithm {
public:
	RandomSearchAlgorithm(ParTestOptions * options, PartitionMap * partitionMap);
	virtual ~RandomSearchAlgorithm();
	virtual Partition * start();
	virtual void update(const ObservableInfo & info, ParTestOptions * run_instance =
			NULL);
private:
	Partition * getRandomPartition(Partition ** partitionsArray, int numberOfPartitions);
	Partition * getRandomPartition(Partition * p0, Partition ** partitionsArray, int numberOfPartitions);
	int numberOfBits;
};

} /* namespace partest */
#endif /* RANDOMSEARCHALGORITHM_H_ */
