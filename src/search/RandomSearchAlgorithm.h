/*
 * RandomSearchAlgorithm.h
 *
 *  Created on: Apr 1, 2013
 *      Author: diego
 */

#ifndef RANDOMSEARCHALGORITHM_H_
#define RANDOMSEARCHALGORITHM_H_

#include "SearchAlgorithm.h"
#include "indata/PartitioningScheme.h"

namespace partest {

class RandomSearchAlgorithm: public SearchAlgorithm {
public:
	RandomSearchAlgorithm(ParTestOptions * options, PartitionMap * partitionMap);
	virtual ~RandomSearchAlgorithm();
	virtual PartitioningScheme * start();
	virtual void update(const ObservableInfo & info, ParTestOptions * run_instance =
			NULL);
private:
	PartitioningScheme * getRandomPartitioningScheme(PartitioningScheme ** schemesArray, int numberOfSchemes);
	PartitioningScheme * getRandomPartitioningScheme(PartitioningScheme * p0, PartitioningScheme ** schemesArray, int numberOfSchemes);
	int numberOfBits;
};

} /* namespace partest */
#endif /* RANDOMSEARCHALGORITHM_H_ */
