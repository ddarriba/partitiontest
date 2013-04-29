/*
 * GreedySearchAlgorithm.h
 *
 *  Created on: Apr 1, 2013
 *      Author: diego
 */

#ifndef GREEDYSEARCHALGORITHM_H_
#define GREEDYSEARCHALGORITHM_H_

#include "SearchAlgorithm.h"
#include "indata/Partition.h"

namespace partest {

class GreedySearchAlgorithm: public SearchAlgorithm {
public:
	GreedySearchAlgorithm(ParTestOptions * options, PartitionMap * partitionMap);
	virtual ~GreedySearchAlgorithm();
	virtual PartitioningScheme * start();
	virtual void update(const ObservableInfo & info, ParTestOptions * run_instance =
			NULL);
private:
	int getNextPartitioningSchemes(PartitioningScheme *currentScheme, PartitioningScheme **nextSchemes);
	int numberOfBits;
};

} /* namespace partest */
#endif /* GreedySearchAlgorithm_H_ */
