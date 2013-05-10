/*
 * HierarchicalSearchAlgorithm.h
 *
 *  Created on: Apr 1, 2013
 *      Author: diego
 */

#ifndef HIERARCHICALSEARCHALGORITHM_H_
#define HIERARCHICALSEARCHALGORITHM_H_

#include "SearchAlgorithm.h"
#include "indata/PartitioningScheme.h"

namespace partest {

class HierarchicalSearchAlgorithm: public SearchAlgorithm {
public:
	HierarchicalSearchAlgorithm(ParTestOptions * options, PartitionMap * partitionMap);
	virtual ~HierarchicalSearchAlgorithm();
	virtual PartitioningScheme * start();
	virtual void update(const ObservableInfo & info, ParTestOptions * run_instance =
			NULL);
private:
	int numberOfBits;
};

} /* namespace partest */
#endif /* HierarchicalSearchAlgorithm_H_ */
