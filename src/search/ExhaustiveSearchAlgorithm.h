/*
 * ExhaustiveSearchAlgorithm.h
 *
 *  Created on: Mar 8, 2013
 *      Author: diego
 */

#ifndef EXHAUSTIVESEARCHALGORITHM_H_
#define EXHAUSTIVESEARCHALGORITHM_H_

#include "SearchAlgorithm.h"

namespace partest {

class ExhaustiveSearchAlgorithm: public SearchAlgorithm {
public:
	ExhaustiveSearchAlgorithm(ParTestOptions * options, PartitionMap * partitionMap);
	virtual ~ExhaustiveSearchAlgorithm();
	virtual PartitioningScheme * start();
	virtual void update(const ObservableInfo & info,
	        ParTestOptions * run_instance = NULL);
};

} /* namespace partest */
#endif /* EXHAUSTIVESEARCHALGORITHM_H_ */
