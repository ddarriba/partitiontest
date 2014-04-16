/*
 * SearchAlgorithm.h
 *
 *  Created on: Apr 10, 2014
 *      Author: diego
 */

#ifndef SEARCHALGORITHM_H_
#define SEARCHALGORITHM_H_

#include "indata/PartitioningScheme.h"
#include "exe/ModelOptimize.h"

namespace partest {

class SearchAlgorithm {
public:
	SearchAlgorithm();
	virtual ~SearchAlgorithm();
	virtual PartitioningScheme * start() = 0;
protected:
	ModelOptimize mo;
};

} /* namespace partest */

#endif /* SEARCHALGORITHM_H_ */
