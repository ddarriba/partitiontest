/*
 * SearchAlgorithm.h
 *
 *  Created on: Apr 10, 2014
 *      Author: diego
 */

#ifndef SEARCHALGORITHM_H_
#define SEARCHALGORITHM_H_

#include "indata/PartitioningScheme.h"

namespace partest {

class SearchAlgorithm {
public:
	virtual PartitioningScheme * start() = 0;
};

} /* namespace partest */

#endif /* SEARCHALGORITHM_H_ */
