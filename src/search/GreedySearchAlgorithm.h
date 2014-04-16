/*
 * GreedySearchAlgorithm.h
 *
 *  Created on: Apr 11, 2014
 *      Author: diego
 */

#ifndef GREEDYSEARCHALGORITHM_H_
#define GREEDYSEARCHALGORITHM_H_

#include "SearchAlgorithm.h"

#include "util/GlobalDefs.h"
#include <vector>

namespace partest {

class GreedySearchAlgorithm: public partest::SearchAlgorithm {
public:
	GreedySearchAlgorithm();
	virtual ~GreedySearchAlgorithm();
	virtual PartitioningScheme * start();
private:
	vector<PartitioningScheme *> getNextSchemes(const t_partitioningScheme * startingScheme);
};

} /* namespace partest */

#endif /* GREEDYSEARCHALGORITHM_H_ */
