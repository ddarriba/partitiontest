/*
 * HierarchicalClusteringSearchAlgorithm.h
 *
 *  Created on: Apr 9, 2014
 *      Author: diego
 */

#ifndef HIERARCHICALCLUSTERING_H_
#define HIERARCHICALCLUSTERING_H_

#include "search/SearchAlgorithm.h"

namespace partest {

class HierarchicalClusteringSearchAlgorithm : public SearchAlgorithm {
public:
	HierarchicalClusteringSearchAlgorithm();
	virtual ~HierarchicalClusteringSearchAlgorithm();
	virtual PartitioningScheme * start(PartitioningScheme * startingPoint = 0);
};

} /* namespace partest */

#endif /* HIERARCHICALCLUSTERING_H_ */
