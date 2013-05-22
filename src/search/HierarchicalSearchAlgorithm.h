/*  PartitionTest, fast selection of the best fit partitioning scheme for
 *  multi-gene data sets.
 *  Copyright May 2013 by Diego Darriba
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 3 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 *
 *  For any other inquiries send an Email to Diego Darriba
 *  ddarriba@udc.es
 */

/**
 * @file HierarchicalSearchAlgorithm.h
 */

#ifndef HIERARCHICALSEARCHALGORITHM_H_
#define HIERARCHICALSEARCHALGORITHM_H_

#include "SearchAlgorithm.h"
#include "indata/PartitioningScheme.h"

namespace partest {

/**
 * @brief Step-by-step hierarchical clustering algorithm.
 *
 * This algorithm starts evaluating the K=N partitioning scheme and merges the
 * two more similar partitions. This process is repeated until a maximum is reached
 * or there is a single partition.
 */
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
