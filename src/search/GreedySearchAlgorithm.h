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
 * @file GreedySearchAlgorithm.h
 */

#ifndef GREEDYSEARCHALGORITHM_H_
#define GREEDYSEARCHALGORITHM_H_

#include "SearchAlgorithm.h"
#include "indata/PartitioningScheme.h"

namespace partest {

/**
 * @brief Multi-stage greedy search algorithm
 *
 * This algorithm starts with each gene-partition separately and analyzes all
 * pairwise possibilities, merging those which make the best-fit partition.
 */
class GreedySearchAlgorithm: public SearchAlgorithm {
public:
	GreedySearchAlgorithm(ParTestOptions * options, PartitionMap * partitionMap);
	virtual ~GreedySearchAlgorithm();
	virtual PartitioningScheme * start();
	virtual PartitioningScheme * start(PartitioningScheme * startingPoint);
	virtual void update(const ObservableInfo & info, ParTestOptions * run_instance =
			NULL);
private:
	int getNextPartitioningSchemes(PartitioningScheme *currentScheme, PartitioningScheme **nextSchemes);
	int numberOfBits;
};

} /* namespace partest */
#endif /* GreedySearchAlgorithm_H_ */
