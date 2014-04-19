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
 * @file ExhaustiveSearchAlgorithm.h
 */

#ifndef EXHAUSTIVESEARCHALGORITHM_H_
#define EXHAUSTIVESEARCHALGORITHM_H_

#include "SearchAlgorithm.h"

namespace partest {

/**
 * @brief Single-step exhaustive search algorithm
 *
 * This algorithm constructs and evaluates all possible partitioning schemes.
 * It is not recommended at all unless the number of gene-partitions is very very low.
 */
class ExhaustiveSearchAlgorithm: public SearchAlgorithm {
public:
	ExhaustiveSearchAlgorithm();
	virtual ~ExhaustiveSearchAlgorithm();
	virtual PartitioningScheme * start();
	virtual PartitioningScheme * start(PartitioningScheme * startingPoint);
};

} /* namespace partest */
#endif /* EXHAUSTIVESEARCHALGORITHM_H_ */
