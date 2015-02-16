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
 * @author Diego Darriba
 * @brief Algorithm for performing a greedy search
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
	virtual PartitioningScheme * start(PartitioningScheme * startingPoint = 0);
private:
	std::vector<PartitioningScheme *> getNextSchemes(const t_partitioningScheme * startingScheme);
};

} /* namespace partest */

#endif /* GREEDYSEARCHALGORITHM_H_ */
