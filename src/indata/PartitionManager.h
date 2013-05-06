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
 *  For any other enquiries send an Email to Diego Darriba
 *  ddarriba@udc.es
 */

/**
 * @file PartitionManager.h
 *
 * @brief Iterative manager for generating the workload in each stage.
 */

#ifndef PARTITIONMANAGER_H_
#define PARTITIONMANAGER_H_

#include "util/GlobalDefs.h"

namespace partest {

/**
 * @brief Iterative manager for generating the workload in each stage.
 */
class PartitionManager {
  public:
    /**
     * Creates a new manager from current execution options.
     */
    PartitionManager();
    virtual ~PartitionManager();

    /**
     * Get all possible permutations from a set of gene groups. For example:
     * \li 1000, 0100, 0010, 0001, 0110, 1100, 0101
     * Would generate:
     * \li 1000, 0100, 0010, 0001
     * \li 1000, 0001, 0110
     * \li 0010, 0001, 1100
     * \li 1000, 0010, 0101
     */
	static void getPermutations(t_partitionElementId * mask, unsigned int size, t_partitionElementId full_mask,
    		t_schemesVector * schemeVector, t_partitionElementId limit = MAX_PARTITIONS, t_partitionElementId current = -1,
        int sum = 0);
};

} /* namespace partest */
#endif /* PARTITIONMANAGER_H_ */
