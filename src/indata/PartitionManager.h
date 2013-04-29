/*
 * PartitionManager.h
 *
 *  Created on: 19/06/2012
 *      Author: diego
 */

#ifndef PARTITIONMANAGER_H_
#define PARTITIONMANAGER_H_

#include "../util/GlobalDefs.h"

namespace partest {

/*!
 * \brief Iterative manager for generating the work in each stage.
 */
class PartitionManager {
  public:
    /*!
     * Creates a new manager from current execution options.
     */
    PartitionManager();
    virtual ~PartitionManager();

    /*!
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
