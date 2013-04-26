/*
 * PartitionManager.cc
 *
 *  Created on: 19/06/2012
 *      Author: diego
 */

#include "PartitionManager.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

namespace partest {

PartitionManager::PartitionManager() {
}

void vprint(const t_partitionElementId list[], t_partitionElementId size) {
	std::cout << "{ ";
	for (t_partitionElementId i = 0; i < size; i++)
		std::cout << list[i] << ", ";
	std::cout << "}" << std::endl;
}

void PartitionManager::get_permutations(t_partitionElementId * mask,
		unsigned int size, t_partitionElementId full_mask,
		t_partition * partitions, t_partitionElementId limit,
		t_partitionElementId current, int sum) {

	if (!sum && current >= limit)
		return;

	if (current >= 0) {
		sum |= mask[current];
		mask[current] = -mask[current];
	}

	if (sum == full_mask) {
		t_partition_elements next_partition;
		for (t_partitionElementId i = 0; i < size; i++) {
			if (mask[i] < 0) {
				next_partition.push_back(-mask[i]);
			}
		}
		partitions->push_back(next_partition);
		return;
	}

	for (unsigned int i = current + 1; i < size; i++) {
		if (sum & mask[i])
			continue;
		get_permutations(mask, size, full_mask, partitions, limit, i, sum);
		mask[i] = -mask[i];
	}

}

PartitionManager::~PartitionManager() {
}

} /* namespace partest */
