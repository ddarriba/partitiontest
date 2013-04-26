/*
 * SearchAlgorithm.cpp
 *
 *  Created on: Mar 8, 2013
 *      Author: diego
 */

#include "SearchAlgorithm.h"

namespace partest {

SearchAlgorithm::SearchAlgorithm(ParTestOptions * options, PartitionMap * partitionMap) :
		options(options), partitionMap(partitionMap) {
}

SearchAlgorithm::~SearchAlgorithm() {
}


PartitionElement * SearchAlgorithm::getPartitionElement(t_partitionElementId id) {
	return partitionMap->getPartitionElement(id);
}

unsigned int SearchAlgorithm::getNumberOfElements() {
	return partitionMap->getNumberOfElements();
}
} /* namespace partest */
