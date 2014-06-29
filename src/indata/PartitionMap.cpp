/*
 * PartitionMap.cpp
 *
 *  Created on: Feb 13, 2013
 *      Author: diego
 */

#include "PartitionMap.h"
#include "util/Utilities.h"
#include <assert.h>
#include <string.h>
#include <iostream>
#include <algorithm>

using namespace std;

namespace partest {

PartitionMap::PartitionMap() {
		_keep.reserve(number_of_genes);
		partitions = new vector<partitionMappingInfo>(pllPartitions->numberOfPartitions);
		numberOfElements = numberOfPartitions = pllPartitions->numberOfPartitions;
		for (unsigned int i = 0; i < numberOfPartitions; i++) {
			t_partitionElementId nextId;
			nextId.push_back(i);
			partitions->at(i).partitionId = nextId;
			partitions->at(i).partitionElement = new PartitionElement(nextId);
		}
}

PartitionMap::~PartitionMap() {
	for (unsigned int i = 0; i < numberOfElements; i++) {
		delete partitions->at(i).partitionElement;
	}
	delete partitions;
}

PartitionElement * PartitionMap::getPartitionElement(
		t_partitionElementId partitionId) {

	if (!partitionId.size()) {
		cerr << "[ERROR] [PartitionMap] Partition ID is empty" << endl;
		exit_partest(EX_SOFTWARE);
	}

	if (partitionId.size() == 1) {
		return partitions->at(partitionId.at(0)).partitionElement;
	}

	sort(partitionId.begin(), partitionId.end());

	/* search for element */
	for (unsigned int i = 0; i < partitions->size(); i++) {
		partitionMappingInfo pInfo = partitions->at(i);
		if (pInfo.partitionId == partitionId) {
			return pInfo.partitionElement;
		}
	}

	/* we need to merge partitions */
	partitionMappingInfo pInfo;
	pInfo.partitionId = partitionId;

	pInfo.partitionElement = new PartitionElement(partitionId);

	partitions->push_back(pInfo);
	numberOfElements++;

	return pInfo.partitionElement;
}

void PartitionMap::deletePartitionElement(t_partitionElementId id) {

		for (unsigned int i = numberOfPartitions; i < numberOfElements; i++) {
			if (partitions->at(i).partitionId == id) {
				delete partitions->at(i).partitionElement;
				partitions->erase(partitions->begin() + i);
				numberOfElements--;
				return;
			}
		}
		cerr << "[ERROR] [PartitionMap] Attempting to delete an inexistent partition element"
				<< endl;
		exit_partest(EX_SOFTWARE);
}

void PartitionMap::purgePartitionMap(t_partitionElementId id) {

	if (id.size() <= 1)
		return;

	/* Loop backwards over complex elements */
	for (unsigned int i = numberOfElements - 1; i >= numberOfPartitions; i--) {
		t_partitionElementId mId = partitions->at(i).partitionId;
		if (mId != id && Utilities::intersec(mId, id) && !Utilities::contains(_keep, mId)) {
			deletePartitionElement(mId);
		}
	}
}

void PartitionMap::keep(t_partitioningScheme id) {
	_keep.clear();
	for (t_partitionElementId eId : id) {
		_keep.push_back(eId);
	}
}

void PartitionMap::keep_add(t_partitioningScheme id) {
	for (t_partitionElementId eId : id) {
		if (!Utilities::contains(_keep, eId)) {
			_keep.push_back(eId);
		}
	}
}

void PartitionMap::keep_remove(t_partitioningScheme id) {
	int position = 0;
	for (t_partitionElementId eId : id) {
		if (Utilities::contains(_keep, eId)) {
			_keep.erase(_keep.begin() + position);
		}
		position++;
	}
}

PartitionMap * PartitionMap::instance = 0;

PartitionMap * PartitionMap::getInstance() {
	if (!instance) {
		instance = new PartitionMap();
	}
	return instance;
}

void PartitionMap::deleteInstance() {
	if (instance)
		delete instance;
	instance = 0;
}

} /* namespace partest */
