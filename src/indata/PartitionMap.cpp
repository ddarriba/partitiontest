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
		partitions = new vector<partitionMappingInfo>((size_t) pllPartitions->numberOfPartitions);
		numberOfElements = numberOfPartitions = (size_t) pllPartitions->numberOfPartitions;
		for (size_t i = 0; i < numberOfPartitions; i++) {
			t_partitionElementId nextId;
			nextId.push_back(i);
			partitions->at(i).partitionId = nextId;
			partitions->at(i).partitionElement = new PartitionElement(nextId);
		}
}

PartitionMap::~PartitionMap() {
	for (size_t i = 0; i < numberOfElements; i++) {
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
	for (size_t i = 0; i < partitions->size(); i++) {
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

		for (size_t i = numberOfPartitions; i < numberOfElements; i++) {
			if (partitions->at(i).partitionId == id) {
				delete partitions->at(i).partitionElement;
				partitions->erase(partitions->begin() + (long) i);
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
	for (size_t i = numberOfElements - 1; i >= numberOfPartitions; i--) {
		t_partitionElementId mId = partitions->at(i).partitionId;
		if (mId != id && Utilities::intersec(mId, id) && !Utilities::contains(_keep, mId)) {
			deletePartitionElement(mId);
		}
	}
}

void PartitionMap::keep(t_partitioningScheme id) {
	_keep.clear();
	for (size_t i=0; i<id.size(); i++) {
		t_partitionElementId eId = id.at(i);
		_keep.push_back(eId);
	}
}

void PartitionMap::keep_add(t_partitioningScheme id) {
	for (size_t i=0; i<id.size(); i++) {
		t_partitionElementId eId = id.at(i);
		if (!Utilities::contains(_keep, eId)) {
			_keep.push_back(eId);
		}
	}
}

void PartitionMap::keep_remove(t_partitioningScheme id) {
	int position = 0;
	for (size_t i=0; i<id.size(); i++) {
		t_partitionElementId eId = id.at(i);
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
