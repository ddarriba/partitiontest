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
#include <sstream>

namespace partest {

PartitionMap::PartitionMap(const char * configFile, Alignment * alignment,
		bitMask rateVariation, DataType dataType) :
		alignment(alignment), rateVariation(rateVariation), dataType(dataType) {
#ifdef DEBUG
	cout << "[TRACE] Instantiated partition map from config file" << endl;
#endif

	if (strcmp(configFile, "")) {

		ConfigParser parser(configFile);

		partitions = new vector<partitionMappingInfo>(
				parser.getNumberOfPartitions());
		numberOfElements = numberOfPartitions = parser.getNumberOfPartitions();
		for (int i = 0; i < parser.getNumberOfPartitions(); i++) {
			struct partitionInfo nextPartition = parser.getPartition(i);
			partitions->at(i).partitionId = nextPartition.partitionId;

			partitions->at(i).partitionElement = new PartitionElement(
					nextPartition.partitionId, nextPartition.name, alignment,
					nextPartition.start, nextPartition.end,
					nextPartition.stride, rateVariation, dataType);
		}

	} else {
		partitions = new vector<partitionMappingInfo>(1);
		string nameStr("UNIQUE");
		int start = 1;
		int end = alignment->getNumSites();
		int stride = 0;

		partitions->at(0).partitionId.push_back(1);
		partitions->at(0).partitionElement = new PartitionElement(
				partitions->at(0).partitionId, nameStr, alignment, start, end,
				stride, rateVariation, dataType);
		numberOfElements = numberOfPartitions = 1;
	}
}

PartitionMap::PartitionMap(Alignment * alignment,
		unsigned int numberOfPartitions, bitMask rateVariation,
		DataType dataType) :
		alignment(alignment), numberOfElements(numberOfPartitions), numberOfPartitions(
				numberOfPartitions), rateVariation(rateVariation), dataType(
				dataType) {

#ifdef DEBUG
	cout << "[TRACE] Instantiated partition map from void" << endl;
#endif

	partitions = partitions = new vector<partitionMappingInfo>(
			numberOfElements);
}

PartitionMap::~PartitionMap() {
	int i;
	for (i = 0; i < numberOfElements; i++) {
		delete partitions->at(i).partitionElement;
	}
	delete partitions;
}

bool PartitionMap::addPartitionElement(unsigned int partitionId, string name,
		unsigned int startPosition, unsigned int endPosition, char stride) {

	/* only gene-partitions can be added */
	if (partitionId >= numberOfElements)
		return false;

	/* partitionId is translated into partition mask */
	partitions->at(partitionId).partitionId.push_back(partitionId);
	partitions->at(partitionId).partitionElement = new PartitionElement(
			partitions->at(partitionId).partitionId, name, alignment,
			startPosition, endPosition, stride, rateVariation, dataType);

	return true;
}

PartitionElement * PartitionMap::getPartitionElement(unsigned int id) {
	assert(id < numberOfPartitions);
	return partitions->at(id).partitionElement;
}

PartitionElement * PartitionMap::getPartitionElement(
		t_partitionElementId partitionId) {

	assert(partitionId.size() > 0);

	if (partitionId.size() == 1) {
		return partitions->at(partitionId.at(0)).partitionElement;
	}

	/* search for element */
	for (int i = 0; i < partitions->size(); i++) {
		partitionMappingInfo pInfo = partitions->at(i);
		if (pInfo.partitionId == partitionId) {
			return pInfo.partitionElement;
		}
	}

	/* we need to merge partitions */
	int numberOfSections = partitionId.size();
	int * start = (int *) malloc(numberOfSections * sizeof(int));
	int * end = (int *) malloc(numberOfSections * sizeof(int));
	int * stride = (int *) malloc(numberOfSections * sizeof(int));
	int curIndex = 0;
	stringstream name;
	name << "( ";
	for (int i = 0; i < numberOfSections; i++) {
		start[curIndex] = getPartitionElement(partitionId.at(i))->getStart();
		end[curIndex] = getPartitionElement(partitionId.at(i))->getEnd();
		stride[curIndex] = getPartitionElement(partitionId.at(i))->getStride();
		name << getPartitionElement(partitionId.at(i))->getName() << " ";
		curIndex++;
	}
	name << ")";

	partitionMappingInfo pInfo;
	pInfo.partitionId = partitionId;

	pInfo.partitionElement = new PartitionElement(partitionId, name.str(),
			alignment, start, end, stride, numberOfSections, rateVariation,
			dataType);

	partitions->push_back(pInfo);
	numberOfElements++;
	return pInfo.partitionElement;
}

void PartitionMap::deletePartitionElement(t_partitionElementId id) {

	if (id.size() > 1) {
		for (int i = numberOfPartitions; i < numberOfElements; i++) {
			if (partitions->at(i).partitionId == id) {
				delete partitions->at(i).partitionElement;
				partitions->erase(partitions->begin() + i);
				numberOfElements--;
				return;
			}
		}
		cerr << "[ERROR] Attempting to delete an inexistent partition element"
				<< endl;
		Utilities::exit_partest(EX_SOFTWARE);
	}
}

} /* namespace partest */
