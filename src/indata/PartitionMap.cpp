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

		partitions->at(0).partitionId = 1;
		partitions->at(0).partitionElement = new PartitionElement(1, nameStr,
				alignment, start, end, stride, rateVariation, dataType);
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

	/* check whether the number of partitions exceed the number of bits of the id mask length */
	assert(numberOfElements < sizeof(t_partitionElementId) * 8);

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

	if (partitionId >= numberOfElements)
		return false;

	/* partitionId is translated into partition mask */
	partitions->at(partitionId).partitionId = Utilities::binaryPow(partitionId);
	partitions->at(partitionId).partitionElement = new PartitionElement(
			partitions->at(partitionId).partitionId, name, alignment,
			startPosition, endPosition, stride, rateVariation, dataType);

	return true;
}

PartitionElement * PartitionMap::getPartitionElement(
		t_partitionElementId partitionId) {

	assert(partitionId > 0);

	int i;
	for (i = 0; i < numberOfElements; i++) {
		// TODO: Hashmap instead of iteration?
		if (partitions->at(i).partitionId == partitionId) {
			return partitions->at(i).partitionElement;
		}
	}

	/* element not found */
	if (Utilities::isPowerOfTwo(partitionId)) {
		/* single elements should be already mapped */
		cerr << "Requested element (" << partitionId << ") does not exist"
				<< endl;
		exit(-1);
	}

	/* we need to merge partitions */
	int numberOfSections = Utilities::setbitsCount(partitionId);
	int * start = (int *) malloc(numberOfSections * sizeof(int));
	int * end = (int *) malloc(numberOfSections * sizeof(int));
	int * stride = (int *) malloc(numberOfSections * sizeof(int));
	int curIndex = 0;
	stringstream name;
	name << "( ";
	for (i = 0; i < Utilities::binaryLog(partitionId); i++) {
		if (partitionId & Utilities::binaryPow(i)) {
			start[curIndex] =
					getPartitionElement(Utilities::binaryPow(i))->getStart();
			end[curIndex] =
					getPartitionElement(Utilities::binaryPow(i))->getEnd();
			stride[curIndex] =
					getPartitionElement(Utilities::binaryPow(i))->getStride();
			name << getPartitionElement(Utilities::binaryPow(i))->getName()
					<< " ";
			curIndex++;
		}
	}
	name << ")";

	partitionMappingInfo pInfo;
	pInfo.partitionId = partitionId;
// t_partitionElementId id, Alignment * alignment,
// int * start, int * end, int * stride, int numberOfSections,
// bitMask rateVariation, DataType dataType

	pInfo.partitionElement = new PartitionElement(partitionId, name.str(),
			alignment, start, end, stride, numberOfSections, rateVariation,
			dataType);
	partitions->push_back(pInfo);
	numberOfElements++;
	return pInfo.partitionElement;
}

} /* namespace partest */
