/*
 * PartitionMap.cpp
 *
 *  Created on: Feb 13, 2013
 *      Author: diego
 */

#include "PartitionMap.h"
#include "../util/Utilities.h"
#include <assert.h>
#include <string.h>
#include <sstream>

namespace partest {

PartitionMap::PartitionMap(const char * configFile, Alignment * alignment,
		bitMask rateVariation, DataType dataType) :
		alignment(alignment), rateVariation(
				rateVariation), dataType(dataType) {
#ifdef DEBUG
	cout << "[TRACE] Instantiated partition map from config file" << endl;
#endif

	FILE *f;
	char *cc = (char *) NULL;
	int nbytes;
	int partitionId = 0;
	numberOfElements = 0;
	//int **partitions;

	f = Utilities::myfopen(configFile, "rb");

	while (Utilities::myGetline(&cc, &nbytes, f) > -1) {

		numberOfElements++;

		if (cc)
			free(cc);
		cc = (char *) NULL;
	}

	rewind(f);

	assert(
			Utilities::binaryPow(numberOfElements)
					< sizeof(t_partitionElementId) * 8);
	partitions = new vector<partitionInfo>(numberOfElements);

	while (Utilities::myGetline(&cc, &nbytes, f) > -1) {

		char * name = strtok(cc, "=");
		string nameStr(name);
		int start = atoi(strtok(NULL, ","));
		int end = atoi(strtok(NULL, "\\"));
		char * strideStr = strtok(NULL, "\\");
		int stride = strideStr ? atoi(strideStr) : 0;

		/* partitionId is translated into partition mask */
		partitions->at(partitionId).partitionId = Utilities::binaryPow(
				partitionId);
		partitions->at(partitionId).partitionElement = new PartitionElement(
				partitions->at(partitionId).partitionId, nameStr, alignment, start, end,
				stride, rateVariation, dataType);

		partitionId++;

		if (cc)
			free(cc);

		cc = (char *) NULL;
	}

	numberOfPartitions = numberOfElements;

}

PartitionMap::PartitionMap(Alignment * alignment, unsigned int numberOfPartitions,
		bitMask rateVariation, DataType dataType) :
		alignment(alignment), numberOfElements(numberOfPartitions), numberOfPartitions(numberOfPartitions), rateVariation(
				rateVariation), dataType(dataType) {

#ifdef DEBUG
	cout << "[TRACE] Instantiated partition map from void" << endl;
#endif

	/* check whether the number of partitions exceed the number of bits of the id mask length */
	assert(numberOfElements < sizeof(t_partitionElementId) * 8);

	partitions = partitions = new vector<partitionInfo>(numberOfElements);
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
			partitions->at(partitionId).partitionId, name, alignment, startPosition,
			endPosition, stride, rateVariation, dataType);

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
			name << getPartitionElement(Utilities::binaryPow(i))->getName() << " ";
			curIndex++;
		}
	}
	name << ")";

	partitionInfo pInfo;
	pInfo.partitionId = partitionId;
// t_partitionElementId id, Alignment * alignment,
// int * start, int * end, int * stride, int numberOfSections,
// bitMask rateVariation, DataType dataType

	pInfo.partitionElement = new PartitionElement(partitionId, name.str(), alignment, start,
			end, stride, numberOfSections, rateVariation, dataType);
	partitions->push_back(pInfo);
	numberOfElements++;
	return pInfo.partitionElement;
}

} /* namespace partest */
