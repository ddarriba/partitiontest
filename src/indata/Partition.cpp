/*
 * Partition.cpp
 *
 *  Created on: Mar 8, 2013
 *      Author: diego
 */

#include "Partition.h"
#include "PartitionMap.h"
#include <stdlib.h>

#define UNDEFINED -1

namespace partest {

Partition::Partition(int numberOfElements) :
		numberOfElements(numberOfElements) {
	currentElement = 0;
	partitions = (PartitionElement **) malloc(
			numberOfElements * sizeof(PartitionElement *));
	numberOfBits = 0;
	code = 0;
}

Partition::Partition(t_partition_elements * partition,
		PartitionMap * partitionMap) {
	currentElement = 0;
	code = 0;
	numberOfElements = partition->size();
	partitions = (PartitionElement **) malloc(
			numberOfElements * sizeof(PartitionElement *));
	int i;
#ifdef DEBUG
	cout << "[TRACE] PARTITION (" << numberOfElements << " elements): [ ";
	for (i = 0; i < numberOfElements; i++) {
		cout << partition->at(i) << " ";
	}
	cout << "]" << endl;
#endif
	numberOfBits = 0;
	for (i = 0; i < numberOfElements; i++) {
		int j;
#ifdef DEBUG
		cout << "[TRACE] ADDING ELEMENT "
		<< i+1 << "/" << numberOfElements << " (" << partition->at(i) << ")" << endl;
#endif
		PartitionElement * newPE = partitionMap->getPartitionElement(
				partition->at(i));
		if (Utilities::binaryLog(newPE->getId() | 1) > numberOfBits) {
			numberOfBits = Utilities::binaryLog(newPE->getId() | 1);
		}
		addElement(newPE);
	}
}

Partition::~Partition() {
	free(partitions);
	if (code)
		delete code;
}

int Partition::addElement(PartitionElement * element) {
	if (currentElement < numberOfElements) {
		partitions[currentElement++] = element;
		if (Utilities::binaryLog(element->getId() | 1) > numberOfBits) {
			numberOfBits = Utilities::binaryLog(element->getId() | 1);
		}
		return 0;
	} else {
		return 1;
	}
}

PartitionElement * Partition::getElement(int id) {
	if (id >= 0 & id < numberOfElements) {
		return partitions[id];
	} else {
		return 0;
	}
}

string Partition::toString() {
	if (!code) {
		int hashmap[numberOfElements];
		int intcode[numberOfBits];
		char charcode[numberOfBits + 1];
		int i, j;
		for (i = 0; i < numberOfElements; i++) {
			t_partitionElementId id = getElement(i)->getId();
			hashmap[i] = UNDEFINED;
			for (j = 0; j < numberOfBits; j++) {
				if (id & Utilities::binaryPow(j)) {
					intcode[j] = i;
				}
			}
		}

		hashmap[intcode[0]] = 0;
		int nextCode = 0;
		charcode[0] = '0';
		for (i = 1; i < numberOfBits; i++) {
			if (hashmap[intcode[i]] == UNDEFINED) {
				hashmap[intcode[i]] = ++nextCode;
			}
			charcode[i] = hashmap[intcode[i]] + '0';
		}
		charcode[numberOfBits] = '\0';
		code = new string(charcode);
	}
	return *code;
}
} /* namespace partest */
