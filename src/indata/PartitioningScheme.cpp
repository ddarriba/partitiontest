/*
 * PartitioningScheme.cpp
 *
 *  Created on: Mar 8, 2013
 *      Author: diego
 */

#include "PartitioningScheme.h"
#include "PartitionMap.h"
#include <stdlib.h>

#define UNDEFINED -1

namespace partest {

struct comparePartitionElements {
	inline bool operator()(PartitionElement * e1,
			PartitionElement * e2) {

		if (e1 == e2) return 0;

		t_partitionElementId i1 = e1->getId();
		t_partitionElementId i2 = e2->getId();

		if (i1 == i2) {
			cerr << "[ERROR] There are 2 different elements in the set with the same ID" << endl;
			Utilities::exit_partest(EX_SOFTWARE);
		}

		while( (i1 & 1) == (i2 & 1)) {
			i1 >>= 1;
			i2 >>= 1;
		}
		return (i1 & 1);
	}
};

PartitioningScheme::PartitioningScheme(int numberOfElements) :
		numberOfElements(numberOfElements) {
	currentElement = 0;
	partitions = new vector<PartitionElement *>(numberOfElements);
	numberOfBits = 1;
	code = 0;
}

PartitioningScheme::PartitioningScheme(t_partitioningScheme * schemeVector,
		PartitionMap * partitionMap) {
	currentElement = 0;
	code = 0;
	numberOfElements = schemeVector->size();
	partitions = new vector<PartitionElement *>(numberOfElements);
	int i;
#ifdef DEBUG
	cout << "[TRACE] PARTITION (" << numberOfElements << " elements): [ ";
	for (i = 0; i < numberOfElements; i++) {
		cout << schemeVector->at(i) << " ";
	}
	cout << "]" << endl;
#endif
	numberOfBits = 1;

	for (i = 0; i < numberOfElements; i++) {
		int j;
#ifdef DEBUG
		cout << "[TRACE] ADDING ELEMENT "
		<< i+1 << "/" << numberOfElements << " (" << schemeVector->at(i) << ")" << endl;
#endif
		PartitionElement * newPE = partitionMap->getPartitionElement(
				schemeVector->at(i));

		if (Utilities::binaryLog(newPE->getId() | 1) > numberOfBits) {
			numberOfBits = Utilities::binaryLog(newPE->getId() | 1);
		}
		addElement(newPE);
	}
}

PartitioningScheme::~PartitioningScheme() {
	delete partitions;
	if (code)
		delete code;
}

int PartitioningScheme::addElement(PartitionElement * element) {
	if (currentElement < numberOfElements) {
		partitions->at(currentElement++) = element;

		if (Utilities::binaryLog(element->getId() | 1) > numberOfBits) {
			numberOfBits = Utilities::binaryLog(element->getId() | 1);
		}

		if (currentElement == numberOfElements) {
			std::sort(partitions->begin(), partitions->end(),
						comparePartitionElements());
			for (int i=0;i<numberOfElements;i++) {
				cout << " " << getElement(i)->getName() << "[" << getElement(i)->getId() << "]";
			}
			cout << endl;
		}
		return 0;
	} else {
		return 1;
	}
}

PartitionElement * PartitioningScheme::getElement(int id) {
	if (id >= 0 & id < numberOfElements) {
		return partitions->at(id);
	} else {
		return 0;
	}
}

int PartitioningScheme::isOptimized(void) {
	int optimized = 1;
	for (int i = 0; i < numberOfElements; i++) {
		optimized &= partitions->at(i)->isOptimized();
	}
	return optimized;
}

string PartitioningScheme::toString() {
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
