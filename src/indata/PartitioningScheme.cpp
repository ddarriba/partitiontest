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
	inline bool operator()(PartitionElement * e1, PartitionElement * e2) {

		if (e1 == e2)
			return 0;

		t_partitionElementId i1 = e1->getId();
		t_partitionElementId i2 = e2->getId();

		if (i1 == i2) {
			cerr
					<< "[ERROR] There are 2 different elements in the set with the same ID"
					<< endl;
			Utilities::exit_partest(EX_SOFTWARE);
		}

		return i1 < i2;
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
	unsigned int i;
#ifdef DEBUG
	cout << "[TRACE] PARTITION (" << numberOfElements << " elements):";
#endif
	numberOfBits = 1;

	for (i = 0; i < numberOfElements; i++) {
#ifdef DEBUG
		cout << "[TRACE] ADDING ELEMENT "
		<< i+1 << "/" << numberOfElements << endl; //" (" << schemeVector->at(i) << ")" << endl;
#endif
		PartitionElement * newPE = partitionMap->getPartitionElement(
				schemeVector->at(i));

		if (newPE->getMaxId() >= numberOfBits) {
			numberOfBits = (newPE->getMaxId() + 1);
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

		t_partitionElementId id(element->getId());
		// Security checks
		for (unsigned int i = 0; i < currentElement; i++) {
			t_partitionElementId id_c(partitions->at(i)->getId());
			if (id_c == id) {
				cerr << "[Error] Duplicated element" << endl;
				Utilities::vprint(cerr, id);
				Utilities::exit_partest(EX_SOFTWARE);
			}
			if (Utilities::intersec(id_c, id)) {
				cerr << "[Error] Overlapped elements" << endl;
				t_partitionElementId id1 = element->getId();
				t_partitionElementId id2 = partitions->at(i)->getId();
				Utilities::vprint(cerr, id1);
				Utilities::vprint(cerr, id2);
				Utilities::exit_partest(EX_SOFTWARE);
			}
		}

		partitions->at(currentElement++) = element;

		if (element->getMaxId() >= numberOfBits) {
			numberOfBits = (element->getMaxId() + 1);
		}
		if (currentElement == numberOfElements) {
			std::sort(partitions->begin(), partitions->end(),
					comparePartitionElements());
		}
		return 0;
	} else {
		return 1;
	}
}

PartitionElement * PartitioningScheme::getElement(unsigned int id) {
	if ((id >= 0) & (id < numberOfElements)) {
		return partitions->at(id);
	} else {
		return 0;
	}
}

bool PartitioningScheme::isOptimized(void) {
	bool optimized = true;
	for (unsigned int i = 0; i < numberOfElements; i++) {
		optimized &= partitions->at(i)->isOptimized();
	}
	return optimized;
}

void PartitioningScheme::resetModelSet() {
	for (unsigned int i = 0; i < numberOfElements; i++) {
		for (unsigned int j = 0;
				j < partitions->at(i)->getModelset()->getNumberOfModels();
				j++) {
			partitions->at(i)->getModelset()->getModel(j)->setLnL(0.0);
		}
	}
}

void PartitioningScheme::buildCompleteModelSet(bool clearAll) {
	for (unsigned int i = 0; i < numberOfElements; i++) {
		partitions->at(i)->buildCompleteModelSet(clearAll);
	}
}

void PartitioningScheme::getClosestPartitions(t_partitionElementId & el1,
		t_partitionElementId & el2) {
	if (!isOptimized()) {
		cerr
				<< "[ERROR] Attempting to get differences of unoptimized partitions"
				<< endl;
		Utilities::exit_partest(EX_SOFTWARE);
	}

	vector<double> alphaValues(numberOfElements * (numberOfElements + 1) / 2);
	vector<double> pinvValues(numberOfElements * (numberOfElements + 1) / 2);
	vector<double> distances(numberOfElements * (numberOfElements + 1) / 2);
	double sumAlpha = 0.0;
	double maxAlpha = 0.0;
	for (unsigned int i = 1; i < numberOfElements; i++) {
		for (unsigned int j = 0; j < i; j++) {
			Model * mi = getElement(i)->getBestModel()->getModel();
			Model * mj = getElement(j)->getBestModel()->getModel();
			int index = (i * (i - 1) / 2) + j;
			pinvValues.at(index) = pow(mi->getpInv() - mj->getpInv(), 2);
			alphaValues.at(index) = pow(mi->getAlpha() - mj->getAlpha(), 2);
			sumAlpha += alphaValues.at(index);
			if (alphaValues.at(index) > maxAlpha) {
				maxAlpha = alphaValues.at(index);
			}
			distances.at(index) = mi->distanceTo(mj);
		}
	}
	double minDistance = DOUBLE_INF;
	for (unsigned int i = 1; i < numberOfElements; i++) {
		for (unsigned int j = 0; j < i; j++) {
			int index = (i * (i - 1) / 2) + j;
			//alphaValues.at(index) /= maxAlpha;
			distances.at(index) += alphaValues[index] + pinvValues[index];
			if (distances.at(index) < minDistance) {
				el1 = getElement(i)->getId();
				el2 = getElement(j)->getId();
				minDistance = distances.at(index);
			}
		}
	}
}

string PartitioningScheme::getName() {
	stringstream ss;
	for (unsigned int i = 0; i < numberOfElements; i++) {
		ss << partitions->at(i)->getName() << " ";
	}
	return ss.str();
}

string PartitioningScheme::getCode() {
	if (!code) {
		int numDigits =
				(numberOfElements > 10) ?
						(Utilities::iDecLog(numberOfElements - 1) + 1) : 1;

		int hashmap[numberOfElements];
		int intcode[numberOfBits];
		char charcode[numDigits * numberOfBits + 1];
		for (unsigned int i = 0; i < numberOfElements; i++) {
			t_partitionElementId id = getElement(i)->getId();
			hashmap[i] = UNDEFINED;
			for (unsigned int j = 0; j < id.size(); j++) {
				intcode[id.at(j)] = i;
			}
		}
		hashmap[intcode[0]] = 0;
		int nextCode = 0;
		for (int i = 0; i < numDigits; i++) {
			charcode[i * numberOfBits + 0] = '0';
		}
		for (unsigned int i = 1; i < numberOfBits; i++) {
			if (hashmap[intcode[i]] == UNDEFINED) {
				hashmap[intcode[i]] = ++nextCode;
			}
			for (int j = 0; j < numDigits; j++) {
				charcode[j * numberOfBits + i] = ((int) floor(
						hashmap[intcode[i]] / pow(10, j)) % 10) + '0';
			}
		}
		for (int j = 0; j < numDigits; j++) {
			charcode[j * numberOfBits + numberOfBits] = '\n';
		}
		charcode[numDigits * numberOfBits] = '\0';

		code = new string(charcode);

	}
	return *code;
}

double PartitioningScheme::getLnL() {
	if (!isOptimized())
		return 0.0;
	double lk = 0.0;
	for (unsigned int i = 0; i < numberOfElements; i++) {
		lk += partitions->at(i)->getBestModel()->getModel()->getLnL();
	}
	return lk;
}
} /* namespace partest */
