/*
 * PartitioningScheme.cpp
 *
 *  Created on: Mar 8, 2013
 *      Author: diego
 */

#include "PartitioningScheme.h"
#include "util/Utilities.h"
#include "indata/PartitionMap.h"

#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>

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
			exit_partest(EX_SOFTWARE);
		}

		return i1 < i2;
	}
};

PartitioningScheme::PartitioningScheme(t_partitioningScheme * schemeVector) :
		partitions(schemeVector->size()), id(schemeVector->size()) {

	eps = 0;
	tree = 0;
	currentElement = 0;
	code = 0;
	sort(schemeVector->begin(), schemeVector->end());
	numberOfElements = schemeVector->size();
	codeLines = Utilities::iDecLog(numberOfElements - 1) + 1;

	for (unsigned int i = 0; i < schemeVector->size(); i++) {
		id.at(i) = schemeVector->at(i);
		partitions.at(i) = PartitionMap::getInstance()->getPartitionElement(
				schemeVector->at(i));
	}
}

PartitioningScheme::~PartitioningScheme() {
	if (code)
		delete code;
	if (tree)
		free(tree);
	if (eps) {
		for (elementPair *ep : *eps) {
			free(ep);
		}
		delete eps;
	}
}

PartitionElement * PartitioningScheme::getElement(unsigned int id) {
	if ((id >= 0) & (id < numberOfElements)) {
		return partitions.at(id);
	} else {
		return 0;
	}
}

bool PartitioningScheme::isOptimized(void) {
	bool optimized = true;
	for (PartitionElement * pe : partitions) {
		optimized &= pe->isOptimized();
	}
	return optimized;
}

/** Functor for sorting the pairs */
struct compareDistancesVector {
	inline bool operator()(elementPair * struct1, elementPair * struct2) {
		return (struct1->distance < struct2->distance);
	}
};

vector<elementPair *> * PartitioningScheme::getElementDistances() {

	if (!eps) {

		if (!isOptimized()) {
			cerr
					<< "[ERROR] Attempting to get differences of unoptimized partitions"
					<< endl;
			exit_partest(EX_SOFTWARE);
		}

		eps = new vector<elementPair *>(
				(numberOfElements * (numberOfElements - 1)) / 2);

		vector<double> alphaValues(
				numberOfElements * (numberOfElements + 1) / 2);
#ifdef _IG_MODELS
		vector<double> pinvValues(
				numberOfElements * (numberOfElements + 1) / 2);
#endif
		vector<double> distances(numberOfElements * (numberOfElements + 1) / 2);
		double sumAlpha = 0.0;
		double maxAlpha = 0.0;
		for (unsigned int i = 1; i < numberOfElements; i++) {
			for (unsigned int j = 0; j < i; j++) {
				Model * mi = getElement(i)->getBestModel()->getModel();
				Model * mj = getElement(j)->getBestModel()->getModel();
				int index = (i * (i - 1) / 2) + j;
#ifdef _IG_MODELS
				pinvValues.at(index) = pow(mi->getpInv() - mj->getpInv(), 2);
#endif
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
				distances.at(index) += alphaValues[index];
#ifdef _IG_MODELS
				distances.at(index) += pinvValues[index]
#endif
				                                                ;
				elementPair * ep = (elementPair *) malloc(sizeof(elementPair));
				ep->e1 = getElement(i);
				ep->e2 = getElement(j);
				ep->distance = distances.at(index);
				eps->at(index) = ep;

				if (distances.at(index) < minDistance) {
					//el1 = getElement(i)->getId();
					//el2 = getElement(j)->getId();
					minDistance = distances.at(index);
				}
			}
		}
		sort(eps->begin(), eps->end(), compareDistancesVector());
	}

	return eps;
}

string PartitioningScheme::getName() {
	stringstream ss;
	for (PartitionElement * pe : partitions) {
		ss << pe->getName() << " ";
	}
	return ss.str();
}

int PartitioningScheme::getCodeLines(void) {
	return codeLines;
}

string PartitioningScheme::getCode(int codeLine) {
	if (!code) {
		if (codeLine >= (int) codeLines)
			exit_partest(EX_SOFTWARE);

		int hashmap[numberOfElements];
		int intcode[number_of_genes];
		char charcode[codeLines * number_of_genes + 1];
		for (unsigned int i = 0; i < numberOfElements; i++) {
			t_partitionElementId id = getElement(i)->getId();
			hashmap[i] = UNDEFINED;
			for (unsigned int j = 0; j < id.size(); j++) {
				intcode[id.at(j)] = i;
			}
		}
		hashmap[intcode[0]] = 0;
		int nextCode = 0;
		for (unsigned int i = 0; i < codeLines; i++) {
			charcode[i * number_of_genes + 0] = '0';
		}
		for (unsigned int i = 1; i < number_of_genes; i++) {
			if (hashmap[intcode[i]] == UNDEFINED) {
				hashmap[intcode[i]] = ++nextCode;
			}
			for (unsigned int j = 0; j < codeLines; j++) {
				charcode[j * number_of_genes + i] = ((int) floor(
						hashmap[intcode[i]] / pow(10, j)) % 10) + '0';
			}
		}
		charcode[codeLines * number_of_genes] = '\0';

		code = new string(charcode);
	}

	if (codeLine == FULL_CODE) {
		return *code;
	} else {
		string the_code = code->substr(codeLine * number_of_genes,
				number_of_genes);
		return the_code;
	}
}

void PartitioningScheme::setTree(char * tree) {
	this->tree = tree;
}

char * PartitioningScheme::getTree(void) const {
	return tree;
}

double PartitioningScheme::getLnL() {
	if (!isOptimized())
		return 0.0;
	double lk = 0.0;
	for (PartitionElement * pe : partitions) {
		lk += pe->getLnL();
	}
	return lk;
}

double PartitioningScheme::getIcValue() {
	if (!isOptimized())
		return 0.0;
	double score = 0.0;
	for (PartitionElement * pe : partitions) {
		score += pe->getBestModel()->getValue();
	}
	return score;
}

void PartitioningScheme::print(ostream & out) {
	out << getName() << endl;
	out << "Code Lines: " << codeLines << endl;
	for (int i = codeLines-1; i >= 0; i--) {
		out << "  " << getCode(i) << endl;
	}
	out << "Num.Elements: " << numberOfElements << endl;
	if (isOptimized()) {
		out << "LnL:          " << getLnL() << endl;
		out << "BIC_score:    " << getIcValue() << endl;
		if (tree > 0) {
			out << "Tree:         " << getTree() << endl;
		}
		out << endl;
		out << setw(6) << "--ID-"
			<< setw(11) << "---Model-- "
							<< setw(5) << " --K-"
							<< setw(21) << " ---------LnL--------"
							<< setw(21) << " ---------BIC--------" << endl;
		for (unsigned int i=0; i<numberOfElements; i++) {
			out << setw(6) << left << i + 1
					<< setw(11) << left << getElement(i)->getBestModel()->getModel()->getName()
					<< setw(5) << right << getElement(i)->getBestModel()->getModel()->getNumberOfFreeParameters()
					<< setw(21) << right << fixed << setprecision(4) << getElement(i)->getBestModel()->getModel()->getLnL()
					<< setw(21) << right << fixed << setprecision(4) << getElement(i)->getBestModel()->getValue() << endl;
		}
	}
}

} /* namespace partest */
