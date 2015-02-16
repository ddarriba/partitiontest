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

#include "PartitioningScheme.h"
#include "util/Utilities.h"
#include "exe/ModelSelector.h"
#include "indata/PartitionMap.h"

#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>

#define UNDEFINED -1

using namespace std;

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
		id(schemeVector->size()), partitions(schemeVector->size()) {

	eps = 0;
	tree = 0;
	currentElement = 0;
	code = 0;
	sort(schemeVector->begin(), schemeVector->end());
	numberOfElements = schemeVector->size();
	codeLines = (size_t) Utilities::iDecLog((int) numberOfElements - 1) + 1;

	for (size_t i = 0; i < schemeVector->size(); i++) {
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
		for (size_t i=0; i<eps->size(); i++) {
			elementPair *ep = eps->at(i);
			free(ep);
		}
		delete eps;
	}
}

PartitionElement * PartitioningScheme::getElement(size_t _id) {
	if (_id < numberOfElements) {
		return partitions.at(_id);
	} else {
		return 0;
	}
}

bool PartitioningScheme::isOptimized(void) {
	bool optimized = true;
	for (size_t i=0; i<partitions.size(); i++) {
		PartitionElement * pe = partitions.at(i);
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
		for (size_t i = 1; i < numberOfElements; i++) {
			for (size_t j = 0; j < i; j++) {
				Model * mi = getElement(i)->getBestModel()->getModel();
				Model * mj = getElement(j)->getBestModel()->getModel();
				size_t index = (i * (i - 1) / 2) + j;
#ifdef _IG_MODELS
				pinvValues.at(index) = pow(mi->getpInv() - mj->getpInv(), 2);
#endif
				alphaValues.at(index) = pow(mi->getAlpha() - mj->getAlpha(), 2);
				sumAlpha += alphaValues.at(index);
				if (alphaValues.at(index) > maxAlpha) {
					maxAlpha = alphaValues.at(index);
				}
				distances.at(index) = mi->distanceTo(mj);
				if (starting_topology == StartTopoFIXED) {
					double treeDist = 0.0;
					double * bl1 = getElement(i)->getBranchLengths();
					double * bl2 = getElement(j)->getBranchLengths();
					for (int cBranch=0; cBranch < Utilities::numberOfBranches((int) num_taxa); cBranch++) {
						treeDist += pow(bl1[cBranch] - bl2[cBranch], 2);
					}
					distances.at(index) += treeDist;
				}
			}
		}
		double minDistance = DOUBLE_INF;
		for (size_t i = 1; i < numberOfElements; i++) {
			for (size_t j = 0; j < i; j++) {
				size_t index = (i * (i - 1) / 2) + j;
				alphaValues.at(index) /= maxAlpha;
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
	for (size_t i=0; i<partitions.size(); i++) {
	PartitionElement * pe = partitions.at(i);
		ss << pe->getName() << " ";
	}
	return ss.str();
}

int PartitioningScheme::getCodeLines(void) {
	return (int) codeLines;
}

string PartitioningScheme::getCode(int codeLine) {
	if (!code) {
		if (codeLine >= (int) codeLines)
			exit_partest(EX_SOFTWARE);

		int * hashmap = (int *) malloc(numberOfElements * sizeof(int));
		int * intcode = (int *) malloc(number_of_genes * sizeof(int));
		char * charcode = (char *) malloc(codeLines * number_of_genes + 1);
		for (size_t i = 0; i < numberOfElements; i++) {
			t_partitionElementId _id = getElement(i)->getId();
			hashmap[i] = UNDEFINED;
			for (size_t j = 0; j < _id.size(); j++) {
				intcode[_id.at(j)] = (int) i;
			}
		}
		hashmap[intcode[0]] = 0;
		int nextCode = 0;
		for (size_t i = 0; i < codeLines; i++) {
			charcode[i * number_of_genes + 0] = '0';
		}
		for (size_t i = 1; i < number_of_genes; i++) {
			if (hashmap[intcode[i]] == UNDEFINED) {
				hashmap[intcode[i]] = ++nextCode;
			}
			for (size_t j = 0; j < codeLines; j++) {
				charcode[j * number_of_genes + i] = ((int) floor(
						hashmap[intcode[i]] / pow(10, j)) % 10) + '0';
			}
		}
		charcode[codeLines * number_of_genes] = '\0';

		code = new string(charcode);

		free(intcode);
		free(hashmap);
		free(charcode);
	}

	if (codeLine == FULL_CODE) {
		return *code;
	} else {
		string the_code = code->substr((size_t) codeLine * number_of_genes,
				number_of_genes);
		return the_code;
	}
}

void PartitioningScheme::setTree(char * _tree) {
	this->tree = _tree;
}

char * PartitioningScheme::getTree(void) const {
	return tree;
}

double PartitioningScheme::getLnL() {
	if (!isOptimized())
		return 0.0;
	double lk = 0.0;
	for (size_t i=0; i<partitions.size(); i++) {
		PartitionElement * pe = partitions.at(i);
		lk += pe->getLnL();
	}
	return lk;
}

unsigned int PartitioningScheme::getNumberOfFreeParameters() {
	if (!isOptimized())
		return 0;
	unsigned int k = reoptimize_branch_lengths?0:(unsigned int) Utilities::numberOfBranches((int) num_taxa);
	for (size_t i=0; i<partitions.size(); i++) {
		PartitionElement * pe = partitions.at(i);
		if (reoptimize_branch_lengths) {
			k += (unsigned int) pe->getBestModel()->getModel()->getNumberOfFreeParameters();
		} else {
			k += (unsigned int) pe->getBestModel()->getModel()->getModelFreeParameters();
		}
	}
	return k;
}

double PartitioningScheme::getIcValue() {
	if (!isOptimized())
		return 0.0;
	if (!reoptimize_branch_lengths) {
		return getLinkedIcValue();
	}
	double score = 0.0;
	for (size_t i = 0; i < partitions.size(); i++) {
		PartitionElement * pe = partitions.at(i);
		score += pe->getBestModel()->getValue();
	}
	return score;
}

double PartitioningScheme::getBicValue() {
	if (!isOptimized())
		return 0.0;
	if (!reoptimize_branch_lengths) {
		return getLinkedBicValue();
	}
	double score = 0.0;
	for (size_t i = 0; i < partitions.size(); i++) {
		PartitionElement * pe = partitions.at(i);
		score += pe->getBestModel()->getBicScore();
	}
	return score;
}

double PartitioningScheme::getAicValue() {
	if (!isOptimized())
		return 0.0;
	if (!reoptimize_branch_lengths) {
		return getLinkedAicValue();
	}
	double score = 0.0;
	for (size_t i = 0; i < partitions.size(); i++) {
		PartitionElement * pe = partitions.at(i);
		score += pe->getBestModel()->getAicScore();
	}
	return score;
}

double PartitioningScheme::getAiccValue() {
	if (!isOptimized())
		return 0.0;
	if (!reoptimize_branch_lengths) {
		return getLinkedAiccValue();
	}
	double score = 0.0;
	for (size_t i = 0; i < partitions.size(); i++) {
		PartitionElement * pe = partitions.at(i);
		score += pe->getBestModel()->getAiccScore();
	}
	return score;
}

double PartitioningScheme::getDTValue() {
	if (!isOptimized())
		return 0.0;
	double score = 0.0;
	for (size_t i = 0; i < partitions.size(); i++) {
		PartitionElement * pe = partitions.at(i);
		score += pe->getBestModel()->getDTScore();
	}
	return score;
}

double PartitioningScheme::getLinkedIcValue() {
	if (!isOptimized())
		return 0.0;
	return ModelSelector::computeIc(ic_type, getLnL(), (int) getNumberOfFreeParameters(),
			seq_len);
}

double PartitioningScheme::getLinkedBicValue() {
	if (!isOptimized())
		return 0.0;
	return ModelSelector::computeIc(BIC, getLnL(), (int) getNumberOfFreeParameters(),
			seq_len);
}

double PartitioningScheme::getLinkedAicValue() {
	if (!isOptimized())
		return 0.0;
	return ModelSelector::computeIc(AIC, getLnL(), (int) getNumberOfFreeParameters(),
			seq_len);
}

double PartitioningScheme::getLinkedAiccValue() {
	if (!isOptimized())
		return 0.0;
	return ModelSelector::computeIc(AICC, getLnL(), (int) getNumberOfFreeParameters(),
			seq_len);
}

double PartitioningScheme::getLinkedDTValue() {
	if (!isOptimized())
		return 0.0;
	return ModelSelector::computeIc(DT, getLnL(), (int) getNumberOfFreeParameters(),
			seq_len);
}

void PartitioningScheme::print(ostream & out) {
	out << getName() << endl;
	out << "Code Lines: " << codeLines << endl;
	for (int i = (int)codeLines-1; i >= 0; i--) {
		out << "  " << getCode(i) << endl;
	}
	out << setw(23) <<  "Num.Elements:" << numberOfElements << endl;
	if (isOptimized()) {
		out << setw(23) << "LnL:" << fixed << setprecision(4) << getLnL() << endl;
		out << setw(23) << "Num.Parameters:" << getNumberOfFreeParameters() << endl;
		out << setw(23) << "Selection score:" << fixed << setprecision(4) << getIcValue() << endl;
		out << setw(23) << "BIC score:" << fixed << setprecision(4) << getBicValue() << endl;
		out << setw(23) << "AIC score:" << fixed << setprecision(4) << getAicValue() << endl;
		out << setw(23) << "AICC score:" << fixed << setprecision(4) << getAiccValue() << endl;
		out << setw(23) << "Linked score:" << fixed << setprecision(4) << getLinkedIcValue() << endl;
		out << setw(23) << "Linked BIC_score:" << fixed << setprecision(4) << getLinkedBicValue() << endl;
		out << setw(23) << "Linked AIC_score:" << fixed << setprecision(4) << getLinkedAicValue() << endl;
		out << setw(23) << "Linked AICC_score:" << fixed << setprecision(4) << getLinkedAiccValue() << endl;
		if (tree > 0) {
			out << "Tree:         " << getTree() << endl;
		}
		out << endl;
		out << setw(6) << "--ID--"
			<< setw(11) << " --Model-- "
							<< setw(5) << " --K-"
							<< setw(15) << " -----LnL------"
							<< setw(15) << " -----BIC------"
							<< setw(15) << " -----AIC------"
							<< setw(15) << " -----AICc-----" << endl;
		for (size_t i=0; i<numberOfElements; i++) {
			out << setw(6) << left << i + 1
					<< setw(11) << left << getElement(i)->getBestModel()->getModel()->getName()
					<< setw(5) << right << getElement(i)->getBestModel()->getModel()->getNumberOfFreeParameters()
					<< setw(15) << right << fixed << setprecision(4) << getElement(i)->getBestModel()->getModel()->getLnL()
					<< setw(15) << right << fixed << setprecision(4) << getElement(i)->getBestModel()->getBicScore()
					<< setw(15) << right << fixed << setprecision(4) << getElement(i)->getBestModel()->getAicScore()
					<< setw(15) << right << fixed << setprecision(4) << getElement(i)->getBestModel()->getAiccScore() << endl;
		}
	}
}

} /* namespace partest */
