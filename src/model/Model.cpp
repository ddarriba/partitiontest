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

#include "Model.h"
#include "util/Utilities.h"
#include <assert.h>
#include <stdlib.h>
#include <iostream>

namespace partest {

Model::Model(bitMask rateVariation, int numberOfTaxa) :
		rateVariation(rateVariation), lnL(0.0), alpha(1.0),
		rates(0), frequencies(0), numberOfFrequencies(0),
		name(), modelFreeParameters(0) {

#ifdef _IG_MODELS
	pInv = 0.0;
#endif

	/* The free parameters are initialized with the number of branches */
	treeFreeParameters = reoptimize_branch_lengths?Utilities::numberOfBranches(numberOfTaxa):0;
}

Model::~Model() {
	if (frequencies)
		free(frequencies);
	if (rates)
		free(rates);
}

string Model::getName() {
	return name;
}

void Model::setName(string name) {
	this->name = name;
}

string Model::getMatrixName() {
	return matrixName;
}

void Model::setMatrixName(string matrixName) {
	this->matrixName = matrixName;
}

bitMask Model::getRateVariation() {
	return rateVariation;
}
string Model::getTree() {
	return tree;
}

void Model::setAlpha(double alpha) {
	assert(isGamma());
	this->alpha = alpha;
}

#ifdef _IG_MODELS
void Model::setpInv(double pInv) {
	assert(isPInv());
	this->pInv = pInv;
}
#endif

#ifdef _IG_MODELS
bool Model::isPInv() {
	return !(~rateVariation & RateVarI);
}
#endif

bool Model::isGamma() {
	return !(~rateVariation & RateVarG);
}

bool Model::isPF() {
	return !(~rateVariation & RateVarF);
}

bool Model::isOptimized() {
	return (lnL < 0.0);
}

void Model::print(ostream& cout, const char * prefix) {
	cout << prefix << "Name:  " << name << endl;
	if (isOptimized()) {
		cout << prefix << "lnL:   " << lnL << endl;
#ifdef _IG_MODELS
		if (isPInv()) {
			cout << prefix << "pInv:  " << pInv << endl;
		}
#endif
		if (isGamma()) {
			cout << prefix << "alpha: " << alpha << endl;
		}
		cout << prefix << "Most Likely Tree: " << tree;
	} else {
		cout << prefix << "State: Unoptimized" << endl;
	}
}

void Model::setTree(string tree) {
	this->tree = tree;
}

void Model::setTree(char * tree) {
	this->tree = string(tree);
}

double * Model::getFrequencies(void) const {
	return frequencies;
}

double * Model::getRates(void) const {
	return rates;
}

} /* namespace partest */
