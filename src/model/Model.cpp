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
 *  For any other enquiries send an Email to Diego Darriba
 *  ddarriba@udc.es
 */

#include "Model.h"
#include "../util/Utilities.h"
#include <assert.h>
#include <stdlib.h>
#include <iostream>

namespace partest {

Model::Model(bitMask rateVariation, int numberOfTaxa) :
		rateVariation(rateVariation), alpha(100.0), pInv(0.0), lnL(0.0),
		numberOfFrequencies(0), frequencies(0), rates(0), modelFreeParameters(0), name() {

	/* The free parameters are initialized with the number of branches */
	treeFreeParameters = Utilities::numberOfBranches(numberOfTaxa);
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
string Model::getMatrixName() {
	return matrixName;
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

void Model::setpInv(double pInv) {
	assert(isPInv());
	this->pInv = pInv;
}

bool Model::isPInv() {
	return !(~rateVariation & RateVarI);
}

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
		if (isPInv()) {
			cout << prefix << "pInv:  " << pInv << endl;
		}
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

void Model::setFrequencies(const double * frequencies) {
// ignore
}

double * Model::getRates(void) const {
	return rates;
}

void Model::setRates(const double * rates) {
// ignore
}

} /* namespace partest */
