/*
 * Model.cpp
 *
 *  Created on: Jan 8, 2013
 *      Author: diego
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

void Model::print() {
	cout << "--------" << endl;
	cout << name << endl;
	if (isOptimized()) {
		cout << "lnL = " << lnL << endl;
		if (isPInv()) {
			cout << "pInv = " << pInv << endl;
		}
		if (isGamma()) {
			cout << "alpha = " << alpha << endl;
		}
		cout << "MOST LIKELY TREE:" << endl;
		cout << tree << endl;
	} else {
		cout << "Unoptimized" << endl;
	}
	cout << "--------" << endl;
}

void Model::setTree(string tree) {
	this->tree = tree;
}

void Model::setTree(char * tree) {
	this->tree = string(tree);
}

void Model::setFrequencies(const double * frequencies) {
// ignore
}

void Model::setRates(const double * rates) {
// ignore
}

} /* namespace partest */
