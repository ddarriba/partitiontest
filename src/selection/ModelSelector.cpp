/*
 * ModelSelector.cpp
 *
 *  Created on: Jan 14, 2013
 *      Author: diego
 */

#include "ModelSelector.h"
#include <math.h>
#include <iomanip>
#include <stdlib.h>
#include <algorithm>

namespace partest {

/** Functor for sorting the selection models */
struct compareSelectionModels {
	inline bool operator()(SelectionModel * struct1, SelectionModel * struct2) {
		return (struct1->getValue() < struct2->getValue());
	}
};

double ModelSelector::computeIc(InformationCriterion ic, double lnL,
		int freeParameters, double sampleSize) {
	double value;
	switch (ic) {
	case AIC:
		value = 0.0;
		break;
	case AICC:
		value = 0.0;
		break;
	case BIC:
		value = computeBic(lnL, freeParameters, sampleSize);
		break;
	case DT:
		value = 0.0;
		break;
	}
	return value;
}
double ModelSelector::computeBic(double lnL, int freeParameters,
		double sampleSize) {
	return (-2 * lnL + freeParameters * log(sampleSize));
}

ModelSelector::ModelSelector(PartitionElement * partitionElement,
		InformationCriterion ic, double sampleSize) :
				ic(ic), sampleSize(sampleSize) {

	doSelection(partitionElement->getModelset(), ic, sampleSize);
	partitionElement->setBestModel(getBestModel());
}

ModelSelector::ModelSelector(ModelSet * modelset, InformationCriterion ic,
		double sampleSize) :
		ic(ic), sampleSize(sampleSize) {

	doSelection(modelset, ic, sampleSize);

	//print();
}

ModelSelector::~ModelSelector() {
	if (selectionModels) {
		for (int i = 0; i < selectionModels->size(); i++) {
			delete selectionModels->at(i);
		}
		delete selectionModels;
	}
}

double ModelSelector::getAlphaImportance() const {
	return alphaImportance;
}

double ModelSelector::getInvImportance() const {
	return invImportance;
}

double ModelSelector::getAlphaInvImportance() const {
	return alphaInvImportance;
}

double ModelSelector::getFImportance() const {
	return fImportance;
}

double ModelSelector::getOverallAlpha() const {
	return overallAlpha;
}

double ModelSelector::getOverallInv() const {
	return overallInv;
}

double ModelSelector::getOverallInvAlpha() const {
	return overallInvAlpha;
}

double ModelSelector::getOverallAlphaInv() const {
	return overallAlphaInv;
}

void ModelSelector::print() {

	/* header */
	cout << endl;
	cout << setw(95) << setfill('-') << "" << setfill(' ') << endl;
	cout << setw(5) << "###" << setw(15) << "Model" << setw(5) << "K" << setw(15) << "lnL"
			<< setw(15) << "Value" << setw(15) << "Delta" << setw(15)
			<< setprecision(4) << "Weight" << setw(15) << "CumWeight" << endl;
	cout << setw(95) << setfill('-') << "" << setfill(' ') << endl;
	for (int i = 0; i < selectionModels->size(); i++) {
		SelectionModel * selectionModel = selectionModels->at(i);

		cout << setw(5) << i + 1 << setw(15)
				<< selectionModel->getModel()->getName()
				<< setw(5) <<  selectionModel->getModel()->getNumberOfFreeParameters()
				<< setw(15)
				<< setprecision(10) << selectionModel->getModel()->getLnL()
				<< setw(15) << selectionModel->getValue() << setw(15)
				<< selectionModel->getDelta() << setw(15) << setprecision(4)
				<< selectionModel->getWeight() << setw(15)
				<< selectionModel->getCumWeight() << endl;
	}
	cout << setw(90) << setfill('-') << "" << setfill(' ') << endl << endl;
	cout << "SAMPLE SIZE: " << sampleSize << endl << endl;
	cout << "PARAMETER IMPORTANCE" << endl;
	cout << setw(15) << "alpha:" << alphaImportance << endl;
	cout << setw(15) << "pInv:" << invImportance << endl;
	cout << setw(15) << "alpha + pInv:" << alphaInvImportance << endl;
	cout << setw(15) << "frequencies:" << fImportance << endl;
	cout << endl << "OVERALLPARAMETER IMPORTANCE" << endl;
	cout << setw(15) << "alpha:" << overallAlpha << endl;
	cout << setw(15) << "pInv:" << overallInv << endl;
	cout << setw(15) << "alpha + pInv:" << overallInvAlpha << endl;
	cout << setw(15) << "pInv + alpha:" << overallAlphaInv << endl;
	cout << setw(90) << setfill('-') << "" << setfill(' ') << endl;
}

void ModelSelector::doSelection(ModelSet * modelset, InformationCriterion ic,
		double sampleSize) {
	selectionModels = new vector<SelectionModel *>(
			modelset->getNumberOfModels());

	int i;
	for (int i = 0; i < modelset->getNumberOfModels(); i++) {
		Model * model = modelset->getModel(i);
		double value = computeIc(ic, model->getLnL(),
				model->getNumberOfFreeParameters(), sampleSize);

		selectionModels->at(i) = new SelectionModel(model, value);
	}

	std::sort(selectionModels->begin(), selectionModels->end(),
			compareSelectionModels());

	bestModel = selectionModels->at(0);
	minValue = selectionModels->at(0)->getValue();
	double cumW = 0.0;
	double sumExp = 0.0;
	for (int i = 0; i < selectionModels->size(); i++) {
		SelectionModel * selectionModel = selectionModels->at(i);
		selectionModel->setIndex(i);
		selectionModel->setDelta(selectionModel->getValue() - minValue);
		sumExp += exp(-0.5 * selectionModel->getDelta());
	}

	alphaImportance = 0.0;
	invImportance = 0.0;
	alphaInvImportance = 0.0;
	fImportance = 0.0;
	overallAlpha = 0.0;
	overallInv = 0.0;
	overallAlphaInv = 0.0;
	overallInvAlpha = 0.0;

	for (int i = 0; i < selectionModels->size(); i++) {
		SelectionModel * selectionModel = selectionModels->at(i);
		Model * model = selectionModel->getModel();

		double weight =
				(selectionModel->getDelta() > 1000) ?
						0.0 : (exp(-0.5 * selectionModel->getDelta()) / sumExp);

		selectionModel->setWeight(weight);
		selectionModel->setCumWeight(cumW += selectionModel->getWeight());

		if (model->isGamma() & model->isPInv()) {
			alphaInvImportance += weight;
			overallAlphaInv += weight * model->getAlpha();
			overallInvAlpha += weight * model->getpInv();
		} else if (model->isGamma()) {
			alphaImportance += weight;
			overallAlpha += weight * model->getAlpha();
		} else if (model->isPInv()) {
			invImportance += weight;
			overallInv += weight * model->getpInv();
		}
		if (model->isPF()) {
			fImportance += weight;
		}
	}

	overallAlpha /= alphaImportance;
	overallInv /= invImportance;
	overallInvAlpha /= alphaInvImportance;
	overallAlphaInv /= alphaInvImportance;
}

} /* namespace partest */
