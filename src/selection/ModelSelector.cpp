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
		value = computeAic(lnL, freeParameters);
		break;
	case AICC:
		value = computeAicc(lnL, freeParameters, sampleSize);
		break;
	case BIC:
		value = computeBic(lnL, freeParameters, sampleSize);
		break;
	case DT:
		value = 0.0;
		cerr << "ERROR: Decision Theory is not implemented yet" << endl;
		Utilities::exit_partest(EX_UNAVAILABLE);
		break;
	default:
		value = 0.0;
		cerr << "ERROR: Undefined Criterion" << endl;
		Utilities::exit_partest(EX_UNAVAILABLE);
		break;
	}
	return value;
}

double ModelSelector::computeBic(double lnL, int freeParameters,
		double sampleSize) {
	return (-2 * lnL + freeParameters * log(sampleSize));
}

double ModelSelector::computeAic(double lnL, int freeParameters) {
	return (2 * (freeParameters - lnL));
}

double ModelSelector::computeAicc(double lnL, int freeParameters,
		double sampleSize) {
	double aicc = computeAic(lnL, freeParameters);
	aicc += (2 * freeParameters * (freeParameters + 1))
			/ (sampleSize - freeParameters - 1);
	return aicc;
}

ModelSelector::ModelSelector(PartitionElement * partitionElement,
		InformationCriterion ic, double sampleSize) :
		ic(ic), sampleSize(sampleSize), partitionElement(partitionElement) {

	dataType = partitionElement->getModelset()->getDataType();
	doSelection(partitionElement->getModelset(), ic, sampleSize);
	partitionElement->setBestModel(getBestModel());
}

ModelSelector::ModelSelector(ModelSet * modelset, InformationCriterion ic,
		double sampleSize) :
		ic(ic), sampleSize(sampleSize) {

	dataType = modelset->getDataType();
	doSelection(modelset, ic, sampleSize);

	//print();
}

ModelSelector::~ModelSelector() {
	if (selectionModels) {
		for (unsigned int i = 0; i < selectionModels->size(); i++) {
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

void ModelSelector::print(ostream& out) {

	/* header */
	out << endl;
	out << "  Partition: " << partitionElement->getName() << endl;
	//out << "         Id: " << partitionElement->getId() << endl;
	out << "  Num.genes: " << partitionElement->getNumberOfSections() << endl;
	out << "  Criterion: ";
	switch (ic) {
	case AIC:
		out << "Akaike Information Criterion" << endl;
		break;
	case AICC:
		out << "Corrected Akaike Information Criterion" << endl;
		break;
	case BIC:
		out << "Bayesian Information Criterion" << endl;
		break;
	case DT:
		out << "Decision Theory" << endl;
		break;
	}
	out << "Sample size: " << sampleSize << endl << endl;
	out << setw(110) << setfill('-') << "" << setfill(' ') << endl;
	out << setw(5) << "###" << setw(15) << "Model" << setw(10) << "K" << setw(20)
			<< "lnL" << setw(20) << "Value" << setw(15) << "Delta" << setw(15)
			<< setprecision(4) << "Weight" << setw(15) << "CumWeight" << endl;
	out << setw(110) << setfill('-') << "" << setfill(' ') << endl;
	for (unsigned int i = 0; i < selectionModels->size(); i++) {
		SelectionModel * selectionModel = selectionModels->at(i);

		out << setw(5) << i + 1 << setw(15)
				<< selectionModel->getModel()->getName() << setw(10)
				<< selectionModel->getModel()->getNumberOfFreeParameters()
				<< setw(20) << setprecision(4)
				<< selectionModel->getModel()->getLnL() << setw(20)
				<< selectionModel->getValue() << setw(15)
				<< selectionModel->getDelta() << setw(15) << setprecision(4)
				<< selectionModel->getWeight() << setw(15)
				<< selectionModel->getCumWeight() << endl;
	}
	out << setw(110) << setfill('-') << "" << setfill(' ') << endl << endl;

	if (selectionModels->size() > 1) {
		out << "PARAMETER IMPORTANCE" << endl;
		out << setw(16) << "alpha: " << alphaImportance << endl;
		out << setw(16) << "pInv: " << invImportance << endl;
		out << setw(16) << "alpha + pInv: " << alphaInvImportance << endl;
		out << setw(16) << "frequencies: " << fImportance << endl;
		out << endl << "OVERALLPARAMETER IMPORTANCE" << endl;
		out << setw(16) << "alpha:" << overallAlpha << endl;
		out << setw(16) << "pInv:" << overallInv << endl;
		out << setw(16) << "alpha + pInv:" << overallInvAlpha << endl;
		out << setw(16) << "pInv + alpha:" << overallAlphaInv << endl;
		out << setw(100) << setfill('-') << "" << setfill(' ') << endl << endl;
	}


	out << setw(122) << setfill('-') << "" << setfill(' ') << endl;
	out << setw(5) << "###";
	if (dataType == DT_NUCLEIC) {
		out << setw(10) << "f(A)" << setw(10) << "f(C)" << setw(10) << "f(G)"
				<< setw(10) << "f(T)" << setw(10) << "R(a->c)" << setw(10)
				<< "R(a->g)" << setw(10) << "R(a->t)" << setw(10) << "R(c->g)"
				<< setw(10) << "R(c->t)" << setw(10) << "R(g->t)";
	}
	out << setw(10) << "alpha" << setw(10) << "pInv" << endl;
	out << setw(122) << setfill('-') << "" << setfill(' ') << endl;
	for (unsigned int i = 0; i < selectionModels->size(); i++) {
		Model * model = selectionModels->at(i)->getModel();

		out << setw(5) << i + 1;
		if (dataType == DT_NUCLEIC) {
			for (int j = 0; j < NUM_NUC_FREQS; j++) {
				out << setw(10) << setprecision(4)
						<< model->getFrequencies()[j];
			}
			for (int j = 0; j < NUM_RATES; j++) {
				out << setw(10) << setprecision(4) << model->getRates()[j];
			}
		}
		out << setw(10) << setprecision(4) << model->getAlpha() << setw(10)
				<< setprecision(4) << model->getpInv();
		out << endl;
	}
	out << setw(122) << setfill('-') << "" << setfill(' ') << endl;
	// bestModel->getModel()->print(out);
	out << endl << "Tree for " << partitionElement->getName() << ": " << bestModel->getModel()->getTree() << endl;
}

void ModelSelector::doSelection(ModelSet * modelset, InformationCriterion ic,
		double sampleSize) {
	selectionModels = new vector<SelectionModel *>(
			modelset->getNumberOfModels());

	for (unsigned int i = 0; i < modelset->getNumberOfModels(); i++) {
		Model * model = modelset->getModel(i);
		double value = computeIc(ic, model->getLnL(),
				model->getNumberOfFreeParameters(), sampleSize);
#ifdef DEBUG
		cout << "[DEBUG_SEL] " << model->getName() << " " << model->getLnL() << " " << model->getNumberOfFreeParameters() << " "
				<< sampleSize << " " << value << endl;
#endif

		selectionModels->at(i) = new SelectionModel(model, value);
	}

	std::sort(selectionModels->begin(), selectionModels->end(),
			compareSelectionModels());

	bestModel = selectionModels->at(0);
	minValue = selectionModels->at(0)->getValue();
	double cumW = 0.0;
	double sumExp = 0.0;
	for (unsigned int i = 0; i < selectionModels->size(); i++) {
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

	for (unsigned int i = 0; i < selectionModels->size(); i++) {
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

	if (alphaImportance > 0.0)
		overallAlpha /= alphaImportance;
	if (invImportance > 0.0)
		overallInv /= invImportance;
	if (alphaInvImportance > 0.0)
		overallInvAlpha /= alphaInvImportance;
	if (alphaInvImportance > 0.0)
		overallAlphaInv /= alphaInvImportance;
}

} /* namespace partest */
