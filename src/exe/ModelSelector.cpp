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

#include "ModelSelector.h"
#include <math.h>
#include <fstream>
#include <iostream>
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

ModelSelector::ModelSelector(PartitionElement * partitionElement,
		InformationCriterion ic, double sampleSize) :
		partitionElement(partitionElement), ic(ic), sampleSize(sampleSize) {

	doSelection(partitionElement->getModels(), ic, sampleSize);
	partitionElement->setBestModel(getBestModel());

	if (outputAvailable && models_logfile) {
		ofstream ofs(models_logfile->c_str(), ios::in | ios::out | ios::app);
		print(ofs);
		ofs.close();
	}
}

ModelSelector::~ModelSelector() {
	if (selectionModels) {
		for (size_t i = 0; i < selectionModels->size(); i++) {
			delete selectionModels->at(i);
		}
		delete selectionModels;
	}
}

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
		cerr << "[ERROR] Decision Theory is not implemented yet" << endl;
		exit_partest(EX_UNAVAILABLE);
		break;
	default:
		value = 0.0;
		cerr << "[ERROR] Undefined Criterion" << endl;
		exit_partest(EX_UNAVAILABLE);
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

double ModelSelector::getAlphaImportance() const {
	return alphaImportance;
}

double ModelSelector::getFImportance() const {
	return fImportance;
}

double ModelSelector::getOverallAlpha() const {
	return overallAlpha;
}

#ifdef _IG_MODELS
double ModelSelector::getInvImportance() const {
	return invImportance;
}

double ModelSelector::getAlphaInvImportance() const {
	return alphaInvImportance;
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
#endif

void ModelSelector::doSelection(vector<Model *> modelset,
		InformationCriterion ic, double sampleSize) {
	selectionModels = new vector<SelectionModel *>(modelset.size());

	for (size_t i = 0; i < modelset.size(); i++) {
		Model * model = modelset.at(i);
		int freeParameters = reoptimize_branch_lengths
				?model->getNumberOfFreeParameters()
				:model->getModelFreeParameters();
		double value = computeIc(ic, model->getLnL(),
				freeParameters, sampleSize);
		double bicScore = computeIc(BIC, model->getLnL(),
				freeParameters, sampleSize);
		double aicScore = computeIc(AIC, model->getLnL(),
				freeParameters, sampleSize);
		double aiccScore = computeIc(AICC, model->getLnL(),
				freeParameters, sampleSize);

#ifdef DEBUG
		cout << "[DEBUG_SEL] " << model->getName() << " " << model->getLnL() << " " << model->getNumberOfFreeParameters() << " "
		<< sampleSize << " " << value << endl;
#endif

		selectionModels->at(i) = new SelectionModel(model, value);
		selectionModels->at(i)->setBicScore(bicScore);
		selectionModels->at(i)->setAicScore(aicScore);
		selectionModels->at(i)->setAiccScore(aiccScore);
	}

	std::sort(selectionModels->begin(), selectionModels->end(),
			compareSelectionModels());

	bestModel = selectionModels->at(0);
	minValue = selectionModels->at(0)->getValue();
	double cumW = 0.0;
	double sumExp = 0.0;
	for (size_t i = 0; i < selectionModels->size(); i++) {
		SelectionModel * selectionModel = selectionModels->at(i);
		selectionModel->setIndex(i);
		selectionModel->setDelta(selectionModel->getValue() - minValue);
		sumExp += exp(-0.5 * selectionModel->getDelta());
	}

	alphaImportance = 0.0;
	fImportance = 0.0;
	overallAlpha = 0.0;
#ifdef _IG_MODELS
	invImportance = 0.0;
	alphaInvImportance = 0.0;
	overallInv = 0.0;
	overallAlphaInv = 0.0;
	overallInvAlpha = 0.0;
#endif

	for (size_t i = 0; i < selectionModels->size(); i++) {
		SelectionModel * selectionModel = selectionModels->at(i);
		Model * model = selectionModel->getModel();

		double weight =
				(selectionModel->getDelta() > 1000) ?
						0.0 : (exp(-0.5 * selectionModel->getDelta()) / sumExp);

		selectionModel->setWeight(weight);
		selectionModel->setCumWeight(cumW += selectionModel->getWeight());

#ifdef _IG_MODELS
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
#else
		if (model->isGamma()) {
			alphaImportance += weight;
			overallAlpha += weight * model->getAlpha();
		}
#endif
		if (model->isPF()) {
			fImportance += weight;
		}
	}

	if (alphaImportance > 0.0)
		overallAlpha /= alphaImportance;
#ifdef _IG_MODELS
	if (invImportance > 0.0)
	overallInv /= invImportance;
	if (alphaInvImportance > 0.0)
	overallInvAlpha /= alphaInvImportance;
	if (alphaInvImportance > 0.0)
	overallAlphaInv /= alphaInvImportance;
#endif
}

void ModelSelector::print(ostream& out) {

	/* header */
	out << endl;
	out << "  Partition: " << partitionElement->getName() << endl;
	//out << "         Id: " << partitionElement->getId() << endl;
	//out << "  Num.genes: " << partitionElement->getNumberOfSections() << endl;
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
	default:
		cerr << "[ERROR] Undefined Criterion" << endl;
		exit_partest(EX_SOFTWARE);
		break;
	}
	out << "Sample size: " << sampleSize << endl << endl;
	out << setw(115) << setfill('-') << "" << setfill(' ') << endl;
	out << setw(5) << "###" << setw(15) << "Model" << setw(10) << "K"
			<< setw(20) << "lnL" << setw(20) << "Value" << setw(15) << "Delta"
			<< setw(15) << setprecision(4) << "Weight" << setw(15)
			<< "CumWeight" << endl;
	out << setw(115 ) << setfill('-') << "" << setfill(' ') << endl;
	for (size_t i = 0; i < selectionModels->size(); i++) {
		SelectionModel * selectionModel = selectionModels->at(i);

		out << setw(5) << i + 1 << setw(15)
				<< selectionModel->getModel()->getName() << setw(10)
				<< selectionModel->getModel()->getNumberOfFreeParameters()
				<< setw(20) << fixed << setprecision(4)
				<< selectionModel->getModel()->getLnL() << setw(20)
				<< selectionModel->getValue() << setw(15)
				<< selectionModel->getDelta() << setw(15) << setprecision(4)
				<< selectionModel->getWeight() << setw(15)
				<< selectionModel->getCumWeight() << endl;
	}
	out << setw(115) << setfill('-') << "" << setfill(' ') << endl;

	out << setw(5) << "###" << setw(20) << "BIC" << setw(20) << "AIC"
			<< setw(20) << "AICc" << endl;
	out << setw(66 ) << setfill('-') << "" << setfill(' ') << endl;
	for (size_t i = 0; i < selectionModels->size(); i++) {
		SelectionModel * selectionModel = selectionModels->at(i);

		out << setw(5) << i + 1
				<< setw(20) << fixed << setprecision(4) << selectionModel->getBicScore()
				<< setw(20) << fixed << setprecision(4) << selectionModel->getAicScore()
				<< setw(20) << fixed << setprecision(4) << selectionModel->getAiccScore() << endl;
	}
	out << setw(66) << setfill('-') << "" << setfill(' ') << endl << endl;

	if (selectionModels->size() > 1) {
		out << "PARAMETER IMPORTANCE" << endl;
		out << setw(16) << "alpha: " << alphaImportance << endl;
#ifdef _IG_MODELS
		out << setw(16) << "pInv: " << invImportance << endl;
		out << setw(16) << "alpha + pInv: " << alphaInvImportance << endl;
#endif
		out << setw(16) << "frequencies: " << fImportance << endl;
		out << endl << "OVERALLPARAMETER IMPORTANCE" << endl;
		out << setw(16) << "alpha:" << overallAlpha << endl;
#ifdef _IG_MODELS
		out << setw(16) << "pInv:" << overallInv << endl;
		out << setw(16) << "alpha + pInv:" << overallInvAlpha << endl;
		out << setw(16) << "pInv + alpha:" << overallAlphaInv << endl;
#endif
		out << setw(115) << setfill('-') << "" << setfill(' ') << endl << endl;
	}

	out << setw(115) << setfill('-') << "" << setfill(' ') << endl;
	out << setw(5) << "###";
	if (data_type == DT_NUCLEIC) {
		out << setw(10) << "f(A)" << setw(10) << "f(C)" << setw(10) << "f(G)"
				<< setw(10) << "f(T)" << setw(10) << "R(a->c)" << setw(10)
				<< "R(a->g)" << setw(10) << "R(a->t)" << setw(10) << "R(c->g)"
				<< setw(10) << "R(c->t)" << setw(10) << "R(g->t)";
	}
	out << setw(10) << "alpha";
#ifdef _IG_MODELS
	out << setw(10) << "pInv";
#endif
	out << endl;
	out << setw(115) << setfill('-') << "" << setfill(' ') << endl;
	for (size_t i = 0; i < selectionModels->size(); i++) {
		Model * model = selectionModels->at(i)->getModel();

		out << setw(5) << i + 1;
		if (data_type == DT_NUCLEIC) {
			for (int j = 0; j < NUM_NUC_FREQS; j++) {
				out << setw(10) << setprecision(4)
						<< model->getFrequencies()[j];
			}
			for (int j = 0; j < NUM_DNA_RATES; j++) {
				out << setw(10) << setprecision(4)
						<< min(model->getRates()[j], 999.9999);
			}
		}
		out << setw(10) << fixed << setprecision(4)
				<< min(model->getAlpha(), 1000.0000);
#ifdef _IG_MODELS
		out << setw(10) << fixed << setprecision(4) << model->getpInv();
#endif
		out << endl;
	}
	out << setw(122) << setfill('-') << "" << setfill(' ') << endl;
	// bestModel->getModel()->print(out);
	out << endl << "Best Tree: " << bestModel->getModel()->getTree() << endl;
}

} /* namespace partest */
