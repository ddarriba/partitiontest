/*
 * PartitionSelector.cpp
 *
 *  Created on: Mar 25, 2013
 *      Author: diego
 */

#include "PartitionSelector.h"
#include "indata/PartitionElement.h"
#include "ModelSelector.h"
#include "SelectionModel.h"
#include <alloca.h>
#include <iomanip>
#include <vector>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <time.h>

#define MAX_SCHEMES_SHOWN 10

namespace partest {

/** Functor for sorting the selection models */
struct compareSelectionPartitions {
	inline bool operator()(SelectionPartitioningScheme * struct1,
			SelectionPartitioningScheme * struct2) {
		return (struct1->value < struct2->value);
	}
};

PartitionSelector::PartitionSelector(PartitioningScheme ** schemesArray,
		int numberOfSchemes, ParTestOptions * options) :
		ic(options->getInformationCriterion()), sampleSize(
				options->getSampleSize()), sampleSizeValue(
				options->getSampleSizeValue()), numberOfSchemes(
				numberOfSchemes), schemesArray(schemesArray) {

	schemesVector = new vector<SelectionPartitioningScheme *>(numberOfSchemes);

// #pragma omp parallel for schedule(dynamic)
	for (int part_id = 0; part_id < numberOfSchemes; part_id++) {

		PartitioningScheme * scheme = schemesArray[part_id];

		double lk = 0;
		int parms = 0;
		double globalSampleSize = 0.0;
		double value = 0.0;
		for (int id = 0; id < scheme->getNumberOfElements(); id++) {
			PartitionElement * pe = scheme->getElement(id);
			SelectionModel * bestModel;
			if (!pe->getBestModel()) {
				double pSampleSize = 0.0;
				switch (sampleSize) {
				case SS_ALIGNMENT:
					pSampleSize = pe->getAlignment()->getNumSites();
					break;
				case SS_SHANNON:
					pSampleSize = pe->getAlignment()->getShannonEntropy();
					break;
				default:
					pSampleSize = sampleSize;
					break;
				}
				globalSampleSize += pSampleSize;

				ModelSelector selector = ModelSelector(pe, ic, pSampleSize);
				bestModel = selector.getBestModel()->clone();
				selector.print(*options->getModelsOutputStream());
			} else {
				bestModel = pe->getBestModel();
			}
#ifdef DEBUG
			cout << "[DEBUG_SEL]    " << pe->getName() << " : " << bestModel->getModel()->getLnL() << " -> " << bestModel->getValue() << endl;
#endif
			value += bestModel->getValue();
			lk += bestModel->getModel()->getLnL();

			/* Since we're optimizing branch lengths for each model, we need to count those parameters everytime */
			/* This means getNumberOfFreeParameters instead of getModelFreeParameters */
			parms += bestModel->getModel()->getNumberOfFreeParameters();
		}
		//double value = ModelSelector::computeIc(ic, lk, parms, globalSampleSize);
		(schemesVector->at(part_id)) = new SelectionPartitioningScheme();
		(schemesVector->at(part_id))->numParameters = parms;
		(schemesVector->at(part_id))->lnL = lk;
		(schemesVector->at(part_id))->scheme = scheme;
		(schemesVector->at(part_id))->value = value;
#ifdef DEBUG
		cout << "[DEBUG_SEL] BIC of " << scheme->getName() << " : " << value << endl;
#endif
	}

	std::sort(schemesVector->begin(), schemesVector->end(),
			compareSelectionPartitions());

	bestSelectionScheme = schemesVector->at(0);

	int limit = min(MAX_SCHEMES_SHOWN, numberOfSchemes);
	print(*options->getSchemesOutputStream(), limit);

	for (int part_id = 0; part_id < limit; part_id++) {
		schemesVector->at(part_id)->print(
				*options->getPartitionsOutputStream());
	}
}

PartitionSelector::~PartitionSelector() {
	for (unsigned int i = 0; i < schemesVector->size(); i++) {
		delete schemesVector->at(i);
	}
	delete schemesVector;
}

void PartitionSelector::print(ostream& out, int limit) {

	int maxSchemes = (limit == -1 ? schemesVector->size() : limit);
	int barLength = 75;

	out << endl;
	out << setw(barLength) << setfill('-') << "" << setfill(' ') << endl;
	out << setw(7) << "###" << setw(5) << " " << "Scheme" << endl;
	out << setw(barLength) << setfill('-') << "" << setfill(' ') << endl;
	for (int id = 0; id < maxSchemes; id++) {
		SelectionPartitioningScheme * sp = schemesVector->at(id);
		int numCodeLines = Utilities::iDecLog(sp->scheme->getNumberOfElements() - 1);
		out << setw(7) << id + 1 << setw(5) << " " << sp->scheme->getCode(numCodeLines) << endl;
		for (int i = (numCodeLines - 1) ; i >= 0; i--) {
			out << setw(12) << " " << sp->scheme->getCode(i) << endl;
		}
	}
	out << setw(barLength) << setfill('-') << "" << setfill(' ') << endl;
	out << setw(7) << "###" << setw(10) << "N" << setw(10)
			<< "K" << setw(20) << "lnL     " << setw(20) << "Value     " << endl;
	out << setw(barLength) << setfill('-') << "" << setfill(' ') << endl;
	for (int id = 0; id < maxSchemes; id++) {
		SelectionPartitioningScheme * sp = schemesVector->at(id);
		out << setw(7) << id + 1 << setw(10)
				<< sp->scheme->getNumberOfElements() << setw(10)
				<< fixed << setprecision(0) << sp->numParameters
				<< setw(20) << fixed << setprecision(4) << sp->lnL << setw(20)
				<< sp->value << endl;
	}
	if (limit > 0 && ((int) schemesVector->size()) > limit) {
		out << "Not showing " << schemesVector->size() - limit
				<< " schemes more..." << endl;
	}
	out << setw(barLength) << setfill('-') << "" << setfill(' ') << endl << endl;

}

void SelectionPartitioningScheme::print(ostream & out) {
	out << endl;
	int numCodeLines = Utilities::iDecLog(scheme->getNumberOfElements() - 1);
	out << "Scheme:     " << scheme->getCode(numCodeLines) << endl;
	for (int i = (numCodeLines - 1) ; i >= 0; i--) {
		out << "            " << scheme->getCode(i) << endl;
	}
	out << "Num.Params: " << numParameters << endl;
	out << "lnL:        " << lnL << endl;
	out << "value:      " << value << endl;
	out << setw(100) << setfill('-') << "" << setfill(' ') << endl;
	out << setw(5) << "###" << setw(15) << "Model" << setw(5) << "K" << setw(15)
			<< "lnL" << setw(15) << "Value" << setw(15) << "Delta" << setw(15)
			<< "Weight" << setw(15) << "CumWeight" << endl;
	out << setw(100) << setfill('-') << "" << setfill(' ') << endl;
	for (int i = 0; i < scheme->getNumberOfElements(); i++) {
		PartitionElement * pe = scheme->getElement(i);
		SelectionModel * selectionModel = pe->getBestModel();

		out << setw(5) << i + 1 << setw(15)
				<< selectionModel->getModel()->getName() << setw(5)
				<< selectionModel->getModel()->getNumberOfFreeParameters()
				<< setw(15) << setprecision(6)
				<< selectionModel->getModel()->getLnL() << setw(15) << fixed
				<< selectionModel->getValue() << setw(15)
				<< selectionModel->getDelta() << setw(15) << setprecision(4)
				<< selectionModel->getWeight() << setw(15)
				<< selectionModel->getCumWeight() << endl;
	}
	out << setw(100) << setfill('-') << "" << setfill(' ') << endl;
	for (int i = 0; i < scheme->getNumberOfElements(); i++) {
		out << setw(5) << i + 1 << ": " << scheme->getElement(i)->getName()
				<< endl;
	}
}
} /* namespace partest */
