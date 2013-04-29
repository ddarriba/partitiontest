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
#include <vector>
#include <iomanip>
#include <sstream>
#include <algorithm>

namespace partest {

/** Functor for sorting the selection models */
struct compareSelectionPartitions {
	inline bool operator()(SelectionPartitioningScheme * struct1, SelectionPartitioningScheme * struct2) {
		return (struct1->value < struct2->value);
	}
};

PartitionSelector::PartitionSelector(PartitioningScheme ** schemesArray,
		int numberOfSchemes, InformationCriterion ic, SampleSize sampleSize, double sampleSizeValue) :
		ic(ic), sampleSize(sampleSize), sampleSizeValue(sampleSizeValue), numberOfSchemes(numberOfSchemes), schemesArray(
				schemesArray) {

	int part_id, id;
	schemesVector = new vector<SelectionPartitioningScheme *>(numberOfSchemes);
	for (part_id = 0; part_id < numberOfSchemes; part_id++) {

		PartitioningScheme * scheme = schemesArray[part_id];
//		ModelSelector * selectors = (ModelSelector *) alloca(
//				partition->getNumberOfElements() * sizeof(ModelSelector));
//		SelectionModel ** bestModels = (SelectionModel **) alloca(
//				partition->getNumberOfElements() * sizeof(SelectionModel *));

		double lk = 0;
		int parms = 0;
		double globalSampleSize = 0.0;
		double value = 0.0;
		for (id = 0; id < scheme->getNumberOfElements(); id++) {
			PartitionElement * pe = scheme->getElement(id);
			double pSampleSize = sampleSize;
			switch (sampleSize) {
			case SS_ALIGNMENT:
				pSampleSize = pe->getAlignment()->getNumSites() * pe->getAlignment()->getNumSeqs();
				break;
			case SS_SHANNON:
				pSampleSize = pe->getAlignment()->getShannonEntropy();
				break;
			}
			globalSampleSize += pSampleSize;

			ModelSelector selector = ModelSelector(pe, ic,
					pSampleSize);
			SelectionModel * bestModel = selector.getBestModel();
			value += selector.getBestModel()->getValue();
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
	}
	std::sort(schemesVector->begin(), schemesVector->end(),
				compareSelectionPartitions());

	bestScheme = schemesVector->at(0)->scheme;
	print();

	for (part_id = 0; part_id < numberOfSchemes; part_id++) {
		delete schemesVector->at(part_id);
	}

}

PartitionSelector::~PartitionSelector() {
	delete schemesVector;
}

void PartitionSelector::print() {
	cout << "*********** PARTITIONS SELECTION *************" << endl;
	cout << " Sample size: ";
	switch(sampleSize) {
	case SS_ALIGNMENT:
	cout << "Alignment size "<< endl;
	break;
	case SS_SHANNON:
	cout << "Shannon entropy "<< endl;
	break;
	case SS_CUSTOM:
		cout << "Custom value "<< sampleSizeValue << endl;
		break;
	}

	int id;
	for (id = 0; id < schemesVector->size(); id++) {
		SelectionPartitioningScheme * sp = schemesVector->at(id);
		cout << sp->scheme->toString() << " " << sp->lnL << " "
				<< sp->numParameters << " " << sp->value << endl;
	}
}

} /* namespace partest */
