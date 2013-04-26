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
	inline bool operator()(SelectionPartition * struct1, SelectionPartition * struct2) {
		return (struct1->value < struct2->value);
	}
};

PartitionSelector::PartitionSelector(Partition ** partitions,
		int numberOfPartitions, InformationCriterion ic, SampleSize sampleSize, double sampleSizeValue) :
		ic(ic), sampleSize(sampleSize), sampleSizeValue(sampleSizeValue), numberOfPartitions(numberOfPartitions), partitions(
				partitions) {

	int part_id, id;
	partitionsVector = new vector<SelectionPartition *>(numberOfPartitions);
	for (part_id = 0; part_id < numberOfPartitions; part_id++) {

		Partition * partition = partitions[part_id];
//		ModelSelector * selectors = (ModelSelector *) alloca(
//				partition->getNumberOfElements() * sizeof(ModelSelector));
//		SelectionModel ** bestModels = (SelectionModel **) alloca(
//				partition->getNumberOfElements() * sizeof(SelectionModel *));

		double lk = 0;
		int parms = 0;
		double globalSampleSize = 0.0;
		double value = 0.0;
		for (id = 0; id < partition->getNumberOfElements(); id++) {
			PartitionElement * pe = partition->getElement(id);
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
		(partitionsVector->at(part_id)) = new SelectionPartition();
		(partitionsVector->at(part_id))->numParameters = parms;
		(partitionsVector->at(part_id))->lnL = lk;
		(partitionsVector->at(part_id))->partition = partition;
		(partitionsVector->at(part_id))->value = value;
	}
	std::sort(partitionsVector->begin(), partitionsVector->end(),
				compareSelectionPartitions());

	bestPartition = partitionsVector->at(0)->partition;
	print();

}

PartitionSelector::~PartitionSelector() {
	delete partitionsVector;
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
	for (id = 0; id < partitionsVector->size(); id++) {
		SelectionPartition * sp = partitionsVector->at(id);
		cout << sp->partition->toString() << " " << sp->lnL << " "
				<< sp->numParameters << " " << sp->value << endl;
	}
}

} /* namespace partest */
