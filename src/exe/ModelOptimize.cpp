/*
 * ModelOptimize.cpp
 *
 *  Created on: Apr 9, 2014
 *      Author: diego
 */

#include "ModelOptimize.h"

#include "exe/ModelSelector.h"
#include "util/Utilities.h"

#include <iostream>
#include <iomanip>
#include <pll.h>
#include <math.h>
#include <string.h>
#include <assert.h>

using namespace std;

namespace partest {

ModelOptimize::ModelOptimize() {

}

ModelOptimize::~ModelOptimize() {

}

int ModelOptimize::optimizePartitioningScheme(PartitioningScheme * scheme, int index, int limit) {

	cout << " - scheme " << setw(Utilities::iDecLog(limit)+1) << setfill('0') << right << index+1 << "/" << limit  << setfill(' ') << endl;

	for (int cur_element = 0; cur_element < scheme->getNumberOfElements();
			cur_element++) {
		PartitionElement * element = scheme->getElement(cur_element);
		optimizePartitionElement(element, cur_element, scheme->getNumberOfElements());
	}

	return EX_OK;
}

int ModelOptimize::optimizePartitionElement(PartitionElement * element, int index, int limit) {

	cout << " - - element " << setw(Utilities::iDecLog(limit)+1) << setfill('0') << right << index+1 << "/" << limit  << setfill(' ') << endl;

	if (element->isOptimized()) return EX_OK;

	element->setupStructures();

	for (int modelIndex = 0; modelIndex < element->getNumberOfModels();
			modelIndex++) {
		optimizeModel(element, modelIndex, element->getNumberOfModels());
	}

	ModelSelector ms(element, BIC, element->getSampleSize());
	//ms.print(cout);

	element->destroyStructures();

	return EX_OK;
}

void ModelOptimize::setModelParameters(Model * _model, pllInstance * _tree,
		partitionList * _partitions, pllAlignmentData * _alignData, int index,
		bool setAlphaFreqs) {
	pInfo * current_part = _partitions->partitionData[index];

	if (data_type == DT_NUCLEIC) {
		current_part->optimizeBaseFrequencies = _model->isPF();
		current_part->alpha = _model->getAlpha();
		current_part->nonGTR = PLL_FALSE;
		current_part->dataType = PLL_DNA_DATA;
		const char * m = _model->getMatrixName().c_str();
		char * symmetryPar = (char *) malloc(12 * sizeof(char));
		symmetryPar[0] = m[0];
		symmetryPar[11] = '\0';
		for (int j = 1; j < 6; j++) {
			symmetryPar[(j - 1) * 2 + 1] = ',';
			symmetryPar[j * 2] = m[j];
		}

		pllSetSubstitutionRateMatrixSymmetries(symmetryPar, _partitions, index);

		if (setAlphaFreqs) {
			memcpy(current_part->frequencies, _model->getFrequencies(),
					4 * sizeof(double));
			memcpy(current_part->substRates, _model->getRates(),
					6 * sizeof(double));
			current_part->alpha = _model->getAlpha();
		} else {
			current_part->optimizeBaseFrequencies = _model->isPF();
			if (!_model->isPF()) {
				for (int i = 0; i < 4; i++) {
					current_part->frequencies[i] = 0.25;
				}
			}
			for (int i = 0; i < 6; i++) {
				current_part->substRates[i] = 1;
			}

			current_part->alpha = 100;
		}

		free(symmetryPar);
	} else {
		ProteicModel * pModel = static_cast<ProteicModel *>(_model);
		current_part->dataType = PLL_AA_DATA;
		current_part->protFreqs = pModel->isPF();
		current_part->optimizeBaseFrequencies = PLL_FALSE;
		current_part->protModels = pModel->getMatrix();
		current_part->alpha = pModel->getAlpha();
	}
	//TODO: This works if partitions has one single partitions
	assert(_partitions->numberOfPartitions == 1);
	double **ef = pllBaseFrequenciesGTR(_partitions, _alignData);
	initModel(_tree, ef, _partitions);
	free(*ef); free(ef);
}

void ModelOptimize::optimizeModel(PartitionElement * element, unsigned int modelIndex, int limit) {
	pllInstance * _tree = element->getTree();
	partitionList * _partitions = element->getPartitions();
	pllAlignmentData * _alignData = element->getAlignData();
	Model * model = element->getModel(modelIndex);

	/* set parameters for single partition element */
	setModelParameters(model, _tree, _partitions, _alignData, 0, false);
	double lk;
	double epsilon = 0.1;

	_tree->thoroughInsertion = PLL_FALSE;

	pllEvaluateLikelihood(_tree, _partitions, _tree->start, PLL_TRUE,
			PLL_FALSE);

	do {
		lk = _tree->likelihood;
		pllOptimizeModelParameters(_tree, _partitions, 1);
		pllOptimizeBranchLengths(_tree, _partitions, 10);
		pllEvaluateLikelihood(_tree, _partitions, _tree->start, PLL_TRUE,
				PLL_FALSE);
	} while (fabs(lk - _tree->likelihood) > epsilon);

	pllTreeToNewick(_tree->tree_string, _tree, _partitions, _tree->start->back,
			PLL_TRUE, PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE,
			PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);
	_tree->tree_string[_tree->treeStringLength-10] = '\0';
	model->setLnL(_tree->likelihood);
	model->setTree(_tree->tree_string);

	model->setFrequencies(_partitions->partitionData[0]->frequencies);
	if (model->isGamma())
		model->setAlpha(_partitions->partitionData[0]->alpha);
	model->setRates(_partitions->partitionData[0]->substRates);

	cout << " - - - " << setw(Utilities::iDecLog(limit)+1) << setfill('0') << right << modelIndex+1 << "/" << limit  << " " << model->getName() << " (" << _tree->likelihood << ")" << setfill(' ') << endl;
}

} /* namespace partest */
