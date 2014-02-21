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

/**
 * @file ModelOptimize.cpp
 *
 * @brief Implementation for common model optimizers interface
 */

#include "ModelOptimize.h"
#include "selection/ModelSelector.h"

namespace partest {

ModelOptimize::ModelOptimize(ParTestOptions * options) :
		options(options) {

}

ModelOptimize::~ModelOptimize() {
}

int ModelOptimize::optimizePartitionElement(PartitionElement * partitionElement,
		int current_index, int max_index) {
	ModelSet * modelset = partitionElement->getModelset();
	//(options->getRateVariation(), options->getDataType());
#pragma omp single
	{
		notify_observers(MT_MODELSET_INIT, partitionElement->getId(), modelset,
				time(NULL), current_index, max_index,
				partitionElement->getName());
	}
#pragma omp for schedule(dynamic)
	for (int i = 0; i < modelset->getNumberOfModels(); i++) {
		//for (int i = modelset->getNumberOfModels() - 1; i >= 0; i--) {
		notify_observers(MT_SINGLE_INIT, partitionElement->getId(),
				modelset->getModel(i), time(NULL), i + 1,
				modelset->getNumberOfModels(),
				modelset->getModel(i)->getName());

		optimizeModel(modelset->getModel(i), partitionElement, i,
				modelset->getNumberOfModels());

		notify_observers(MT_SINGLE_END, partitionElement->getId(),
				modelset->getModel(i), time(NULL), i + 1,
				modelset->getNumberOfModels(),
				modelset->getModel(i)->getName());
	}
#pragma omp single
	{
		double pSampleSize = 0.0;
		switch (options->getSampleSize()) {
		case SS_ALIGNMENT:
			pSampleSize = partitionElement->getAlignment()->getNumSites()
					* partitionElement->getAlignment()->getNumSeqs();
			break;
		case SS_SHANNON:
			pSampleSize = partitionElement->getAlignment()->getShannonEntropy();
			break;
		default:
			pSampleSize = options->getSampleSize();
			break;
		}
		ModelSelector selector = ModelSelector(partitionElement,
				options->getInformationCriterion(), pSampleSize);

		notify_observers(MT_MODELSET_END, partitionElement->getId(),
				partitionElement->getBestModel()->getModel(), time(NULL),
				current_index, max_index, partitionElement->getName());
	}

	return 0;
}

int ModelOptimize::optimizePartitioningScheme(PartitioningScheme * scheme,
		bool forceRecomputation, int current_index, int max_index) {
	if (!scheme->isOptimized() || forceRecomputation) {
		string pString(scheme->getName());

		t_partitionElementId nullId;
		notify_observers(MT_SCHEME_INIT, nullId, time(NULL), current_index,
				max_index, pString);
		for (int i = 0; i < scheme->getNumberOfElements(); i++) {
			PartitionElement * currentElement = scheme->getElement(i);
#pragma omp parallel
			if (!currentElement->isOptimized() || forceRecomputation) {
				optimizePartitionElement(currentElement, i + 1,
						scheme->getNumberOfElements());
			}
#pragma omp end parallel
		}

		//pString = scheme->getLk();
		notify_observers(MT_SCHEME_END, nullId, time(NULL), current_index,
				max_index, pString);
	}

	return 0;
}

} /* namespace partest */
