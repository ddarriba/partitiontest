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

namespace partest {

ModelOptimize::ModelOptimize(ParTestOptions * options) :
    options(options)
{

}

ModelOptimize::~ModelOptimize()
{
}

//int ModelOptimize::optimizeModel(Model * model,
//		PartitionElement * partitionElement, int index, int groupCount) {
//  // MOCK!
//  return 0;
//}

int ModelOptimize::optimizePartitionElement(PartitionElement * partitionElement) {
		ModelSet * modelset = partitionElement->getModelset();
		//(options->getRateVariation(), options->getDataType());
		notify_observers(MT_MODELSET_INIT, 0, modelset, time(NULL), partitionElement->getId(),
				partitionElement->getId());
#pragma omp for schedule(dynamic)
		for (int i = 0; i < modelset->getNumberOfModels(); i++) {
			optimizeModel(modelset->getModel(i), partitionElement, i, modelset->getNumberOfModels());
		}

//		ModelSelector selector(partitionElement->getModelset(), BIC, 113.0);
//
//		SelectionModel * best = selector.getBestModel();
//
//		cout << " ------------- BEST MODEL --------------- " << endl;
//		cout << best->getModel()->getName() << endl;
//		cout << " ------------- ---------- --------------- " << endl;

		notify_observers(MT_MODELSET_END, 0, modelset, time(NULL), 0,
				modelset->getNumberOfModels());

		return 0;
}

int ModelOptimize::optimizePartitioningScheme(PartitioningScheme * scheme,
		bool forceRecomputation) {
	string * pString = new string(scheme->toString());
	notify_observers(MT_SCHEME_INIT, 0, time(NULL), 0,
			scheme->getNumberOfElements(), pString);
// #pragma omp parallel for schedule(dynamic)

	for (int i = 0; i < scheme->getNumberOfElements(); i++) {
		PartitionElement * currentElement = scheme->getElement(i);
#pragma omp parallel
		if (!currentElement->isOptimized() || forceRecomputation) {
			optimizePartitionElement(currentElement);
		}
#pragma omp end parallel
	}

	notify_observers(MT_SCHEME_END, 0, time(NULL), 0,
			scheme->getNumberOfElements(), pString);
	delete pString;

	return 0;
}

} /* namespace partest */
