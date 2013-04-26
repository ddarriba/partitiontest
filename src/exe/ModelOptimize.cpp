/*
 * ModelOptimize.cpp
 *
 *  Created on: Jan 8, 2013
 *      Author: diego
 */

#include "ModelOptimize.h"
//#include "selection/ModelSelector.h"

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
#pragma omp parallel for schedule(dynamic) //private(partition)
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

int ModelOptimize::optimizePartition(Partition * partition,
		bool forceRecomputation) {
	string * pString = new string(partition->toString());
	notify_observers(MT_PARTITION_INIT, 0, time(NULL), 0,
			partition->getNumberOfElements(), pString);
	for (int i = 0; i < partition->getNumberOfElements(); i++) {
		PartitionElement * currentElement = partition->getElement(i);
		if (!currentElement->isOptimized() || forceRecomputation) {
			optimizePartitionElement(currentElement);
		}

	}
	notify_observers(MT_PARTITION_END, 0, time(NULL), 0,
			partition->getNumberOfElements(), pString);
	delete pString;

	return 0;
}

} /* namespace partest */
