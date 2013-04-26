/*
 * ModelOptimize.h
 *
 *  Created on: Jan 8, 2013
 *      Author: Diego Darriba
 */

/** @file ModelOptimize.h
 * @brief Generic abstract class for optimize models
 */
#ifndef MODELOPTIMIZE_H_
#define MODELOPTIMIZE_H_

#include "../model/Model.h"
#include "../indata/Partition.h"
#include "../indata/PartitionElement.h"
#include "../options/ParTestOptions.h"
#include "../observer/Observable.h"

namespace partest {

class ModelOptimize: public Observable {
public:
	ModelOptimize(ParTestOptions * options);
	virtual int optimizePartition(Partition * partition,
			bool forceRecomputation = false);
	virtual int optimizePartitionElement(PartitionElement * partitionElement);
	virtual int optimizeModel(Model * model,
			PartitionElement * partitionElement, int index, int groupCount) = 0;
	virtual ~ModelOptimize();
protected:
	ParTestOptions * options; /** Genomictest options*/
};

} /* namespace partest */
#endif /* MODELOPTIMIZE_H_ */
