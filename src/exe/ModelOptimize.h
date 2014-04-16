/*
 * ModelOptimize.h
 *
 *  Created on: Apr 9, 2014
 *      Author: diego
 */

#ifndef MODELOPTIMIZE_H_
#define MODELOPTIMIZE_H_

#include "util/GlobalDefs.h"
#include "indata/PartitioningScheme.h"
#include "indata/PartitionElement.h"

namespace partest {

class ModelOptimize {
public:
	ModelOptimize();
	virtual ~ModelOptimize();

	int optimizePartitioningScheme(PartitioningScheme * scheme, int index=0, int limit=1);
	int optimizePartitionElement(PartitionElement * scheme, int index=0, int limit=1);
private:
	void optimizeModel(PartitionElement * element, unsigned int modelIndex, int limit);
	void setModelParameters(Model * _model, pllInstance * _tree,
			partitionList * _partitions, pllAlignmentData * _alignData, int index, bool setAlphaFreqs);
};

} /* namespace partest */

#endif /* MODELOPTIMIZE_H_ */
