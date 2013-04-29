/*
 * SearchAlgorithm.h
 *
 *  Created on: Mar 8, 2013
 *      Author: diego
 */

#ifndef SEARCHALGORITHM_H_
#define SEARCHALGORITHM_H_

#include "../indata/PartitionMap.h"
#include "../indata/Partition.h"
#include "../indata/PartitionElement.h"
#include "../model/ModelSet.h"
#include "../options/ParTestOptions.h"
#include "../util/GlobalDefs.h"
#include "../indata/Alignment.h"
#include "../observer/Observer.h"

namespace partest {

class SearchAlgorithm : public Observer {
public:
	SearchAlgorithm(ParTestOptions * options, PartitionMap * partitionMap);
	virtual ~SearchAlgorithm();
	virtual PartitioningScheme * start() = 0;
	virtual void update(const ObservableInfo & info,
	        ParTestOptions * run_instance = NULL) = 0;
protected:
	PartitionElement * getPartitionElement(t_partitionElementId id);
	unsigned int getNumberOfElements();
	ParTestOptions * options;
	PartitionMap * partitionMap;
};

} /* namespace partest */
#endif /* SEARCHALGORITHM_H_ */
