/*
 * ParTestFactory.h
 *
 *  Created on: Mar 8, 2013
 *      Author: diego
 */

#ifndef PARTESTFACTORY_H_
#define PARTESTFACTORY_H_

#include "options/ParTestOptions.h"
#include "exe/ModelOptimize.h"
#include "search/SearchAlgorithm.h"
#include "util/GlobalDefs.h"
#ifdef _PLL
#include "exe/PLLModelOptimize.h"
#else
#include "exe/PhymlModelOptimize.h"
#endif

namespace partest {

class ParTestFactory {
public:
	ParTestFactory();
	virtual ~ParTestFactory();
	static ModelOptimize * createModelOptimize( ParTestOptions * options );
	static SearchAlgorithm * createSearchAlgorithm( ParTestOptions * options, PartitionMap * partitionMap );
};

} /* namespace partest */
#endif /* PARTESTFACTORY_H_ */
