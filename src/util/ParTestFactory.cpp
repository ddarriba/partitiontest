/*
 * ParTestFactory.cpp
 *
 *  Created on: Mar 8, 2013
 *      Author: diego
 */

#include "ParTestFactory.h"
#ifdef _PLL
#include "exe/PLLModelOptimize.h"
#else
#include "exe/PhymlModelOptimize.h"
#endif
#include "search/ExhaustiveSearchAlgorithm.h"
#include "search/RandomSearchAlgorithm.h"
#include "search/GreedySearchAlgorithm.h"
#include "search/HierarchicalSearchAlgorithm.h"

namespace partest {

ParTestFactory::ParTestFactory() {
	// TODO Auto-generated constructor stub

}

ParTestFactory::~ParTestFactory() {
	// TODO Auto-generated destructor stub
}

ModelOptimize * ParTestFactory::createModelOptimize(ParTestOptions * options) {
	ModelOptimize * mo;
#ifdef _PLL
	mo = new PLLModelOptimize(options);
#else
	mo = new PhymlModelOptimize(options);
#endif
	return mo;
}

SearchAlgorithm * ParTestFactory::createSearchAlgorithm(
		ParTestOptions * options, PartitionMap * partitionMap) {
	SearchAlgorithm * searchAlgorithm = 0;
	switch (options->getSearchAlgorithm()) {
	case SearchExhaustive:
		searchAlgorithm = new ExhaustiveSearchAlgorithm(options, partitionMap);
		break;
	case SearchRandom:
		searchAlgorithm = new RandomSearchAlgorithm(options, partitionMap);
		break;
	case SearchGreedy:
		searchAlgorithm = new GreedySearchAlgorithm(options, partitionMap, false);
		break;
	case SearchGreedyExtended:
			searchAlgorithm = new GreedySearchAlgorithm(options, partitionMap, true);
			break;
	case SearchHCluster:
		searchAlgorithm = new HierarchicalSearchAlgorithm(options,
				partitionMap);
		break;
	}
	return searchAlgorithm;
}

} /* namespace partest */
