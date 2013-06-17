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
 *  For any other inquiries send an Email to Diego Darriba
 *  ddarriba@udc.es
 */

/**
 * @file ExhaustiveSearchAlgorithm.cpp
 */

#include "ExhaustiveSearchAlgorithm.h"

#include "indata/PartitionManager.h"
#include "util/ParTestFactory.h"
#include "exe/ModelOptimize.h"
#include "observer/ConsoleObserver.h"
#include "selection/ModelSelector.h"
#include "selection/PartitionSelector.h"
#include <algorithm>    // std::adjacent_find
#include <iostream>     // std::cout
#include <vector>

namespace partest {

ExhaustiveSearchAlgorithm::ExhaustiveSearchAlgorithm(ParTestOptions * options,
		PartitionMap * partitionMap) :
		SearchAlgorithm(options, partitionMap) {

}

ExhaustiveSearchAlgorithm::~ExhaustiveSearchAlgorithm() {
	// TODO Auto-generated destructor stub
}

PartitioningScheme * ExhaustiveSearchAlgorithm::start(PartitioningScheme * startingPoint) {
	cerr << "[ERROR] Not implemented yet" << endl;
	Utilities::exit_partest(EX_UNAVAILABLE);
	return 0;
}

PartitioningScheme * ExhaustiveSearchAlgorithm::start() {
	int i, j;
	int n = getNumberOfElements();
	int numberOfPartitions;

	//TODO: UNINPLEMENTED
	cerr << "TO BE IMPLEMENTED" << endl;
	Utilities::exit_partest(EX_UNAVAILABLE);
//#ifdef DEBUG
//	cout << "[TRACE] START SEARCH FOR " << n << " ELEMENTS"<< endl;
//#endif
//	PartitionElement * part_element = getPartitionElement(1);
//
//	numberOfPartitions = Utilities::binaryPow(n);
//	t_partitionElementId * mask = (t_partitionElementId *) malloc(
//			numberOfPartitions * sizeof(t_partitionElementId));
//
//	for (i = 0; i < numberOfPartitions; i++) {
//		mask[i] = i;
//	}
//
//	t_schemesVector * schemeVector = new t_schemesVector;
//	PartitionManager::getPermutations(mask, numberOfPartitions,
//			numberOfPartitions - 1, schemeVector);
//	free(mask);
//
//	std::sort(schemeVector->begin(), schemeVector->end());
//
//	std::vector<t_partitioningScheme>::iterator it;
//	it = schemeVector->begin();
//	while ((it = std::adjacent_find(it, schemeVector->end())) != schemeVector->end()) {
//		schemeVector->erase(it);
//	}
//	//std::unique ( partitions->begin(), partitions->end());
//
//#ifdef DEBUG
//	cout << "[TRACE] PARTITIONING SCHEMES = " << schemeVector->size() << endl;
//#endif
//
//	for (i = 0; i < schemeVector->size(); i++) {
//		cout << "( ";
//		for (int j = 0; j < schemeVector->at(i).size(); j++) {
//			cout << schemeVector->at(i)[j] << " ";
//		}
//		cout << ")";
//	}
//
//#ifdef DEBUG
//	cout << "[TRACE] Created partitioning scheme" << endl;
//#endif
//
//	ModelOptimize * mo = ParTestFactory::createModelOptimize(options);
//
//#ifdef DEBUG
//	cout << "[TRACE] Created optimizer" << endl;
//#endif
//
//	ConsoleObserver * observer = new ConsoleObserver();
//	mo->attach(observer);
//	mo->attach(this);
//#ifdef DEBUG
//	cout << "[TRACE] Attached observer" << endl;
//#endif
//
//	numberOfPartitions = schemeVector->size();
//	PartitioningScheme * scheme;
//	PartitioningScheme ** schemesArray = (PartitioningScheme **) malloc(
//			numberOfPartitions * sizeof(PartitioningScheme *));
//	for (i = 0; i < numberOfPartitions; i++) {
//		scheme = new PartitioningScheme(&(schemeVector->at(i)), partitionMap);
//		mo->optimizePartitioningScheme(scheme);
//		schemesArray[i] = scheme;
//
//	}
//
//#ifdef DEBUG
//	cout << "[TRACE] Done optimization" << endl;
//#endif
//
//	PartitionSelector partSelector(schemesArray, numberOfPartitions, options);
//
//#ifdef DEBUG
//	cout << "[TRACE] Done selection" << endl;
//#endif
//
//	for (i = 0; i < numberOfPartitions; i++) {
//		if (schemesArray[i] != partSelector.getBestScheme())
//			delete schemesArray[i];
//	}
//
//	free(schemesArray);
//	delete schemeVector;
//	delete mo;
//	delete observer;
//#ifdef DEBUG
//	cout << "[TRACE] Deleted" << endl;
//	cout << "[TRACE] Returning " << (partSelector.getBestScheme())->toString() << endl;
//#endif
//	return partSelector.getBestScheme();
}

void ExhaustiveSearchAlgorithm::update(const ObservableInfo & info,
		ParTestOptions * run_instance) {
	switch (info.type) {
	case MT_MODELSET_INIT:
		break;
	case MT_MODELSET_END:
		break;
	}
}
} /* namespace partest */
