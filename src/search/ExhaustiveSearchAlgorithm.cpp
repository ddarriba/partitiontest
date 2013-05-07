/*
 * ExhaustiveSearchAlgorithm.cpp
 *
 *  Created on: Mar 8, 2013
 *      Author: diego
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

PartitioningScheme * ExhaustiveSearchAlgorithm::start() {
	int i, j;
	int n = getNumberOfElements();
	int numberOfPartitions;

#ifdef DEBUG
	cout << "[TRACE] START SEARCH FOR " << n << " ELEMENTS"<< endl;
#endif
	PartitionElement * part_element = getPartitionElement(1);

	numberOfPartitions = Utilities::binaryPow(n);
	t_partitionElementId * mask = (t_partitionElementId *) malloc(
			numberOfPartitions * sizeof(t_partitionElementId));

	for (i = 0; i < numberOfPartitions; i++) {
		mask[i] = i;
	}

	t_schemesVector * schemeVector = new t_schemesVector;
	PartitionManager::getPermutations(mask, numberOfPartitions,
			numberOfPartitions - 1, schemeVector);
	free(mask);

	std::sort(schemeVector->begin(), schemeVector->end());

	std::vector<t_partitioningScheme>::iterator it;
	it = schemeVector->begin();
	while ((it = std::adjacent_find(it, schemeVector->end())) != schemeVector->end()) {
		schemeVector->erase(it);
	}
	//std::unique ( partitions->begin(), partitions->end());

#ifdef DEBUG
	cout << "[TRACE] PARTITIONING SCHEMES = " << schemeVector->size() << endl;
#endif

	for (i = 0; i < schemeVector->size(); i++) {
		cout << "( ";
		for (int j = 0; j < schemeVector->at(i).size(); j++) {
			cout << schemeVector->at(i)[j] << " ";
		}
		cout << ")";
	}

#ifdef DEBUG
	cout << "[TRACE] Created partitioning scheme" << endl;
#endif

	ModelOptimize * mo = ParTestFactory::createModelOptimize(options);

#ifdef DEBUG
	cout << "[TRACE] Created optimizer" << endl;
#endif

	ConsoleObserver * observer = new ConsoleObserver();
	mo->attach(observer);
	mo->attach(this);
#ifdef DEBUG
	cout << "[TRACE] Attached observer" << endl;
#endif

	numberOfPartitions = schemeVector->size();
	PartitioningScheme * scheme;
	PartitioningScheme ** schemesArray = (PartitioningScheme **) malloc(
			numberOfPartitions * sizeof(PartitioningScheme *));
	for (i = 0; i < numberOfPartitions; i++) {
		scheme = new PartitioningScheme(&(schemeVector->at(i)), partitionMap);
		mo->optimizePartitioningScheme(scheme);
		schemesArray[i] = scheme;

	}
	PartitionSelector partSelector(schemesArray, numberOfPartitions, options);

	for (i = 0; i < numberOfPartitions; i++) {
		delete schemesArray[i];
	}

	free(schemesArray);
	delete schemeVector;
	delete mo;
	delete observer;

	return partSelector.getBestScheme();
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
