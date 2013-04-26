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

Partition * ExhaustiveSearchAlgorithm::start() {
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

	t_partition * partitions = new t_partition;
	PartitionManager::get_permutations(mask, numberOfPartitions,
			numberOfPartitions - 1, partitions);
	free(mask);

	std::sort(partitions->begin(), partitions->end());

	std::vector<t_partition_elements>::iterator it;
	it = partitions->begin();
	while ((it = std::adjacent_find(it, partitions->end())) != partitions->end()) {
		partitions->erase(it);
	}
	//std::unique ( partitions->begin(), partitions->end());

#ifdef DEBUG
	cout << "[TRACE] PARTITIONS = " << partitions->size() << endl;
#endif

	for (i = 0; i < partitions->size(); i++) {
		cout << "( ";
		for (int j = 0; j < partitions->at(i).size(); j++) {
			cout << partitions->at(i)[j] << " ";
		}
		cout << ")";
	}

#ifdef DEBUG
	cout << "[TRACE] Created partition" << endl;
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

	numberOfPartitions = partitions->size();
	Partition * part;
	Partition ** partitionsArray = (Partition **) malloc(
			numberOfPartitions * sizeof(Partition *));
	for (i = 0; i < numberOfPartitions; i++) {
		part = new Partition(&(partitions->at(i)), partitionMap);
		mo->optimizePartition(part);
		partitionsArray[i] = part;

	}
	PartitionSelector partSelector(partitionsArray, numberOfPartitions, options->getInformationCriterion(),
			options->getSampleSize(), options->getSampleSizeValue());

	for (i = 0; i < numberOfPartitions; i++) {
		delete partitionsArray[i];
	}

	free(partitionsArray);
	delete partitions;
	delete mo;
	delete observer;

	return partSelector.getBestPartition();
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
