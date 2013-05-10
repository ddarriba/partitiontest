/*
 * HierarchicalSearchAlgorithm.cpp
 *
 *  Created on: Apr 1, 2013
 *      Author: diego
 */

#include "HierarchicalSearchAlgorithm.h"
#include "util/ParTestFactory.h"
#include "util/Utilities.h"
#include "exe/ModelOptimize.h"
#include "observer/ConsoleObserver.h"
#include "selection/ModelSelector.h"
#include "selection/PartitionSelector.h"
#include "indata/PartitioningScheme.h"
#include "indata/PartitionElement.h"
#include "indata/PartitionManager.h"

namespace partest {

HierarchicalSearchAlgorithm::HierarchicalSearchAlgorithm(
		ParTestOptions * options, PartitionMap * partitionMap) :
		SearchAlgorithm(options, partitionMap) {

	numberOfBits = this->getNumberOfElements();

}

HierarchicalSearchAlgorithm::~HierarchicalSearchAlgorithm() {
	// TODO Auto-generated destructor stub
}

PartitioningScheme * HierarchicalSearchAlgorithm::start() {

	int i;

	ModelOptimize * mo = ParTestFactory::createModelOptimize(options);
	ConsoleObserver * observer = new ConsoleObserver();
	mo->attach(observer);
	mo->attach(this);

	/* 1. start with k=n groups */
	PartitioningScheme * firstScheme = new PartitioningScheme(
			partitionMap->getNumberOfPartitions());
	for (i = 0; i < partitionMap->getNumberOfPartitions(); i++) {
		firstScheme->addElement(
				partitionMap->getPartitionElement(Utilities::binaryPow(i)));
	}

	mo->optimizePartitioningScheme(firstScheme);
	PartitionSelector partSelector(&firstScheme, 1, options);
	double bestCriterionValue = partSelector.getBestSelectionScheme()->value;

	PartitioningScheme * bestScheme = firstScheme;
	bool reachedMaximum = true;
	double distance;
	if (firstScheme->getNumberOfElements() > 1) {
		distance =
				bestScheme->getElement(0)->getBestModel()->getModel()->distanceTo(
						bestScheme->getElement(1)->getBestModel()->getModel());
		reachedMaximum = false;
	}
	while (!reachedMaximum) {
		t_partitionElementId * bestMatch0, bestMatch1;
		/* 2. evaluate next-step partitions */

		/*    Look for the minimum distance */
		reachedMaximum = true;
		for (int i = 1; i < bestScheme->getNumberOfElements(); i++) {
			Model * m1 = bestScheme->getElement(i)->getBestModel()->getModel();
			for (int j = 0; j < i; j++) {
				Model * m2 =
						bestScheme->getElement(j)->getBestModel()->getModel();
				double newDist = m1->distanceTo(m2);
				if (newDist < distance) {
					bestMatch0 = bestScheme->getElement(i)->getId();
					bestMatch1 = bestScheme->getElement(j)->getId();
					distance = newDist;
					reachedMaximum = false;
				}
			}
		}

		/* Merge 2 best elements */

	}

	delete mo;

	return bestScheme;
}

void HierarchicalSearchAlgorithm::update(const ObservableInfo& info,
		ParTestOptions* run_instance) {

}

} /* namespace partest */
