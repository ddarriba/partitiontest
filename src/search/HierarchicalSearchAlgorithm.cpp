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

#define DOUBLE_INF 1e140

HierarchicalSearchAlgorithm::HierarchicalSearchAlgorithm(
		ParTestOptions * options, PartitionMap * partitionMap) :
		SearchAlgorithm(options, partitionMap) {

	numberOfBits = this->getNumberOfElements();

}

HierarchicalSearchAlgorithm::~HierarchicalSearchAlgorithm() {
	// TODO Auto-generated destructor stub
}

PartitioningScheme * HierarchicalSearchAlgorithm::start() {

#ifdef DEBUG
		cout << "[TRACE] Hcluster - START" << endl;
#endif

	int i;

	ModelOptimize * mo = ParTestFactory::createModelOptimize(options);
	ConsoleObserver * observer = new ConsoleObserver();
	mo->attach(observer);
	mo->attach(this);

	/* 1. start with k=n groups */
	PartitioningScheme * nextScheme = new PartitioningScheme(
			partitionMap->getNumberOfPartitions());
	for (i = 0; i < partitionMap->getNumberOfPartitions(); i++) {
		nextScheme->addElement(
				partitionMap->getPartitionElement(Utilities::binaryPow(i)));
	}

	PartitioningScheme * bestScheme = nextScheme;
	bool reachedMaximum = false;
	double distance;
	double bestCriterionValue = DOUBLE_INF;

	while (!reachedMaximum) {
		reachedMaximum = true;

#ifdef DEBUG
		cout << "[TRACE] Hcluster - Next scheme: " << nextScheme->toString() << endl;
#endif
		mo->optimizePartitioningScheme(nextScheme);
		PartitionSelector partSelector(&nextScheme, 1, options);
		double criterionValue = partSelector.getBestSelectionScheme()->value;
		if (criterionValue < bestCriterionValue) {
			bestCriterionValue = criterionValue;
			bestScheme = partSelector.getBestScheme();
			reachedMaximum = false;
		}

#ifdef DEBUG
		cout << "[TRACE] Hcluster - Computing distances" << endl;
#endif

		if (nextScheme->getNumberOfElements() > 1 && !reachedMaximum) {
			t_partitionElementId bestMatch0, bestMatch1;
			distance = DOUBLE_INF;
			/*    Look for the minimum distance */
			for (int i = 1; i < bestScheme->getNumberOfElements(); i++) {
				Model * m1 =
						bestScheme->getElement(i)->getBestModel()->getModel();
				for (int j = 0; j < i; j++) {
					Model * m2 =
							bestScheme->getElement(j)->getBestModel()->getModel();
					double newDist = m1->distanceTo(m2);
					if (newDist < distance) {
						bestMatch0 = bestScheme->getElement(i)->getId();
						bestMatch1 = bestScheme->getElement(j)->getId();
						distance = newDist;
					}
				}
			}

#ifdef DEBUG
		cout << "[TRACE] Hcluster - Building next scheme" << endl;
#endif

			PartitioningScheme * prevScheme = nextScheme;
			int numPartitions = prevScheme->getNumberOfElements() - 1;
			nextScheme = new PartitioningScheme(numPartitions);
			nextScheme->addElement(partitionMap->getPartitionElement(bestMatch0 + bestMatch1));
			for (int i = 0; i < prevScheme->getNumberOfElements(); i++) {
				PartitionElement * element = prevScheme->getElement(i);
				if (element->getId() != bestMatch0 && element->getId() != bestMatch1) {
					nextScheme->addElement(element);
				}
			}
			delete prevScheme;
		}
	}

	delete mo;
#ifdef DEBUG
		cout << "[TRACE] Hcluster - END Return: " << bestScheme->toString() << endl;
#endif
	return bestScheme;
}

void HierarchicalSearchAlgorithm::update(const ObservableInfo& info,
		ParTestOptions* run_instance) {

}

} /* namespace partest */
