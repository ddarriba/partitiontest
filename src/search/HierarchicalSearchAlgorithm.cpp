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
 * @file HierarchicalSearchAlgorithm.cpp
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

PartitioningScheme * HierarchicalSearchAlgorithm::start(
		PartitioningScheme * startingPoint) {
	cerr << "[ERROR] Not implemented yet" << endl;
	Utilities::exit_partest(EX_UNAVAILABLE);
	return 0;
}

PartitioningScheme * HierarchicalSearchAlgorithm::start() {

#ifdef DEBUG
	cout << "[TRACE] Hcluster - START" << endl;
#endif

	int current_scheme = 0;

	ModelOptimize * mo = ParTestFactory::createModelOptimize(options);
	ConsoleObserver * observer = new ConsoleObserver();
	this->attach(observer);
	mo->attach(observer);
	mo->attach(this);

	/* 1. start with k=n groups */
	vector<PartitioningScheme *> nextSchemes;
	nextSchemes.push_back(
			new PartitioningScheme(partitionMap->getNumberOfPartitions()));
#ifdef DEBUG
	cout << "[TRACE]            Creating first scheme" << endl;
#endif
	for (int i = 0; i < partitionMap->getNumberOfPartitions(); i++) {
		PartitionElement * nextElement = partitionMap->getPartitionElement(i);
		nextSchemes.at(0)->addElement(nextElement);
	}
#ifdef DEBUG
	cout << "[TRACE]            Created first scheme" << endl;
#endif

	PartitioningScheme * bestScheme = nextSchemes.at(0);
	bool reachedMaximum = false;
	double distance;
	double bestCriterionValue = DOUBLE_INF;
	vector<t_partitionElementId> bestMatch0, bestMatch1;
	t_partitionElementId nullId;
	int currentStep = 1;
	while (!reachedMaximum) {
		reachedMaximum = true;
		notify_observers(MT_NEXT_STEP, nullId, time(NULL), currentStep++,
				numberOfBits);
		for (int i = 0; i < nextSchemes.size(); i++) {
#ifdef DEBUG
			cout << "[TRACE] Hcluster - Next scheme: " << nextSchemes.at(i)->toString() << endl;
#endif
			mo->optimizePartitioningScheme(nextSchemes.at(i), false,
					++current_scheme, partitionMap->getNumberOfPartitions());
#ifdef DEBUG
			cout << "[TRACE] Hcluster - Done scheme: " << nextSchemes.at(i)->toString() << endl;
#endif
			PartitionSelector partSelector(&(nextSchemes.at(i)), 1, options);
			double criterionValue = partSelector.getBestSelectionScheme()->value;
			if (criterionValue < bestCriterionValue) {
				bestCriterionValue = criterionValue;
				bestScheme = partSelector.getBestScheme();
				reachedMaximum = false;
				for (int j = 0; j < bestMatch0.size(); j++) {
					partitionMap->deletePartitionElement(bestMatch0.at(j));
					partitionMap->deletePartitionElement(bestMatch1.at(j));
				}
				bestMatch0.clear();
				bestMatch1.clear();
			}
		}

#ifdef DEBUG
		cout << "[TRACE] Hcluster - Computing distances" << endl;
#endif

		if (nextSchemes.at(0)->getNumberOfElements() > 1 && !reachedMaximum) {

			t_partitionElementId element1, element2;
			bestScheme->getClosestPartitions(element1, element2);
			bestMatch0.push_back(element1);
			bestMatch1.push_back(element2);

#ifdef DEBUG
			cout << "[TRACE] Hcluster - Building next scheme " << endl;
#endif

			PartitioningScheme * prevScheme = bestScheme; //nextSchemes.at(0);
			int size = nextSchemes.size();
			for (int i = 0; i < size; i++) {
				nextSchemes.pop_back();
			}
			int numPartitions = prevScheme->getNumberOfElements() - 1;
			for (int i = 0; i < bestMatch0.size(); i++) {
				nextSchemes.push_back(new PartitioningScheme(numPartitions));
				t_partitionElementId nextId;
				Utilities::mergeIds(nextId, bestMatch0.at(i), bestMatch1.at(i));
				PartitionElement * newElement =
						partitionMap->getPartitionElement(nextId);
				nextSchemes.at(nextSchemes.size() - 1)->addElement(newElement);
				for (int j = 0; j < prevScheme->getNumberOfElements(); j++) {

					PartitionElement * element = prevScheme->getElement(j);
					if (element->getId() != bestMatch0.at(i)
							&& element->getId() != bestMatch1.at(i)) {
						nextSchemes.at(nextSchemes.size() - 1)->addElement(
								element);
					}
				}

			}
			if (prevScheme != bestScheme) {
				delete prevScheme;
			}
		}
	}

	delete mo;
	delete observer;

	return bestScheme;
}

void HierarchicalSearchAlgorithm::update(const ObservableInfo& info,
		ParTestOptions* run_instance) {

}

} /* namespace partest */
