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

	int i;

	ModelOptimize * mo = ParTestFactory::createModelOptimize(options);
	ConsoleObserver * observer = new ConsoleObserver();
	mo->attach(observer);
	mo->attach(this);

	/* 1. start with k=n groups */
	vector<PartitioningScheme *> nextSchemes;
	nextSchemes.push_back(
			new PartitioningScheme(partitionMap->getNumberOfPartitions()));
#ifdef DEBUG
	cout << "[TRACE]            Creating first scheme" << endl;
#endif
	for (i = 0; i < partitionMap->getNumberOfPartitions(); i++) {
		PartitionElement * nextElement = partitionMap->getPartitionElement(
				Utilities::binaryPow(i));
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

	while (!reachedMaximum) {
		reachedMaximum = true;
		for (int i = 0; i < nextSchemes.size(); i++) {
#ifdef DEBUG
			cout << "[TRACE] Hcluster - Next scheme: " << nextSchemes.at(i)->toString() << endl;
#endif
			mo->optimizePartitioningScheme(nextSchemes.at(i));
#ifdef DEBUG
			cout << "[TRACE] Hcluster - Done scheme: " << nextSchemes.at(i)->toString() << endl;
#endif
			PartitionSelector partSelector(&(nextSchemes.at(i)), 1, options);
			double criterionValue = partSelector.getBestSelectionScheme()->value;
			if (criterionValue < bestCriterionValue) {
				bestCriterionValue = criterionValue;
				bestScheme = partSelector.getBestScheme();
				reachedMaximum = false;
				for (int j=0; j<bestMatch0.size(); j++) {
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

			t_partitionElementId * closestElements = bestScheme->getClosestPartitions();
			bestMatch0.push_back(closestElements[0]);
			bestMatch1.push_back(closestElements[1]);
			free(closestElements);

#ifdef DEBUG
			cout << "[TRACE] Hcluster - Building next scheme " << endl;
#endif

			PartitioningScheme * prevScheme = bestScheme; //nextSchemes.at(0);

			for (int i = 0; i < nextSchemes.size(); i++)
				nextSchemes.pop_back();
			int numPartitions = prevScheme->getNumberOfElements() - 1;
			for (int i = 0; i < bestMatch0.size(); i++) {
				nextSchemes.push_back(new PartitioningScheme(numPartitions));
				nextSchemes.at(nextSchemes.size() - 1)->addElement(
						partitionMap->getPartitionElement(
								bestMatch0.at(i) + bestMatch1.at(i)));
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

	bestScheme->buildCompleteModelSet();
#ifdef DEBUG
	cout << "[TRACE] Hcluster - OPTIMIZING BEST SCHEME: " << bestScheme->toString() << endl;
#endif
	mo->optimizePartitioningScheme(bestScheme);
	PartitionSelector partSelector(&bestScheme, 1, options);
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
