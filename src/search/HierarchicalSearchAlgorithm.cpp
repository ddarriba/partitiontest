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

namespace partest {

char convert2(unsigned char c) {
	switch (c) {
	case 1:
		return 'A';
	case 2:
		return 'C';
	case 4:
		return 'G';
	case 8:
		return 'T';
	default:
		return 'X';
	}
}

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

	PLLModelOptimize * mo =
			static_cast<PLLModelOptimize *>(ParTestFactory::createModelOptimize(
					options));
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
		for (unsigned int i = 0; i < partitionMap->getNumberOfPartitions(); i++) {
			PartitionElement * nextElement = partitionMap->getPartitionElement(i);
			nextSchemes.at(0)->addElement(nextElement);
		}
	#ifdef DEBUG
		cout << "[TRACE]            Created first scheme" << endl;
	#endif

	if (options->getStartingTopology() == StartTopoFIXED) {

		/* Starting topology */
		notify_observers(MT_FTREE_INIT, time(NULL));

		PLLAlignment * alignment =
				static_cast<PLLAlignment *>(options->getAlignment());

		mo->initializeStructs(alignment->getTree(), alignment->getPartitions(),
				alignment->getPhylip());

		pllComputeRandomizedStepwiseAdditionParsimonyTree(alignment->getTree(),
				alignment->getPartitions());

		mo->evaluateSPR(alignment->getTree(), alignment->getPartitions(), true, true);

		pllTreeToNewick(alignment->getTree()->tree_string, alignment->getTree(),
				alignment->getPartitions(), alignment->getTree()->start->back,
				PLL_TRUE, PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE,
				PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);
		char * startingTree = alignment->getTree()->tree_string;
		startingTree[strlen(startingTree)-1] = '\0';
		options->setTreeString(startingTree);

		notify_observers(MT_FTREE_END, time(NULL), startingTree);

	}


	PartitioningScheme * bestScheme = nextSchemes.at(0);
	bool reachedMaximum = false;
	double bestCriterionValue = DOUBLE_INF;
	t_partitionElementId nullId, toMerge1, toMerge2;
	int currentStep = 1;

	while (!reachedMaximum) {
		reachedMaximum = true;
		notify_observers(MT_NEXT_STEP, nullId, time(NULL), currentStep++,
				numberOfBits);
		for (unsigned int i = 0; i < nextSchemes.size(); i++) {
			mo->optimizePartitioningScheme(nextSchemes.at(i), false,
					++current_scheme, partitionMap->getNumberOfPartitions());
		}

		PartitionSelector partSelector(&(nextSchemes.at(0)), nextSchemes.size(), options);
					double criterionValue = partSelector.getBestSelectionScheme()->value;

					if (criterionValue < bestCriterionValue) {
						bestCriterionValue = criterionValue;
						bestScheme = partSelector.getBestScheme();
						reachedMaximum = false;
					}

		for (unsigned int i = 0; i < nextSchemes.size(); i++) {
			for (unsigned int j = 0; j < nextSchemes.at(i)->getNumberOfElements(); j++) {
				t_partitionElementId t1 = nextSchemes.at(i)->getElement(j)->getId();
				if (t1.size() > 1) {
					for (unsigned int k=0; k<bestScheme->getNumberOfElements(); k++) {
						t_partitionElementId t2 = bestScheme->getElement(k)->getId();
						if (t1 == t2) break;
						if (Utilities::intersec(t1, t2)) {
							partitionMap->deletePartitionElement(t1);
							break;
						}
					}
				}
			}
		}

#ifdef DEBUG
		cout << "[TRACE] Hcluster - Computing distances" << endl;
#endif

		if (bestScheme->getNumberOfElements() > 1 && !reachedMaximum) {

			//t_partitionElementId element1, element2;
			vector<elementPair *> * eps = bestScheme->getElementDistances();

#ifdef DEBUG
			cout << "[TRACE] Hcluster - Building next scheme " << endl;
#endif

			PartitioningScheme * prevScheme = bestScheme; //nextSchemes.at(0);
			int size = nextSchemes.size();
			for (int i = 0; i < size; i++) {
				nextSchemes.pop_back();
			}

			int numPartitions = prevScheme->getNumberOfElements() - 1;
			for (int i = 0; i < min(options->getMaxSamples(), (int)eps->size()); i++) {
				nextSchemes.push_back(new PartitioningScheme(numPartitions));
				t_partitionElementId nextId;
				Utilities::mergeIds(nextId, eps->at(i)->e1->getId(), eps->at(i)->e2->getId());
				PartitionElement * newElement =
						partitionMap->getPartitionElement(nextId);
				nextSchemes.at(nextSchemes.size() - 1)->addElement(newElement);
				for (int j = 0; j < prevScheme->getNumberOfElements(); j++) {

					PartitionElement * element = prevScheme->getElement(j);
					if (element->getId() != eps->at(i)->e1->getId()
							&& element->getId() != eps->at(i)->e2->getId()) {
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
