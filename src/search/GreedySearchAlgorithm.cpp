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
 * @file GreedySearchAlgorithm.cpp
 */
#include "GreedySearchAlgorithm.h"
#include "util/ParTestFactory.h"
#include "util/Utilities.h"
#include "exe/PLLModelOptimize.h"
#include "observer/ConsoleObserver.h"
#include "selection/ModelSelector.h"
#include "selection/PartitionSelector.h"
#include "indata/PartitioningScheme.h"
#include "indata/PartitionElement.h"

namespace partest {

/* functor for getting the next partitioning scheme in a group */
struct nextSchemeFunctor {
	nextSchemeFunctor(PartitioningScheme * currentScheme,
			PartitionMap * partitionMap, bool extended) :
			currentScheme(currentScheme), partitionMap(partitionMap), extended(
					extended) {
		numberOfElements = currentScheme->getNumberOfElements();
		i = 0;
		j = 1;
		schemeVector = 0;
		mask = 0;
		maskIndex = 0;
//		t_partitionElementId maxIndex = 0;

		if (extended) {
			// TODO: NOT IMPLEMENTED YET
//			mask = (t_partitionElementId *) malloc(
//					(numberOfElements * (1 + numberOfElements)) / 2
//							* sizeof(t_partitionElementId));
//			for (maskIndex = 0; maskIndex < numberOfElements; maskIndex++) {
//				t_partitionElementId nextId = currentScheme->getElement(
//						maskIndex)->getId();
//				mask[maskIndex] = nextId;
//				if (nextId > maxIndex)
//					maxIndex = nextId;
//			}
//
//			for (i = 0; i < numberOfElements; i++) {
//				for (j = i + 1; j < numberOfElements; j++) {
//
//					t_partitionElementId nextId = mask[i] + mask[j];
//					mask[maskIndex++] = nextId;
//					if (nextId > maxIndex)
//						maxIndex = nextId;
//				}
//			}
//			i = 0;
//			schemeVector = new t_schemesVector;
//			unsigned int bitMask;
//			if (Utilities::isPowerOfTwo(maskIndex)) {
//				bitMask = maskIndex - 1;
//			} else {
//				bitMask = Utilities::binaryPow(Utilities::binaryLog(maxIndex))
//						- 1;
//			}
//			PartitionManager::getPermutations(mask, maskIndex, bitMask,
//					schemeVector);
		}
	}

	~nextSchemeFunctor() {
		if (mask)
			free(mask);
		delete schemeVector;
	}

	PartitioningScheme * operator()(void) {
		if (extended) {

			if (i == schemeVector->size()) {
				return 0;
			}
			t_partitioningScheme elements = schemeVector->at(i);
			PartitioningScheme * nextScheme = new PartitioningScheme(
					elements.size());

			if (!elements.size()) {
				return 0;
			}

			for (j = 0; j < elements.size(); j++) {

				PartitionElement * nextElement =
						partitionMap->getPartitionElement(elements.at(j));
				nextScheme->addElement(nextElement);
			}

			i++;
			return nextScheme;

		} else {

			if (i == (numberOfElements - 1)) {
				return 0;
			}

			PartitioningScheme * nextScheme = new PartitioningScheme(
					numberOfElements - 1);

			for (unsigned int k = 0; k < numberOfElements; k++) {
				if (k != i && k != j) {
					nextScheme->addElement(currentScheme->getElement(k));
				}
			}

			t_partitionElementId nextId;
			Utilities::mergeIds(nextId, currentScheme->getElement(i)->getId(),
					currentScheme->getElement(j)->getId());

			PartitionElement * nextElement = partitionMap->getPartitionElement(
					nextId);

			nextScheme->addElement(nextElement);

			/* increment loop iterators */
			j++;
			if (j == numberOfElements) {
				i++;
				j = i + 1;
			}
			return nextScheme;

		}
	}

	t_partitionElementId * getMask() {
		return mask;
	}

	unsigned int size() {
		if (extended)
			return schemeVector->size();
		else
			return (numberOfElements * (numberOfElements - 1) / 2);
	}

private:
	t_schemesVector * schemeVector;
	t_partitionElementId * mask;
	int maskIndex;
	unsigned int i, j;
	unsigned int numberOfElements;
	PartitioningScheme * currentScheme;
	PartitionMap * partitionMap;
	bool extended;
}
;

GreedySearchAlgorithm::GreedySearchAlgorithm(ParTestOptions * options,
		PartitionMap * partitionMap, bool extended) :
		SearchAlgorithm(options, partitionMap), extended(extended) {

	numberOfBits = this->getNumberOfElements();

}

GreedySearchAlgorithm::~GreedySearchAlgorithm() {
}

PartitioningScheme * GreedySearchAlgorithm::start(
		PartitioningScheme * startingPoint) {
	cerr << "[ERROR] Not implemented yet" << endl;
	Utilities::exit_partest(EX_UNAVAILABLE);
	return 0;
}

PartitioningScheme * GreedySearchAlgorithm::start() {

#ifdef DEBUG
	cout << "[TRACE] Greedy - START" << endl;
#endif

	int i;

	PLLModelOptimize * mo =
			static_cast<PLLModelOptimize *>(ParTestFactory::createModelOptimize(
					options));
	ConsoleObserver * observer = new ConsoleObserver();
	this->attach(observer);
	mo->attach(observer);
	mo->attach(this);

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

	/* 1. start with k=n groups */
	PartitioningScheme * firstScheme = new PartitioningScheme(
			partitionMap->getNumberOfPartitions());
	for (unsigned int i = 0; i < partitionMap->getNumberOfPartitions(); i++) {
		firstScheme->addElement(partitionMap->getPartitionElement(i));
	}

#ifdef DEBUG
	cout << "[TRACE] Greedy - Optimizing first scheme" << endl;
#endif

	mo->optimizePartitioningScheme(firstScheme, false, 1, 1);
	PartitionSelector partSelector(&firstScheme, 1, options);
	double bestCriterionValue = partSelector.getBestSelectionScheme()->value;

	PartitioningScheme * bestScheme = firstScheme;

#ifdef DEBUG
	cout << "[TRACE] Greedy - First Scheme Done " << firstScheme->getName() << endl;
#endif

	bool reachedMaximum = firstScheme->getNumberOfElements() == 1;

	while (!reachedMaximum) {
		/* 2. evaluate next-step partitions */

#ifdef DEBUG
		cout << "[TRACE] Greedy - Next iteration" << endl;
#endif

		nextSchemeFunctor nextScheme(bestScheme, partitionMap, extended);

#ifdef DEBUG
		cout << "[TRACE] Greedy - Created functor" << endl;
#endif

		int numberOfSchemes = nextScheme.size();
		PartitioningScheme **partitioningSchemesVector =
				(PartitioningScheme **) malloc(
						numberOfSchemes * sizeof(PartitioningScheme *));
		i = 0;
		while (PartitioningScheme * currentScheme = nextScheme()) {
			partitioningSchemesVector[i++] = currentScheme;
#ifdef DEBUG
			cout << "[TRACE] Greedy - Optimizing next scheme: " << currentScheme->getName() << endl;
#endif
			mo->optimizePartitioningScheme(currentScheme, false, i,
					nextScheme.size());
		}

		PartitionSelector partSelector(partitioningSchemesVector,
				numberOfSchemes, options);
#ifdef DEBUG
		cout << "[TRACE] Greedy - Best scheme: " << partSelector.getBestScheme()->getName() << endl;
#endif

		reachedMaximum = (bestCriterionValue
				<= partSelector.getBestSelectionScheme()->value);
		if (!reachedMaximum) {
			bestCriterionValue = partSelector.getBestSelectionScheme()->value;
			bestScheme = partSelector.getBestScheme();
			reachedMaximum |= bestScheme->getNumberOfElements() == 1;
		}
#ifdef DEBUG
		cout << "[TRACE] Greedy - Purging schemes" << endl;
#endif
		for (int i = 0; i < bestScheme->getNumberOfElements(); i++) {
			partitionMap->purgePartitionMap(bestScheme->getElement(i)->getId());
//			for (int j=0; j<bestScheme->getNumberOfElements(); j++) {
//				if (bestScheme->getElement(j)->getNumberOfSections() > 0) {
//					for (int k=0; k<bestScheme->getElement(j)->getNumberOfSections()) {
//						bestScheme->getElement(j)->getId().at(k)
//						t_partitionElementId k;
//					}
//				}
//			}
		}
#ifdef DEBUG
		cout << "[TRACE] Greedy - Free schemes vector" << endl;
#endif
		free(partitioningSchemesVector);
	}

#ifdef DEBUG
	cout << "[TRACE] Greedy - DONE " << endl;
#endif

	delete mo;

	return bestScheme;
}

void GreedySearchAlgorithm::update(const ObservableInfo& info,
		ParTestOptions* run_instance) {

}

int GreedySearchAlgorithm::getNextPartitioningSchemes(
		PartitioningScheme *currentScheme, PartitioningScheme **nextSchemes) {

//	int numberOfPartitions = currentScheme->getNumberOfElements()
//			* (currentScheme->getNumberOfElements() - 1);

	return 0;
}

} /* namespace partest */
