/*
 * GreedySearchAlgorithm.cpp
 *
 *  Created on: Apr 1, 2013
 *      Author: diego
 */

#include "GreedySearchAlgorithm.h"
#include "util/ParTestFactory.h"
#include "util/Utilities.h"
#include "exe/ModelOptimize.h"
#include "observer/ConsoleObserver.h"
#include "selection/ModelSelector.h"
#include "selection/PartitionSelector.h"
#include "indata/PartitioningScheme.h"
#include "indata/PartitionElement.h"
#include "indata/PartitionManager.h"

#define GREEDY_EXTENDED true

namespace partest {

/* functor for getting the next partitioning scheme in a group */
struct nextSchemeFunctor {
	nextSchemeFunctor(PartitioningScheme * currentScheme,
			PartitionMap * partitionMap) :
			currentScheme(currentScheme), partitionMap(partitionMap) {
		numberOfElements = currentScheme->getNumberOfElements();
		i = 0;
		j = 1;
#ifdef GREEDY_EXTENDED
		maskIndex = 0;
		t_partitionElementId maxIndex = 0;
		mask = (t_partitionElementId *) malloc(
				(numberOfElements * (1 + numberOfElements)) / 2
						* sizeof(t_partitionElementId));
		for (maskIndex = 0; maskIndex < numberOfElements; maskIndex++) {
			t_partitionElementId nextId =
					currentScheme->getElement(maskIndex)->getId();
			mask[maskIndex] = nextId;
			if (nextId > maxIndex)
				maxIndex = nextId;
		}

		for (i = 0; i < numberOfElements; i++) {
			for (j = i + 1; j < numberOfElements; j++) {

				t_partitionElementId nextId = mask[i] + mask[j];
				mask[maskIndex++] = nextId;
				if (nextId > maxIndex)
					maxIndex = nextId;
			}
		}
		i = 0;
		schemeVector = new t_schemesVector;
		unsigned int bitMask;
		if (Utilities::isPowerOfTwo(maskIndex)) {
			bitMask = maskIndex - 1;
		} else {
			bitMask = Utilities::binaryPow(Utilities::binaryLog(maxIndex)) - 1;
		}
		PartitionManager::getPermutations(mask, maskIndex, bitMask,
				schemeVector);
#endif
	}

	~nextSchemeFunctor() {
		free(mask);
		delete schemeVector;
	}
#ifdef GREEDY_EXTENDED
	PartitioningScheme * operator()(void) {
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

			PartitionElement * nextElement = partitionMap->getPartitionElement(
					elements.at(j));
			nextScheme->addElement(nextElement);
		}

		i++;
		return nextScheme;
	}
#else
	PartitioningScheme * operator()(void) {
		if (i == (numberOfElements - 1)) {
			return 0;
		}
		PartitioningScheme * nextScheme = new PartitioningScheme(numberOfElements - 1);
		int k;

		for (k = 0; k < numberOfElements; k++) {
			if (k != i && k != j) {
				nextScheme->addElement(currentScheme->getElement(k));
			}
		}

		t_partitionElementId nextId = currentScheme->getElement(i)->getId()
		+ currentScheme->getElement(j)->getId();
		mask[maskIndex++] = nextId;
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
#endif

	t_partitionElementId * getMask() {
		return mask;
	}

	unsigned int size() {
		return schemeVector->size();
	}

private:
	t_schemesVector * schemeVector;
	t_partitionElementId * mask;
	int i, j, maskIndex;
	int numberOfElements;
	PartitioningScheme * currentScheme;
	PartitionMap * partitionMap;
};

GreedySearchAlgorithm::GreedySearchAlgorithm(ParTestOptions * options,
		PartitionMap * partitionMap) :
		SearchAlgorithm(options, partitionMap) {

	numberOfBits = this->getNumberOfElements();

}

GreedySearchAlgorithm::~GreedySearchAlgorithm() {
	// TODO Auto-generated destructor stub
}

PartitioningScheme * GreedySearchAlgorithm::start() {

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
	bool reachedMaximum = firstScheme->getNumberOfElements() == 1;
	while (!reachedMaximum) {
		/* 2. evaluate next-step partitions */

		nextSchemeFunctor nextScheme(bestScheme, partitionMap);

		int numberOfSchemes = nextScheme.size();
		PartitioningScheme **partitioningSchemesVector =
				(PartitioningScheme **) malloc(
						numberOfSchemes * sizeof(PartitioningScheme *));
		i = 0;
		while (PartitioningScheme * currentScheme = nextScheme()) {
			partitioningSchemesVector[i++] = currentScheme;
			mo->optimizePartitioningScheme(currentScheme);
		}
		PartitionSelector partSelector(partitioningSchemesVector,
				numberOfSchemes, options);
		reachedMaximum = (bestCriterionValue
				<= partSelector.getBestSelectionScheme()->value);
		if (!reachedMaximum) {
			bestCriterionValue = partSelector.getBestSelectionScheme()->value;
			bestScheme = partSelector.getBestScheme();
			reachedMaximum |= bestScheme->getNumberOfElements() == 1;
		}

		free(partitioningSchemesVector);
	}

	delete mo;

	return bestScheme;
}

void GreedySearchAlgorithm::update(const ObservableInfo& info,
		ParTestOptions* run_instance) {

}

int GreedySearchAlgorithm::getNextPartitioningSchemes(
		PartitioningScheme *currentScheme, PartitioningScheme **nextSchemes) {
	int numberOfPartitions = currentScheme->getNumberOfElements()
			* (currentScheme->getNumberOfElements() - 1);

}

} /* namespace partest */
