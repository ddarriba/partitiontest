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
#include "indata/Partition.h"
#include "indata/PartitionElement.h"
#include "indata/PartitionManager.h"

#define GREEDY_EXTENDED true

namespace partest {

/* functor for getting the next partitioning scheme in a group */
struct nextPartitionFunctor {
	nextPartitionFunctor(Partition * currentPartition,
			PartitionMap * partitionMap) :
			currentPartition(currentPartition), partitionMap(partitionMap) {
		numberOfElements = currentPartition->getNumberOfElements();
		i = 0;
		j = 1;
#ifdef GREEDY_EXTENDED
		maskIndex = 0;
		t_partitionElementId maxIndex = 0;
		mask = (t_partitionElementId *) malloc(
				(numberOfElements * (1 + numberOfElements)) / 2
						* sizeof(t_partitionElementId));
		for (maskIndex = 0; maskIndex < numberOfElements; maskIndex++) {
			t_partitionElementId nextId = currentPartition->getElement(maskIndex)->getId();
			mask[maskIndex] = nextId;
			if (nextId > maxIndex) maxIndex = nextId;
		}

		for (i = 0; i < numberOfElements; i++) {
			for (j = i + 1; j < numberOfElements; j++) {

				t_partitionElementId nextId = mask[i] + mask[j];
				mask[maskIndex++] = nextId;
				if (nextId > maxIndex) maxIndex = nextId;
			}
		}
		i=0;
		partitions = new t_partition;
		unsigned int bitMask;
		if (Utilities::isPowerOfTwo(maskIndex)) {
			bitMask = maskIndex-1;
		} else {
			bitMask = Utilities::binaryPow(Utilities::binaryLog(maxIndex)) - 1;
		}
		PartitionManager::get_permutations(mask, maskIndex, bitMask,
				partitions);
#endif
	}

	~nextPartitionFunctor() {
		free(mask);
	}
#ifdef GREEDY_EXTENDED
	Partition * operator()(void) {
		if (i == partitions->size()) {
			return 0;
		}
		t_partition_elements elements = partitions->at(i);
		Partition * nextPartition = new Partition(elements.size());

		for (j = 0; j < elements.size(); j++) {

			PartitionElement * nextElement = partitionMap->getPartitionElement(
					elements.at(j));
			nextPartition->addElement(nextElement);
		}

		i++;
		return nextPartition;
	}
#else
	Partition * operator()(void) {
		if (i == (numberOfElements - 1)) {
			return 0;
		}
		Partition * nextPartition = new Partition(numberOfElements - 1);
		int k;

		for (k = 0; k < numberOfElements; k++) {
			if (k != i && k != j) {
				nextPartition->addElement(currentPartition->getElement(k));
			}
		}

		t_partitionElementId nextId = currentPartition->getElement(i)->getId()
		+ currentPartition->getElement(j)->getId();
		mask[maskIndex++] = nextId;
		PartitionElement * nextElement = partitionMap->getPartitionElement(
				nextId);
		nextPartition->addElement(nextElement);

		/* increment loop iterators */
		j++;
		if (j == numberOfElements) {
			i++;
			j = i + 1;
		}
		return nextPartition;
	}
#endif

	t_partitionElementId * getMask() {
		return mask;
	}

	unsigned int size() {
		return partitions->size();
	}

private:
	t_partition * partitions;
	t_partitionElementId * mask;
	int i, j, maskIndex;
	int numberOfElements;
	Partition * currentPartition;
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

Partition * GreedySearchAlgorithm::start() {

	int i;

	ModelOptimize * mo = ParTestFactory::createModelOptimize(options);
	ConsoleObserver * observer = new ConsoleObserver();
	mo->attach(observer);
	mo->attach(this);

	/* 1. start with k=n groups */
	Partition * firstPartition = new Partition(
			partitionMap->getNumberOfPartitions());
	for (i = 0; i < partitionMap->getNumberOfPartitions(); i++) {
		firstPartition->addElement(
				partitionMap->getPartitionElement(Utilities::binaryPow(i)));
	}

	mo->optimizePartition(firstPartition);
	PartitionSelector partSelector(&firstPartition, 1,
			options->getInformationCriterion(), options->getSampleSize(),
			options->getSampleSizeValue());
	double bestCriterionValue = partSelector.getBestSelectionPartition()->value;

	Partition * bestPartition = firstPartition;
	bool reachedMaximum = false;
	while (!reachedMaximum) {
		/* 2. evaluate next-step partitions */
		nextPartitionFunctor nextPartition(bestPartition, partitionMap);
		int numberOfPartitions = nextPartition.size();
		Partition **partitionsVector = (Partition **) malloc(
				numberOfPartitions * sizeof(Partition *));
		i = 0;

		while (Partition * currentPartition = nextPartition()) {
			partitionsVector[i++] = currentPartition;
			mo->optimizePartition(currentPartition);
		}

		PartitionSelector partSelector(partitionsVector, numberOfPartitions,
				options->getInformationCriterion(), options->getSampleSize(),
				options->getSampleSizeValue());
		reachedMaximum = (bestCriterionValue
				<= partSelector.getBestSelectionPartition()->value);
		if (!reachedMaximum) {
			bestCriterionValue =
					partSelector.getBestSelectionPartition()->value;
			bestPartition = partSelector.getBestPartition();
			reachedMaximum |= bestPartition->getNumberOfElements() == 1;
		}

		free(partitionsVector);
	}

	return bestPartition;
}

void GreedySearchAlgorithm::update(const ObservableInfo& info,
		ParTestOptions* run_instance) {

}

int GreedySearchAlgorithm::getNextPartitions(Partition *currentPartition,
		Partition **nextPartitions) {
	int numberOfPartitions = currentPartition->getNumberOfElements()
			* (currentPartition->getNumberOfElements() - 1);

}

} /* namespace partest */
