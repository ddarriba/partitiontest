/*
 * RandomSearchAlgorithm.cpp
 *
 *  Created on: Apr 1, 2013
 *      Author: diego
 */

#include "RandomSearchAlgorithm.h"
#include "util/ParTestFactory.h"
#include "exe/ModelOptimize.h"
#include "observer/ConsoleObserver.h"
#include "selection/ModelSelector.h"
#include "selection/PartitionSelector.h"
#include "indata/Partition.h"
#include "indata/PartitionElement.h"

#define RND_THRESHOLD 0.5
#define NUM_STEPS 3
#define NUM_PARTITIONS 20

namespace partest {

RandomSearchAlgorithm::RandomSearchAlgorithm(ParTestOptions * options,
		PartitionMap * partitionMap) :
		SearchAlgorithm(options, partitionMap) {

	numberOfBits = this->getNumberOfElements();

}

RandomSearchAlgorithm::~RandomSearchAlgorithm() {
	// TODO Auto-generated destructor stub
}

Partition * RandomSearchAlgorithm::start() {

	int cur_partition, cur_step;

	ModelOptimize * mo = ParTestFactory::createModelOptimize(options);
	ConsoleObserver * observer = new ConsoleObserver();
	mo->attach(observer);
	mo->attach(this);

	Partition * bestPartition = 0;
	int numPartitions =
			(NUM_PARTITIONS < Utilities::bell(numberOfBits)) ?
					NUM_PARTITIONS : Utilities::bell(numberOfBits);
	for (cur_step = 0; cur_step < NUM_STEPS; cur_step++) {
		Partition ** partitionsArray = (Partition **) malloc(
				numPartitions * sizeof(Partition *));
		if (bestPartition) {
			if (bestPartition->getNumberOfElements() == 1)
				break;
			numPartitions =
					(NUM_PARTITIONS
							< Utilities::bell(
									bestPartition->getNumberOfElements())) ?
							NUM_PARTITIONS :
							Utilities::bell(
									bestPartition->getNumberOfElements());
		}
		for (cur_partition = 0; cur_partition < numPartitions;
				cur_partition++) {

			if (!(partitionsArray[cur_partition] = getRandomPartition(
					bestPartition, partitionsArray, cur_partition))) {
				break;
			}
			mo->optimizePartition(partitionsArray[cur_partition]);
		}
		PartitionSelector partSelector(partitionsArray, numPartitions, options->getInformationCriterion(),
				options->getSampleSize(), options->getSampleSizeValue());

		bestPartition = partSelector.getBestPartition();

		for (cur_partition = 0; cur_partition < numPartitions;
				cur_partition++) {
			if (partitionsArray[cur_partition] != bestPartition)
				delete partitionsArray[cur_partition];
		}
		free(partitionsArray);
	}
	delete mo;
	delete observer;

	return bestPartition;
}

void RandomSearchAlgorithm::update(const ObservableInfo& info,
		ParTestOptions* run_instance) {

}

bool existPartition(t_partitionElementId * classes, int numberOfClasses,
		Partition ** partitionsArray, int numberOfPartitions) {
	for (int i = 0; i < numberOfPartitions; i++) {
		if (partitionsArray[i]->getNumberOfElements() == numberOfClasses) {
			bool equals = true;
			for (int j = 0; j < numberOfClasses; j++) {
				bool existElement = false;
				t_partitionElementId id =
						partitionsArray[i]->getElement(j)->getId();
				for (int k = 0; k < numberOfClasses; k++) {
					existElement |= (id == classes[j]);
				}
				equals &= existElement;
			}
			if (equals) {
				return true;
			}
		}
	}
	return false;
}

Partition * RandomSearchAlgorithm::getRandomPartition(
		Partition ** partitionsArray, int numberOfPartitions) {

	//Partition * p = new Partition();

	bool gotPartition = false;
	int numberOfClasses;
	int i;
	t_partitionElementId classes[numberOfBits];

	while (!gotPartition) {

		numberOfClasses = 1;
		// first element to the first
		classes[0] = 1;
		for (i = 1; i < numberOfBits; i++) {
			int currentClass = 0;
			int assigned = 0;

			//randomly assign each bit to each class
			while (!assigned) {
				double rndNumber = (double) (rand() % 500 / 500.0);
				if (rndNumber < RND_THRESHOLD) {
					classes[currentClass] += Utilities::binaryPow(i);
					assigned = 1;
				} else {
					currentClass++;
					if (currentClass == numberOfClasses) {
						classes[currentClass] = Utilities::binaryPow(i);
						numberOfClasses++;
						assigned = 1;
					}
				}
			}
		}

		gotPartition = !existPartition(classes, numberOfClasses,
				partitionsArray, numberOfPartitions);
	}
	Partition * p = new Partition(numberOfClasses);
	for (i = 0; i < numberOfClasses; i++) {
		p->addElement(partitionMap->getPartitionElement(classes[i]));
	}

	return p;
}

Partition * RandomSearchAlgorithm::getRandomPartition(Partition * p0,
		Partition ** partitionsArray, int numberOfPartitions) {

	if (p0 == 0)
		return getRandomPartition(partitionsArray, numberOfPartitions);

	int maxClasses = p0->getNumberOfElements() - 1;

	bool gotPartition = false;
	int i;
	int numberOfClasses;
	t_partitionElementId classes[maxClasses];

	while (!gotPartition) {

		numberOfClasses = 1;
		// first element to the first
		classes[0] = p0->getElement(0)->getId();

		for (i = 1; i < maxClasses; i++) {
			int currentClass = 0;
			int assigned = 0;

			//randomly assign each bit to each class
			while (!assigned) {
				double rndNumber = (double) (rand() % 500 / 500.0);
				if (rndNumber < RND_THRESHOLD) {

					classes[currentClass] += p0->getElement(i)->getId();
					assigned = 1;
				} else {
					currentClass++;
					if (currentClass == numberOfClasses) {
						classes[currentClass] = p0->getElement(i)->getId();
						numberOfClasses++;
						assigned = 1;
					}
				}
			}
		}

		gotPartition = !existPartition(classes, numberOfClasses,
				partitionsArray, numberOfPartitions);
	}

	Partition * p = new Partition(numberOfClasses);
	for (i = 0; i < numberOfClasses; i++) {
		p->addElement(partitionMap->getPartitionElement(classes[i]));
	}

	return p;
}

} /* namespace partest */
