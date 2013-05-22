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
 * @file RandomSearchAlgorithm.cpp
 */


#include "RandomSearchAlgorithm.h"
#include "util/ParTestFactory.h"
#include "exe/ModelOptimize.h"
#include "observer/ConsoleObserver.h"
#include "selection/ModelSelector.h"
#include "selection/PartitionSelector.h"
#include "indata/PartitioningScheme.h"
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

PartitioningScheme * RandomSearchAlgorithm::start() {

	int currentScheme, currentStep;

	ModelOptimize * mo = ParTestFactory::createModelOptimize(options);
	ConsoleObserver * observer = new ConsoleObserver();
	mo->attach(observer);
	mo->attach(this);

	PartitioningScheme * bestScheme = 0;
	int numSchemes =
			(NUM_PARTITIONS < Utilities::bell(numberOfBits)) ?
					NUM_PARTITIONS : Utilities::bell(numberOfBits);
	for (currentStep = 0; currentStep < NUM_STEPS; currentStep++) {
		PartitioningScheme ** schemesArray = (PartitioningScheme **) malloc(
				numSchemes * sizeof(PartitioningScheme *));
		if (bestScheme) {
			if (bestScheme->getNumberOfElements() == 1)
				break;
			numSchemes =
					(NUM_PARTITIONS
							< Utilities::bell(bestScheme->getNumberOfElements())) ?
							NUM_PARTITIONS :
							Utilities::bell(bestScheme->getNumberOfElements());
		}
		for (currentScheme = 0; currentScheme < numSchemes; currentScheme++) {

			if (!(schemesArray[currentScheme] = getRandomPartitioningScheme(
					bestScheme, schemesArray, currentScheme))) {
				break;
			}
			mo->optimizePartitioningScheme(schemesArray[currentScheme]);
		}
		PartitionSelector partSelector(schemesArray, numSchemes, options);

		bestScheme = partSelector.getBestScheme();

		for (currentScheme = 0; currentScheme < numSchemes; currentScheme++) {
			if (schemesArray[currentScheme] != bestScheme)
				delete schemesArray[currentScheme];
		}
		free(schemesArray);
	}
	delete mo;
	delete observer;

	return bestScheme;
}

void RandomSearchAlgorithm::update(const ObservableInfo& info,
		ParTestOptions* run_instance) {

}

bool existPartition(t_partitionElementId * classes, int numberOfClasses,
		PartitioningScheme ** schemesArray, int numberOfSchemes) {
	for (int i = 0; i < numberOfSchemes; i++) {
		if (schemesArray[i]->getNumberOfElements() == numberOfClasses) {
			bool equals = true;
			for (int j = 0; j < numberOfClasses; j++) {
				bool existElement = false;
				t_partitionElementId id =
						schemesArray[i]->getElement(j)->getId();
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

PartitioningScheme * RandomSearchAlgorithm::getRandomPartitioningScheme(
		PartitioningScheme ** schemesArray, int numberOfSchemes) {

	bool gotScheme = false;
	int numberOfClasses;
	int i;
	t_partitionElementId classes[numberOfBits];

	while (!gotScheme) {

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

		gotScheme = !existPartition(classes, numberOfClasses, schemesArray,
				numberOfSchemes);
	}
	PartitioningScheme * p = new PartitioningScheme(numberOfClasses);
	for (i = 0; i < numberOfClasses; i++) {
		p->addElement(partitionMap->getPartitionElement(classes[i]));
	}

	return p;
}

PartitioningScheme * RandomSearchAlgorithm::getRandomPartitioningScheme(
		PartitioningScheme * p0, PartitioningScheme ** schemesArray,
		int numberOfSchemes) {

	if (p0 == 0)
		return getRandomPartitioningScheme(schemesArray, numberOfSchemes);

	int maxClasses = p0->getNumberOfElements() - 1;

	bool gotScheme = false;
	int i;
	int numberOfClasses;
	t_partitionElementId classes[maxClasses];

	while (!gotScheme) {

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

		gotScheme = !existPartition(classes, numberOfClasses, schemesArray,
				numberOfSchemes);
	}

	PartitioningScheme * p = new PartitioningScheme(numberOfClasses);
	for (i = 0; i < numberOfClasses; i++) {
		p->addElement(partitionMap->getPartitionElement(classes[i]));
	}

	return p;
}

} /* namespace partest */
