/*
 * SearchAlgorithm.cpp
 *
 *  Created on: Apr 16, 2014
 *      Author: diego
 */

#include "SearchAlgorithm.h"

#include "indata/PartitionMap.h"
#include <iostream>
#include <pthread.h>
#include <thread>
#include <memory>
#include <unistd.h>

namespace partest {

SearchAlgorithm::SchemeManager::SchemeManager() {
	nextSchemes = new vector<PartitioningScheme *>();
}

SearchAlgorithm::SchemeManager::~SchemeManager() {
	delete nextSchemes;
}

int SearchAlgorithm::SchemeManager::addSchemes(
		vector<PartitioningScheme *> schemesToAdd) {
	(*nextSchemes) = schemesToAdd;
	return nextSchemes->size();
}

int SearchAlgorithm::SchemeManager::addScheme(
		PartitioningScheme * schemeToAdd) {
	nextSchemes->push_back(schemeToAdd);
	return nextSchemes->size();
}

#ifdef _MPI
void * distribute(void * arg) {
	vector<PartitioningScheme *> * nextSchemes =
			(vector<PartitioningScheme *> *) arg;

	if (numProcs > 1) {
		int buf;
		MPI_Status targetStatus;
		for (unsigned int i = 0; i < nextSchemes->size(); i++) {
			PartitioningScheme * scheme = nextSchemes->at(i);
			for (unsigned int j = 0; j < scheme->getNumberOfElements(); j++) {
				PartitionElement * element = scheme->getElement(j);
				if (!(element->isOptimized() || element->isTagged())) {
					// wait for request
					element->setTagged(true);
					MPI_Recv(&buf, 1, MPI_INT, MPI_ANY_SOURCE, 0,
							MPI_COMM_WORLD, &targetStatus);
					buf = element->getId().size(); //getNumberOfSections();
					// send element
					MPI_Ssend(&buf, 1, MPI_INT, targetStatus.MPI_SOURCE, 1,
							MPI_COMM_WORLD);
					MPI_Ssend(&(element->getId().front()), buf, MPI_INT,
							targetStatus.MPI_SOURCE, 2, MPI_COMM_WORLD);
				}
			}
		}
		for (int i = 1; i < numProcs; i++) {
			MPI_Recv(&buf, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD,
					&targetStatus);
			buf = 0;
			MPI_Ssend(&buf, 1, MPI_INT, targetStatus.MPI_SOURCE, 1,
					MPI_COMM_WORLD);
		}
	}

	return 0;
}
#endif

int SearchAlgorithm::SchemeManager::optimize(ModelOptimize &mo) {
	t_partitionElementId id(3);
#ifdef _MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Status status;
	if (I_AM_ROOT) {
		pthread_t t1;
		pthread_create(&t1, NULL, &distribute, (void *) nextSchemes);
		int nextItem = 1;
		while (nextItem) {
			nextItem = 0;
			for (unsigned int i = 0; i < nextSchemes->size(); i++) {
				PartitioningScheme * scheme = nextSchemes->at(i);
				for (unsigned int j = 0; j < scheme->getNumberOfElements(); j++) {
					PartitionElement * element = scheme->getElement(j);
					if (!(element->isOptimized() || element->isTagged())) {
						element->setTagged(true);
						nextItem = 1;
						mo.optimizePartitionElement(element);
					}
				}
			}
		}
		pthread_join(t1, NULL);
	} else {
		int nextItem = 1;
		while (nextItem > 0) {
			MPI_Ssend(&nextItem, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
			MPI_Recv(&nextItem, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
			if (nextItem) {
				t_partitionElementId id(nextItem);
				MPI_Recv(&(id.front()), nextItem, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);
				PartitionElement * element = PartitionMap::getInstance()->getPartitionElement(id);
				mo.optimizePartitionElement(element);
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (I_AM_ROOT) {
		for (unsigned int i = 0; i < nextSchemes->size(); i++) {
			PartitioningScheme * scheme = nextSchemes->at(i);
			for (unsigned int j = 0; j < scheme->getNumberOfElements(); j++) {
				PartitionElement * element = scheme->getElement(j);
				if (element->isTagged() && !element->isOptimized()) {
					if (element->loadData()) {
						exit_partest(EX_IOERR);
					}
				}
			}
		}
	}
#else
	for (size_t i = 0; i < nextSchemes->size(); i++) {
		PartitioningScheme * scheme = nextSchemes->at(i);
		mo.optimizePartitioningScheme(scheme, i, nextSchemes->size());
	}
#endif
	nextSchemes->clear();
	return 0;
}

SearchAlgorithm::SearchAlgorithm() {
	if (starting_topology == StartTopoFIXED) {
		mo.buildStartingTree();
	}
}

SearchAlgorithm::~SearchAlgorithm() {

}

}
