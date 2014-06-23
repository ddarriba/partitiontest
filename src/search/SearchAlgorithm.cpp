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
//#include <mpi.h>

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

void * distribute(void * arg) {

	vector<PartitioningScheme *> * nextSchemes =
			(vector<PartitioningScheme *> *) arg;
	//vector<int *> * vec = (vector<int *> *)arg;
	int buf;
	MPI_Status targetStatus, statusRecv;

	for (unsigned int i = 0; i < nextSchemes->size(); i++) {
		PartitioningScheme * scheme = nextSchemes->at(i);
		for (unsigned int j = 0; j < scheme->getNumberOfElements(); j++) {
			PartitionElement * element = scheme->getElement(j);
			if (!element->isOptimized()) {
				// wait for request
//				cerr << "MASTER: wait for request" << endl;
				MPI_Recv(&buf, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,
						MPI_COMM_WORLD, &targetStatus);
				MPI_Barrier (MPI_COMM_WORLD);
				buf = element->getId().size(); //getNumberOfSections();
				// send element
//				cerr << "MASTER: " << buf << " sending request" << endl;
				MPI_Ssend(&buf, 1, MPI_INT, targetStatus.MPI_SOURCE, 1,
						MPI_COMM_WORLD);
				MPI_Barrier(MPI_COMM_WORLD);
//				cerr << "MASTER: " << buf << " sent request" << endl;
				MPI_Ssend(&(element->getId().front()), buf, MPI_INT,
						targetStatus.MPI_SOURCE, 2, MPI_COMM_WORLD);
//				cerr << "MASTER: " << element->getId().size() << " sent objects" << endl;
			}
		}
	}
	cerr << "FINALIZING FOR " << numProcs << endl;
	for (int i = 0; i < numProcs; i++) {
		// wait for request
		cout << "Going to finish!" << endl;
		MPI_Recv(&buf, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
				&targetStatus);
		// send termination
		cout << "Send FINALIZING!" << endl;
		buf = 0;
		MPI_Ssend(&buf, 1, MPI_INT, targetStatus.MPI_SOURCE, 1, MPI_COMM_WORLD);
	}
	return 0;
}

int SearchAlgorithm::SchemeManager::optimize(ModelOptimize &mo) {
	t_partitionElementId id(3);
#ifdef _MPI
	MPI_Status status;
	if (I_AM_ROOT) {
		pthread_t t1;
		cout << "LAUNCH NEW THREAD" << endl;
		pthread_create(&t1, NULL, &distribute, (void *) nextSchemes);

	} else {
		int nextItem = 1;
		while (nextItem > 0) {
			//			cerr << "SLAVE: send request" << endl;
			MPI_Ssend(&nextItem, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
			//			cerr << "SLAVE: receiving request" << endl;
			MPI_Recv(&nextItem, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
			//			cerr << "SLAVE: " << nextItem << " received request" << endl;
			if (nextItem) {
				t_partitionElementId id(nextItem);
				MPI_Recv(&(id.front()), nextItem, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);
				//				cerr << "SLAVE: " << id.size() << " received objects" << endl;
				PartitionElement * element = PartitionMap::getInstance()->getPartitionElement(id);
				mo.optimizePartitionElement(element);
			} else {
				cout << "FINALIZING!" << endl;
			}
		}
		nextSchemes->clear();
//#ifdef _MPI
//				int numElemsInScheme = scheme->getId().size();
//				MPI_Ssend(&numElemsInScheme, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
//				for (int i=0; i<numElemsInScheme; i++) {
//					int vecSize = scheme->getId().at(i).size();
//					MPI_Ssend(&vecSize, 1, MPI_INT, 1, i+1, MPI_COMM_WORLD);
//					MPI_Ssend(&(scheme->getId().at(i).front()), scheme->getId().at(i).size(), MPI_INT, 1, i+1, MPI_COMM_WORLD);
//				}
//#endif
//#ifdef _MPI
//			int zeroSend = 0;
//			MPI_Ssend(&zeroSend, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
//#endif
	}
#else
	for (int i = 0; i < nextSchemes.size(); i++) {
		PartitioningScheme * scheme = nextSchemes.at(i);
		mo.optimizePartitioningScheme(scheme, i, nextSchemes.size());
	}
#endif
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
