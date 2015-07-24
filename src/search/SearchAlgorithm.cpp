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
 * @file SearchAlgorithm.cpp
 * @author Diego Darriba
 */

#include "SearchAlgorithm.h"

#include "indata/PartitionMap.h"
#include <iostream>
#include <iomanip>
#include <pthread.h>
#include <memory>
#include <cmath>
#include <unistd.h>

using namespace std;

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
	return (int) nextSchemes->size();
}

int SearchAlgorithm::SchemeManager::addScheme(
		PartitioningScheme * schemeToAdd) {
	nextSchemes->push_back(schemeToAdd);
	return (int) nextSchemes->size();
}

void SearchAlgorithm::printStepLog(int id, PartitioningScheme *bestScheme) {
	(*ofs) << id << "\t" << fixed << setprecision(4) << bestScheme->getLnL() << "\t"
			<< bestScheme->getNumberOfFreeParameters() << "\t"
			<< fixed << setprecision(4) << bestScheme->getBicValue() << "\t"
			<< fixed << setprecision(4) << bestScheme->getAicValue() << "\t"
			<< fixed << setprecision(4) << bestScheme->getAiccValue() << "\t"
			<< fixed << setprecision(4) << bestScheme->getLinkedBicValue() << "\t"
			<< fixed << setprecision(4) << bestScheme->getLinkedAicValue() << "\t"
			<< fixed << setprecision(4) << bestScheme->getLinkedAiccValue() << endl;
}

void SearchAlgorithm::printStep(SearchAlgo algo, double nextScore) {
	cout << timestamp();
	switch (algo) {
	case SearchK1:
		cout << " [K=1] ";
		break;
	case SearchKN:
		cout << " [K=N] ";
		break;
	case SearchHCluster:
		cout << " [HCL] ";
		break;
	case SearchGreedy:
	case SearchGreedyExtended:
		cout << " [GRE] ";
		break;
	case SearchExhaustive:
		cout << " [EXH] ";
		break;
	case SearchRandom:
		cout << " [RND] ";
		break;
	case SearchAuto:
		cout << " [XXX] ";
		break;
	}
	if (fabs(bestScore) < 1e-10) {
		bestScore = nextScore;
		cout << "Initial score is " << fixed << setprecision(4) << nextScore
				<< "." << endl;
	} else {
		if (nextScore < bestScore) {
			cout << "Improving " << fixed << setprecision(4)
					<< bestScore - nextScore << " score units (" << fixed
					<< setprecision(4) << nextScore << ")." << endl;
			bestScore = nextScore;
		} else {
			cout << "Scheme is " << fixed << setprecision(4)
					<< nextScore - bestScore
					<< " score units ahead the best score (" << fixed
					<< setprecision(4) << nextScore << ")." << endl;
		}
	}
}

#ifdef HAVE_MPI
void * distribute(void * arg) {
	vector<PartitioningScheme *> * nextSchemes =
			(vector<PartitioningScheme *> *) arg;

	if (numProcs > 1) {
		int buf[3];
		MPI_Status targetStatus;
		for (size_t i = 0; i < nextSchemes->size(); i++) {
			PartitioningScheme * scheme = nextSchemes->at(i);
			buf[2] = scheme->getNumberOfElements();
			for (size_t j = 0; j < scheme->getNumberOfElements(); j++) {
				PartitionElement * element = scheme->getElement(j);
				if (!(element->isOptimized() || element->isTagged())) {
					// wait for request
					element->setTagged(true);
					MPI_Recv(buf, 1, MPI_INT, MPI_ANY_SOURCE, 0,
							MPI_COMM_WORLD, &targetStatus);
					buf[0] = element->getId().size(); //getNumberOfSections();
					buf[1] = j;
					// send element
					MPI_Ssend(buf, 3, MPI_INT, targetStatus.MPI_SOURCE, 1,
							MPI_COMM_WORLD);
					MPI_Ssend(&(element->getId().front()), buf[0], MPI_INT,
							targetStatus.MPI_SOURCE, 2, MPI_COMM_WORLD);
				}
			}
		}
		for (int i = 1; i < numProcs; i++) {
			MPI_Recv(buf, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD,
					&targetStatus);
			buf[0] = 0;
			MPI_Ssend(buf, 3, MPI_INT, targetStatus.MPI_SOURCE, 1,
					MPI_COMM_WORLD);
		}
	}

	return 0;
}
#endif

int SearchAlgorithm::SchemeManager::optimize(ModelOptimize &_mo) {
	t_partitionElementId id(3);
#ifdef HAVE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Status status;
	if (I_AM_ROOT) {
		pthread_t t1;
		pthread_create(&t1, NULL, &distribute, (void *) nextSchemes);
		int nextItem = 1;
		while (nextItem) {
			nextItem = 0;
			for (size_t i = 0; i < nextSchemes->size(); i++) {
				PartitioningScheme * scheme = nextSchemes->at(i);
				size_t numElements = scheme->getNumberOfElements();
				for (size_t j = 0; j < numElements; j++) {
					PartitionElement * element = scheme->getElement(j);
					if (!(element->isOptimized() || element->isTagged())) {
						element->setTagged(true);
						nextItem = 1;
						_mo.optimizePartitionElement(element, j, numElements);
					}
				}
			}
		}
		pthread_join(t1, NULL);
	} else {
		int nextItem[3];
		nextItem[0] = 1;
		while (nextItem[0] > 0) {
			MPI_Ssend(nextItem, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
			MPI_Recv(nextItem, 3, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
			if (nextItem[0]) {
				t_partitionElementId id(nextItem[0]);
				MPI_Recv(&(id.front()), nextItem[0], MPI_INT, 0, 2, MPI_COMM_WORLD, &status);
				PartitionElement * element = PartitionMap::getInstance()->getPartitionElement(id);
				_mo.optimizePartitionElement(element, nextItem[1], nextItem[2]);
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (I_AM_ROOT) {
		for (size_t i = 0; i < nextSchemes->size(); i++) {
			PartitioningScheme * scheme = nextSchemes->at(i);
			for (size_t j = 0; j < scheme->getNumberOfElements(); j++) {
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
		_mo.optimizePartitioningScheme(scheme, (int) i, (int) nextSchemes->size());
	}
#endif
	nextSchemes->clear();
	return 0;
}

SearchAlgorithm::SearchAlgorithm() {
	if (starting_topology == StartTopoFIXED) {
		mo.buildStartingTree();
	}
	if (outputAvailable && log_logfile) {
		ofs = new ofstream(log_logfile->c_str(), ios::in | ios::out | ios::app);
	} else {
		ofs = 0;
	}
	bestScore = 0.0;
}

SearchAlgorithm::~SearchAlgorithm() {
	if (ofs) {
		ofs->close();
		delete ofs;
	}
}

}
