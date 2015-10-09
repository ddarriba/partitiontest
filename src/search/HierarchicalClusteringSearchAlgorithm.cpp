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
 * @file HierarchicalClusteringSearch.cpp
 * @author Diego Darriba
 */

#include "HierarchicalClusteringSearchAlgorithm.h"

#include "util/GlobalDefs.h"
#include "util/Utilities.h"
#include "exe/ModelOptimize.h"
#include "exe/PartitionSelector.h"
#include "indata/PartitioningScheme.h"
#include "indata/PartitionMap.h"

#include <iostream>
#include <cmath>

using namespace std;

namespace partest {

HierarchicalClusteringSearchAlgorithm::HierarchicalClusteringSearchAlgorithm() {
}

HierarchicalClusteringSearchAlgorithm::~HierarchicalClusteringSearchAlgorithm() {
}

PartitioningScheme * HierarchicalClusteringSearchAlgorithm::start(
		PartitioningScheme * startingPoint) {

	SchemeManager schemeManager;
	vector<PartitioningScheme *> nextSchemes;

	PartitioningScheme *bestScheme = 0, *localBestScheme = 0;
	double bestScore, score;

	ModelOptimize * modelOptimize = new ModelOptimize();
	size_t numberOfPartitions = number_of_genes;

	bool continueExec = true;
	size_t maxSteps;
	if (I_AM_ROOT) {
		if (startingPoint) {
			nextSchemes.push_back(startingPoint);
			maxSteps = startingPoint->getNumberOfElements();
		} else {
			/* building first scheme */
			t_partitioningScheme * firstSchemeId = new t_partitioningScheme(
					number_of_genes);
			for (size_t gene = 0; gene < number_of_genes; gene++) {
				t_partitionElementId geneId(1);
				geneId.at(0) = gene;
				firstSchemeId->at(gene) = geneId;
			}
			nextSchemes.push_back(new PartitioningScheme(firstSchemeId));
			maxSteps = firstSchemeId->size();
			delete firstSchemeId;
		}
		bestScore = DOUBLE_INF;

		int currentStep = 0;

		while (continueExec) {

			cout << timestamp() << " [HCL] Step " << ++currentStep << "/"
					<< maxSteps << endl;

			for (size_t i=0; i<nextSchemes.size(); i++) {
				PartitioningScheme * scheme = nextSchemes.at(i);
				if (!scheme->isOptimized())
					schemeManager.addScheme(scheme);
			}
			schemeManager.optimize(mo);

			PartitionSelector ps(nextSchemes);
			localBestScheme = ps.getBestScheme();
			score = ps.getBestScheme()->getIcValue();

			printStepLog(currentStep, localBestScheme);

			if (score < bestScore) {
				//bestScheme->print(cout);
				delete bestScheme;
				bestScheme = localBestScheme;
				PartitionMap::getInstance()->keep(bestScheme->getId());
				for (size_t i = 0;
						i < localBestScheme->getNumberOfElements(); i++) {
					PartitionMap::getInstance()->purgePartitionMap(
							localBestScheme->getElement(i)->getId());
				}
				bestScore = score;
			}
			printStep(SearchHCluster, score);

			continueExec = ((non_stop || fabs(bestScore - score) < 1e-10)
					&& (numberOfPartitions > 1));
#ifdef HAVE_MPI
			MPI_Bcast(&continueExec, 1, MPI_INT, 0, MPI_COMM_WORLD );
#endif
			for (size_t i=0; i<nextSchemes.size(); i++) {
				PartitioningScheme * scheme = nextSchemes.at(i);
				if (scheme != localBestScheme)
					delete scheme;
			}
			nextSchemes.clear();

			if (continueExec) {
				vector<elementPair *> * eps =
						localBestScheme->getElementDistances();
				int sample_size;
				if (max_samples) {
					sample_size = min(max_samples, (int) eps->size());
				} else {
					sample_size = max((int) (eps->size() * samples_percent), 1);
				}
				for (size_t i = 0; i < (size_t)sample_size; i++) {
					t_partitionElementId nextId;
					t_partitioningScheme nextScheme;
					Utilities::mergeIds(nextId, eps->at(i)->e1->getId(),
							eps->at(i)->e2->getId());
					nextScheme.push_back(nextId);
					for (size_t j = 0;
							j < localBestScheme->getNumberOfElements(); j++) {

						PartitionElement * element =
								localBestScheme->getElement(j);
						if (element->getId() != eps->at(i)->e1->getId()
								&& element->getId()
										!= eps->at(i)->e2->getId()) {
							nextScheme.push_back(element->getId());
						}
					}
					nextSchemes.push_back(new PartitioningScheme(&nextScheme));
				}
				numberOfPartitions = nextSchemes.at(0)->getNumberOfElements();
			}
		}
	}
#ifdef HAVE_MPI
	else {
		while(continueExec) {
			schemeManager.optimize(mo);
			MPI_Bcast(&continueExec, 1, MPI_INT, 0, MPI_COMM_WORLD );
		}
	}
#endif

	if (bestScheme != localBestScheme)
		delete localBestScheme;

	delete modelOptimize;
	return bestScheme;
}

} /* namespace partest */
