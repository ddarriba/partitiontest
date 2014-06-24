/*
 * HierarchicalClusteringSearchAlgorithm.cpp
 *
 *  Created on: Apr 9, 2014
 *      Author: diego
 */

#include "HierarchicalClusteringSearchAlgorithm.h"

#include "util/GlobalDefs.h"
#include "util/Utilities.h"
#include "exe/ModelOptimize.h"
#include "exe/PartitionSelector.h"
#include "indata/PartitioningScheme.h"
#include "indata/PartitionMap.h"

#include <iostream>

using namespace std;

namespace partest {

HierarchicalClusteringSearchAlgorithm::HierarchicalClusteringSearchAlgorithm() {
}

HierarchicalClusteringSearchAlgorithm::~HierarchicalClusteringSearchAlgorithm() {
}

PartitioningScheme * HierarchicalClusteringSearchAlgorithm::start() {

	SchemeManager schemeManager;
	vector<PartitioningScheme *> nextSchemes;

	PartitioningScheme *bestScheme = 0, *localBestScheme = 0;
	double bestScore, score;

	ModelOptimize * modelOptimize = new ModelOptimize();
	int numberOfPartitions = number_of_genes;

	bool continueExec = true;
	if (I_AM_ROOT) {
		/* building first scheme */
		t_partitioningScheme * firstSchemeId = new t_partitioningScheme(
				number_of_genes);
		for (unsigned int gene = 0; gene < number_of_genes; gene++) {
			t_partitionElementId geneId(1);
			geneId.at(0) = gene;
			firstSchemeId->at(gene) = geneId;
		}
		nextSchemes.push_back(new PartitioningScheme(firstSchemeId));

		bestScore = DOUBLE_INF;

		int currentStep = 1;
		int maxSteps = firstSchemeId->size();
		delete firstSchemeId;

		while (continueExec) {

			cout << timestamp() << " [HCL] Step " << currentStep++ << "/"
					<< maxSteps << endl;

			for (PartitioningScheme * scheme : nextSchemes) {
				if (!scheme->isOptimized())
					schemeManager.addScheme(scheme);
			}
			schemeManager.optimize(mo);

			PartitionSelector ps(nextSchemes);
			localBestScheme = ps.getBestScheme();
			score = ps.getBestScheme()->getIcValue();

			if (score < bestScore) {
				//bestScheme->print(cout);
				delete bestScheme;
				bestScheme = localBestScheme;
				PartitionMap::getInstance()->keep(bestScheme->getId());
				for (unsigned int i = 0;
						i < localBestScheme->getNumberOfElements(); i++) {
					PartitionMap::getInstance()->purgePartitionMap(
							localBestScheme->getElement(i)->getId());
				}
				if (currentStep > 1) {
					cout << timestamp() << " [HCL] Improving "
							<< bestScore - score << " score units." << endl;
				}
				bestScore = score;
			} else {
				cout << timestamp() << " [HCL] Scheme is " << score - bestScore
						<< " score units ahead the best score." << endl;
			}

			continueExec = ((non_stop || bestScore == score)
					&& (numberOfPartitions > 1));
#ifdef _MPI
			MPI_Bcast(&continueExec, 1, MPI_INT, 0, MPI_COMM_WORLD );
#endif
			for (PartitioningScheme * scheme : nextSchemes) {
				if (scheme != localBestScheme)
					delete scheme;
			}
			nextSchemes.clear();

			if (continueExec) {
				vector<elementPair *> * eps =
						localBestScheme->getElementDistances();
				for (int i = 0; i < min(max_samples, (int) eps->size()); i++) {
					t_partitionElementId nextId;
					t_partitioningScheme nextScheme;
					Utilities::mergeIds(nextId, eps->at(i)->e1->getId(),
							eps->at(i)->e2->getId());
					nextScheme.push_back(nextId);
					for (unsigned int j = 0;
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
#ifdef _MPI
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
