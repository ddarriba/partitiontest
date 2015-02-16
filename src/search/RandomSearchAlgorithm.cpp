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
 * @author Diego Darriba
 */

#include "RandomSearchAlgorithm.h"

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

RandomSearchAlgorithm::RandomSearchAlgorithm() {
}

RandomSearchAlgorithm::~RandomSearchAlgorithm() {
}

PartitioningScheme * RandomSearchAlgorithm::start(
		PartitioningScheme * startingPoint) {

	SchemeManager schemeManager;
	vector<PartitioningScheme *> nextSchemes;

	PartitioningScheme *bestScheme = 0, *localBestScheme = 0;
	double bestScore, score;

	ModelOptimize * modelOptimize = new ModelOptimize();

	bool continueExec = true;
	size_t maxSteps;
	if (I_AM_ROOT) {
		if (startingPoint) {
			nextSchemes.push_back(startingPoint);
			maxSteps = startingPoint->getId().size();
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

			cout << timestamp() << " [RND] Step " << ++currentStep << "/"
					<< maxSteps << endl;

			for (size_t i = 0; i < nextSchemes.size(); i++) {
				PartitioningScheme * scheme = nextSchemes.at(i);
				if (!scheme->isOptimized())
					schemeManager.addScheme(scheme);
			}
			schemeManager.optimize(mo);

			PartitionSelector ps(nextSchemes);
			localBestScheme = ps.getBestScheme();
			score = ps.getBestScheme()->getIcValue();

			if (score < bestScore) {
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
			printStep(SearchRandom, score);

			continueExec = ((non_stop || fabs(bestScore - score) < 1e-10)
					&& (localBestScheme->getNumberOfElements() > 1));
#ifdef HAVE_MPI
			MPI_Bcast(&continueExec, 1, MPI_INT, 0, MPI_COMM_WORLD );
#endif
			for (size_t i = 0; i < nextSchemes.size(); i++) {
				PartitioningScheme * scheme = nextSchemes.at(i);
				if (scheme != localBestScheme)
					delete scheme;
			}
			nextSchemes.clear();

			if (continueExec) {
				getRandomPartitioningScheme(nextSchemes, max_samples,
						localBestScheme->getId());
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

int RandomSearchAlgorithm::getRandomPartitioningScheme(
		vector<PartitioningScheme *> & nextSchemes, int numberOfSchemes,
		t_partitioningScheme p0) {

	size_t maxClasses = p0.size();

	for (int schemeIndex = 0; schemeIndex < numberOfSchemes; schemeIndex++) {
		t_partitionElementId * classes = (t_partitionElementId *) malloc(maxClasses * sizeof(t_partitionElementId));

		int numberOfClasses = 1;
		// first element to the first
		for (size_t k=0; k<p0.at(0).size(); k++) {
			size_t singleP = p0.at(0).at(k);
			classes[0].push_back(singleP);
		}
		for (size_t i = 1; i < maxClasses; i++) {
			int currentClass = 0;

			//randomly assign each bit to each class
			int rndNumber = rand() % numberOfClasses;
			currentClass = rndNumber;
			if (currentClass < numberOfClasses) {
				for (size_t k=0; k<p0.at(i).size(); k++) {
					size_t singleP = p0.at(i).at(k);
					classes[(size_t)currentClass].push_back(singleP);
				}
			} else {
				currentClass = numberOfClasses;
				for (size_t k=0; k<p0.at(i).size(); k++) {
					size_t singleP = p0.at(i).at(k);
					classes[currentClass].push_back(singleP);
				}
				numberOfClasses++;
			}

			free(classes);
		}

		t_partitioningScheme * newSchemeId = new t_partitioningScheme(
				(size_t) numberOfClasses);
		for (size_t i = 0; i < (size_t) numberOfClasses; i++) {
			for (size_t j = 0; j < classes[i].size(); j++) {
				newSchemeId->at(i).push_back(classes[i].at(j));
			}
		}
		nextSchemes.push_back(new PartitioningScheme(newSchemeId));
		delete newSchemeId;
	}

	return 0;
}

} /* namespace partest */
