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
 * @file GreedySearchAlgorithm.cpp
 * @author Diego Darriba
 */

#include "GreedySearchAlgorithm.h"

#include "util/Utilities.h"
#include "exe/ModelOptimize.h"
#include "indata/PartitionMap.h"
#include "exe/PartitionSelector.h"
#include <iostream>

using namespace std;

namespace partest {

struct nextSchemeFunctor {
	nextSchemeFunctor(PartitioningScheme * _scheme) :
			elementIndex1(1), elementIndex2(0), scheme(_scheme) {
		numberOfElements = scheme->getNumberOfElements();
		numberOfSchemes = (numberOfElements * (numberOfElements - 1)) / 2;
		currentScheme = 0;
	}

	~nextSchemeFunctor() {

	}
	PartitioningScheme * operator()(void) {

		if (currentScheme >= numberOfSchemes)
			return 0;

		t_partitioningScheme nextSchemeId;
		nextSchemeId.reserve(numberOfElements - 1);

		t_partitionElementId e1 = scheme->getElement(elementIndex1)->getId();
		t_partitionElementId e2 = scheme->getElement(elementIndex2)->getId();
		t_partitionElementId nextId;
		Utilities::mergeIds(nextId, e1, e2);
		nextSchemeId.push_back(nextId);

		/* fill */
		for (size_t k = 0; k < numberOfElements; k++) {

			t_partitionElementId eX = scheme->getElement(k)->getId();
			if (eX != e1 && eX != e2) {
				nextSchemeId.push_back(eX);
			}
		}

		currentScheme++;
		elementIndex2++;
		if (elementIndex2 == elementIndex1) {
			elementIndex1++;
			elementIndex2 = 0;
		}

		return (new PartitioningScheme(&nextSchemeId));
	}
	size_t size() {
		return ((numberOfElements * (numberOfElements - 1)) / 2);
	}
private:
	size_t elementIndex1, elementIndex2;
	size_t numberOfElements, numberOfSchemes, currentScheme;
	PartitioningScheme * scheme;
};

GreedySearchAlgorithm::GreedySearchAlgorithm() {
}

GreedySearchAlgorithm::~GreedySearchAlgorithm() {
}

vector<PartitioningScheme *> GreedySearchAlgorithm::getNextSchemes(
		const t_partitioningScheme * startingScheme) {

	vector<PartitioningScheme *> nextSchemes;
	for (size_t i = 1; i < startingScheme->size(); i++) {
		t_partitionElementId e1 = startingScheme->at(i);
		for (size_t j = 0; j < i; j++) {

			t_partitioningScheme nextSchemeId(startingScheme->size() - 1);

			/* merge elements i and j */
			t_partitionElementId e2 = startingScheme->at(j);
			t_partitionElementId nextId;
			Utilities::mergeIds(nextId, e1, e2);

			/* fill */
			for (size_t k = 0; k < startingScheme->size(); k++) {

				t_partitionElementId eX = startingScheme->at(k);
				if (eX != e1 && eX != e2) {
					nextSchemeId.push_back(eX);
				}
			}
			nextSchemes.push_back(new PartitioningScheme(&nextSchemeId));
		}
	}
	return nextSchemes;

}

PartitioningScheme * GreedySearchAlgorithm::start(
		PartitioningScheme * startingPoint) {

	SchemeManager schemeManager;

	PartitioningScheme *bestScheme = 0, *localBestScheme = 0;
	double bestScore, score;
	size_t numberOfPartitions = number_of_genes;
	vector<PartitioningScheme *> nextSchemes;
	int currentStep = 0;
	size_t maxSteps = number_of_genes;

	bool continueExec = (numberOfPartitions > 1);
	if (I_AM_ROOT) {

		/* building first scheme */
		cout << timestamp() << " [GRE] Step " << ++currentStep << "/" << maxSteps
				<< endl;
		if (startingPoint) {
			nextSchemes.push_back(startingPoint);
			maxSteps = startingPoint->getNumberOfElements();
			localBestScheme = bestScheme = startingPoint;
		} else {
			t_partitioningScheme * firstSchemeId = new t_partitioningScheme(
					number_of_genes);
			for (size_t gene = 0; gene < number_of_genes; gene++) {
				t_partitionElementId geneId(1);
				geneId.at(0) = gene;
				firstSchemeId->at(gene) = geneId;
			}

			localBestScheme = bestScheme = new PartitioningScheme(
					firstSchemeId);
			nextSchemes.push_back(bestScheme);
			delete firstSchemeId;
		}

		schemeManager.addScheme(bestScheme);
		schemeManager.optimize(mo);
#ifdef HAVE_MPI
		MPI_Bcast(&continueExec, 1, MPI_INT, 0, MPI_COMM_WORLD );
#endif

		printStepLog(currentStep, localBestScheme);

		//mo.optimizePartitioningScheme(bestScheme);
		PartitionSelector ps(nextSchemes);
		nextSchemes.clear();

		bestScore = score = ps.getBestScheme()->getIcValue();

		while (continueExec) {
			cout << timestamp() << " [GRE] Step " << ++currentStep << "/" << maxSteps
					<< endl;
			nextSchemeFunctor nextScheme(localBestScheme);
			nextSchemes.reserve(nextScheme.size());

			while (PartitioningScheme * scheme = nextScheme()) {
				if (!scheme->isOptimized())
					schemeManager.addScheme(scheme);
				nextSchemes.push_back(scheme);
			}
			schemeManager.optimize(mo);

			numberOfPartitions = localBestScheme->getNumberOfElements() - 1;

			PartitionSelector _ps(nextSchemes);
			//ps.print(cout);
			localBestScheme = _ps.getBestScheme();
			score = localBestScheme->getIcValue();

			printStepLog(currentStep, localBestScheme);

			if (score < bestScore) {
				delete bestScheme;
				bestScheme = localBestScheme;
				PartitionMap::getInstance()->keep(bestScheme->getId());
				bestScore = score;
			}
			printStep(SearchGreedy, score);

			continueExec = ((non_stop || fabs(bestScore - score) < 1e-10)
					&& (numberOfPartitions > 1));

#ifdef HAVE_MPI
			MPI_Bcast(&continueExec, 1, MPI_INT, 0, MPI_COMM_WORLD );
#endif
			for (size_t i=0; i<nextSchemes.size(); i++) {
				PartitioningScheme * scheme = nextSchemes.at(i);
				if (scheme != localBestScheme) {
					delete scheme;
				}
			}

			if (continueExec) {
				for (size_t i = 0;
						i < localBestScheme->getNumberOfElements(); i++) {
					PartitionMap::getInstance()->purgePartitionMap(
							localBestScheme->getElement(i)->getId());
				}
			}

			nextSchemes.clear();
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

	return bestScheme;
}

} /* namespace partest */
