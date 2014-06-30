/*
 * GreedySearchAlgorithm.cpp
 *
 *  Created on: Apr 11, 2014
 *      Author: diego
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
	nextSchemeFunctor(PartitioningScheme * scheme) :
			elementIndex1(1), elementIndex2(0), scheme(scheme) {
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
	int size() {
		return ((numberOfElements * (numberOfElements - 1)) / 2);
	}
private:
	size_t elementIndex1, elementIndex2;
	size_t numberOfElements, numberOfSchemes, currentScheme;
	PartitioningScheme * scheme;
};

GreedySearchAlgorithm::GreedySearchAlgorithm() {
	// TODO Auto-generated constructor stub

}

GreedySearchAlgorithm::~GreedySearchAlgorithm() {
	// TODO Auto-generated destructor stub
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
	int numberOfPartitions = number_of_genes;
	vector<PartitioningScheme *> nextSchemes;
	int step = 1;
	int maxSteps = number_of_genes;

	bool continueExec = (numberOfPartitions > 1);
	if (I_AM_ROOT) {

		/* building first scheme */
		cout << timestamp() << " [GRE] Step " << step++ << "/" << maxSteps
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
#ifdef _MPI
		MPI_Bcast(&continueExec, 1, MPI_INT, 0, MPI_COMM_WORLD );
#endif

		//mo.optimizePartitioningScheme(bestScheme);
		PartitionSelector ps(nextSchemes);
		nextSchemes.clear();

		bestScore = score = ps.getBestScheme()->getIcValue();

		while (continueExec) {
			cout << timestamp() << " [GRE] Step " << step++ << "/" << maxSteps
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

			PartitionSelector ps(nextSchemes);
			//ps.print(cout);
			localBestScheme = ps.getBestScheme();
			score = localBestScheme->getIcValue();

			if (score < bestScore) {
				delete bestScheme;
				bestScheme = localBestScheme;
				PartitionMap::getInstance()->keep(bestScheme->getId());
				cout << timestamp() << " [GRE] Improving " << bestScore - score
						<< " score units." << endl;
				bestScore = score;
			} else {
				cout << timestamp() << " [GRE] Scheme is " << score - bestScore
						<< " score units ahead the best score." << endl;
			}

			continueExec = ((non_stop || bestScore == score)
					&& (numberOfPartitions > 1));

#ifdef _MPI
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

#ifdef _MPI
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
