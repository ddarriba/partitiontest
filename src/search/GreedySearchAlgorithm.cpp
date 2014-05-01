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
			scheme(scheme) {
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

		t_partitionElementId e1 = scheme->getElement(i)->getId();
		t_partitionElementId e2 = scheme->getElement(j)->getId();
		t_partitionElementId nextId;
		Utilities::mergeIds(nextId, e1, e2);
		nextSchemeId.push_back(nextId);

		/* fill */
		for (unsigned int k = 0; k < numberOfElements; k++) {

			t_partitionElementId eX = scheme->getElement(k)->getId();
			if (eX != e1 && eX != e2) {
				nextSchemeId.push_back(eX);
			}
		}

		currentScheme++;
		j++;
		if (j == i) {
			i++; j=0;
		}

		return (new PartitioningScheme(&nextSchemeId));
	}
	int size() {
		return ((numberOfElements * (numberOfElements - 1)) / 2);
	}
private:
	unsigned int i = 1, j = 0;
	unsigned int numberOfElements, numberOfSchemes, currentScheme;
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
	for (unsigned int i = 1; i < startingScheme->size(); i++) {
		t_partitionElementId e1 = startingScheme->at(i);
		for (unsigned int j = 0; j < i; j++) {

			t_partitioningScheme nextSchemeId(startingScheme->size() - 1);

			/* merge elements i and j */
			t_partitionElementId e2 = startingScheme->at(j);
			t_partitionElementId nextId;
			Utilities::mergeIds(nextId, e1, e2);

			/* fill */
			for (unsigned int k = 0; k < startingScheme->size(); k++) {

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

PartitioningScheme * GreedySearchAlgorithm::start() {

	PartitioningScheme *bestScheme = 0, *localBestScheme = 0;
	double bestScore, score;
	ModelOptimize * modelOptimize = new ModelOptimize();
	int numberOfPartitions = number_of_genes;
	vector<PartitioningScheme *> nextSchemes;
	int step = 1;
	int maxSteps = number_of_genes;

	/* building first scheme */
	cout << timestamp() << " [GRE] Step " << step++ << "/" << maxSteps << endl;
	t_partitioningScheme * firstSchemeId = new t_partitioningScheme(
			number_of_genes);
	for (unsigned int gene = 0; gene < number_of_genes; gene++) {
		t_partitionElementId geneId(1);
		geneId.at(0) = gene;
		firstSchemeId->at(gene) = geneId;
	}

	localBestScheme = bestScheme = new PartitioningScheme(firstSchemeId);
	nextSchemes.push_back(bestScheme);
	delete firstSchemeId;

	modelOptimize->optimizePartitioningScheme(bestScheme);

	PartitionSelector ps(nextSchemes);
	nextSchemes.clear();

	bestScore = score = ps.getBestScheme()->getIcValue();

	bool improving = true;
	while (improving) {
		cout << timestamp() << " [GRE] Step " << step++ << "/" << maxSteps << endl;
		nextSchemeFunctor nextScheme(localBestScheme);
		nextSchemes.reserve(nextScheme.size());

		int schemeIndex = 0;
		while (PartitioningScheme * scheme = nextScheme()) {
			modelOptimize->optimizePartitioningScheme(scheme, schemeIndex++, nextScheme.size());
			nextSchemes.push_back(scheme);
		}
		numberOfPartitions = localBestScheme->getNumberOfElements()-1;

		PartitionSelector ps(nextSchemes);
		//ps.print(cout);
		localBestScheme = ps.getBestScheme();
		score = localBestScheme->getIcValue();

		if (score < bestScore) {
			delete bestScheme;
			bestScheme = localBestScheme;
			PartitionMap::getInstance()->keep(bestScheme->getId());
			cout << timestamp() << " [GRE] Improving " << bestScore - score << " score units." << endl;
			bestScore = score;
		} else {
			cout << timestamp() << " [GRE] Scheme is " << score - bestScore << " score units ahead the best score." << endl;
		}
		improving = ((non_stop || bestScore == score) && (numberOfPartitions > 1));

		for (PartitioningScheme * scheme : nextSchemes) {
			if (scheme != localBestScheme) {
				delete scheme;
			}
		}

		if (improving) {
			for (unsigned int i = 0; i < localBestScheme->getNumberOfElements(); i++) {
				PartitionMap::getInstance()->purgePartitionMap(
						localBestScheme->getElement(i)->getId());
			}
		}

		nextSchemes.clear();
	}

	return bestScheme;
}

} /* namespace partest */
