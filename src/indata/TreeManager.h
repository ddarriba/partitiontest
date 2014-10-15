/*
 * TreeManager.h
 *
 *  Created on: Sep 18, 2014
 *      Author: diego
 */

#ifndef TREEMANAGER_H_
#define TREEMANAGER_H_

#include "model/Model.h"
#include <cstdlib>

namespace partest {

class TreeManager {
public:

	TreeManager(const t_partitionElementId _id, size_t numberOfSites, size_t numberOfPatterns);
	virtual ~TreeManager();

	size_t getNumberOfPatterns(void) {
		return numberOfPatterns;
	}

	size_t getNumberOfSites(void) {
		return numberOfSites;
	}

	/**
	 * Allocates and returns a vector with the branch lengths
	 */
	virtual double * getBranchLengths(void) = 0;

	virtual void setModelParameters(Model * _model, int index,	bool setAlphaFreqs) = 0;
	virtual double searchMlTopology(bool estimateModel) = 0;
	virtual double getLikelihood() = 0;
	virtual void optimizeBranchLengths(int smoothIterations) = 0;
	virtual void optimizeModelParameters(double epsilon) = 0;
	virtual void optimizeBaseFreqs(double epsilon) = 0;
	virtual void optimizeRates(double epsilon) = 0;
	virtual void optimizeAlphas(double epsilon) = 0;
	virtual double evaluateLikelihood(bool fullTraversal) = 0;
	virtual const char * getNewickTree() = 0;

	virtual double * getFrequencies(size_t partition = 0) = 0;
	virtual double * getRates(size_t partition = 0) = 0;
	virtual double getAlpha(size_t partition = 0) = 0;

	virtual int getAutoProtModel(size_t partition = 0) = 0;

protected:
	t_partitionElementId _id;
	size_t numberOfSites;
	size_t numberOfPatterns;
};

} /* namespace partest */

#endif /* TREEMANAGER_H_ */
