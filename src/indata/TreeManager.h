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

	size_t getNumberOfTaxa(void) {
			return numberOfTaxa;
		}

	size_t getNumberOfPatterns(void) {
		return numberOfPatterns;
	}

	size_t getNumberOfSites(void) {
		return numberOfSites;
	}

	/**
	 * Allocates and returns a vector with the branch lengths
	 */
	virtual double * getBranchLengths(bool update = true) = 0;
	virtual void setBranchLengths(double * bls) = 0;

	/**
	 * Applies the parameters contained in a Model instance to the tree
	 *
	 * @param model The model to take the parameters
	 * @param index The partition index
	 * @param setAlphaFreqs Whether to initialize alpha and frequencies
	 */
	virtual void setModelParameters(const Model * model, int index,	bool setAlphaFreqs) = 0;

	/**
	 * Performs a Maximum-Likelihood Search
	 *
	 * @param estimateModel Whether to also optimize model parameters
	 */
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

	double getBranchLengthMultiplier( void ) {
		return branchLengthMultiplier;
	}

protected:
	t_partitionElementId _id;
	size_t numberOfTaxa;
	size_t numberOfSites;
	size_t numberOfPatterns;
	double * branchLengths;
	double branchLengthMultiplier;
};

} /* namespace partest */

#endif /* TREEMANAGER_H_ */
