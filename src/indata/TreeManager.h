/*
 * TreeManager.h
 *
 *  Created on: Sep 18, 2014
 *      Author: diego
 */

#ifndef TREEMANAGER_H_
#define TREEMANAGER_H_

#include <cstdlib>

namespace partest {

class TreeManager {
public:

	TreeManager(size_t numberOfSites, size_t numberOfPatterns);
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

protected:
	size_t numberOfSites;
	size_t numberOfPatterns;
};

} /* namespace partest */

#endif /* TREEMANAGER_H_ */
