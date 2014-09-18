/*
 * PLLTreeManager.h
 *
 *  Created on: Sep 18, 2014
 *      Author: diego
 */

#include "util/GlobalDefs.h"
#include "indata/TreeManager.h"

#include <pll.h>
#include <vector>

#ifndef PLLTREEMANAGER_H_
#define PLLTREEMANAGER_H_

using namespace std;

namespace partest {

class PllTreeManager : public TreeManager {
public:
	PllTreeManager(const pllAlignmentData * phylip,
			const vector<PEsection> & sections, size_t numberOfSites);
	virtual ~PllTreeManager();

	virtual double * getBranchLengths(void);

	pllInstance * _tree;
	pllAlignmentData * _alignData;
	partitionList * _partitions;

};

}

#endif /* PLLTREEMANAGER_H_ */
