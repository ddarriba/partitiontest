/*
 * TreeManager.cpp
 *
 *  Created on: Sep 18, 2014
 *      Author: diego
 */

#include "TreeManager.h"

namespace partest {

TreeManager::TreeManager(const t_partitionElementId id, size_t numberOfSites, size_t numberOfPatterns) :
	_id(id), numberOfSites(numberOfSites), numberOfPatterns(numberOfPatterns), branchLengthMultiplier(1.0) {
}

TreeManager::~TreeManager() {
}

} /* namespace partest */
