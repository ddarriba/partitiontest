/*
 * TreeManager.cpp
 *
 *  Created on: Sep 18, 2014
 *      Author: diego
 */

#include "TreeManager.h"

namespace partest {

TreeManager::TreeManager(const t_partitionElementId id, size_t _numberOfSites, size_t _numberOfPatterns) :
	_id(id), numberOfSites(_numberOfSites), numberOfPatterns(_numberOfPatterns), branchLengthMultiplier(1.0) {
}

TreeManager::~TreeManager() {
}

} /* namespace partest */
