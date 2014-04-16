/*
 * SearchAlgorithm.cpp
 *
 *  Created on: Apr 16, 2014
 *      Author: diego
 */

#include "SearchAlgorithm.h"

#include <iostream>

namespace partest {

SearchAlgorithm::SearchAlgorithm() {
	if (starting_topology == StartTopoFIXED) {
		mo.buildStartingTree();
	}
}

SearchAlgorithm::~SearchAlgorithm() {

}

}
