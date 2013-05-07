/*
 * SearchAlgorithm.cpp
 *
 *  Created on: Mar 8, 2013
 *      Author: diego
 */

#include "SearchAlgorithm.h"
#include "util/PrintMeta.h"
#include <iostream>
#include <fstream>

namespace partest {

SearchAlgorithm::SearchAlgorithm(ParTestOptions * options,
		PartitionMap * partitionMap) :
		options(options), partitionMap(partitionMap) {
	ofstream modelsOutputStream;
	modelsOutputStream.open(options->getOutputFileModels().c_str());
	PrintMeta::print_header(modelsOutputStream);
	PrintMeta::print_options(modelsOutputStream, *options);
	modelsOutputStream.close();
}

SearchAlgorithm::~SearchAlgorithm() {
}

PartitionElement * SearchAlgorithm::getPartitionElement(
		t_partitionElementId id) {
	return partitionMap->getPartitionElement(id);
}

unsigned int SearchAlgorithm::getNumberOfElements() {
	return partitionMap->getNumberOfElements();
}
} /* namespace partest */
