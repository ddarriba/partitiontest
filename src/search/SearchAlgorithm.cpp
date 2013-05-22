/*  PartitionTest, fast selection of the best fit partitioning scheme for
 *  multi-gene data sets.
 *  Copyright May 2013 by Diego Darriba
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 3 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 *
 *  For any other inquiries send an Email to Diego Darriba
 *  ddarriba@udc.es
 */

/**
 * @file SearchAlgorithm.cpp
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
