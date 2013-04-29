/*
 * GenomicTest 1.0 source
 *
 * File: genomictest.cpp
 *
 * Copyright (c) 2012 Diego Darriba
 *
 * GenomicTest is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License (the file gpl.txt included with this
 * distribution or http://www.gnu.org/licenses/gpl.txt for more
 * details.
 */
#include <iostream>
#include "util/PrintMeta.h"
#include "util/GlobalDefs.h"
#include "util/ParTestFactory.h"
#include "options/ParTestOptions.h"
#include "parser/ArgumentParser.h"
#include "indata/PartitionMap.h"
#include "selection/ModelSelector.h"

using namespace std;
using namespace partest;

//namespace partest {

int main(int argc, char *argv[]) {

	PrintMeta::print_header (cout);

#ifdef DEBUG
	cout << "[TRACE] Creating argument parser" << endl;
#endif
	ArgumentParser* parser = new ArgumentParser();
#ifdef DEBUG
	cout << "[TRACE] Building options" << endl;
#endif
	ParTestOptions *options = new ParTestOptions();
#ifdef DEBUG
	cout << "[TRACE] Parsing application arguments" << endl;
#endif
	parser->fill_options(argc, argv, options);

	delete parser;

	PrintMeta::print_options(cout, *options);

//	PartitionMap * partitionMap = new PartitionMap(options->getAlignment(), 3,
//			options->getRateVariation(), options->getDataType());
//	partitionMap->addPartitionElement(0, 1, 50, 1);
//	partitionMap->addPartitionElement(1, 51, 114, 1);
//	partitionMap->addPartitionElement(2, 1, 114, 1);

	PartitionMap * partitionMap = new PartitionMap(options->getConfigFile(), options->getAlignment(),
				options->getRateVariation(), options->getDataType());

	SearchAlgorithm * searchAlgo = ParTestFactory::createSearchAlgorithm( options, partitionMap );

	PartitioningScheme * partition = searchAlgo->start();

	cout << "***** BEST PARTITIONING SCHEME *****" << endl;
	cout << "Partitioning scheme: " << partition->toString() << endl;
	cout << "Number of elements:  " << partition->getNumberOfElements() << endl;
	for (int i=0; i<partition->getNumberOfElements(); i++) {
		PartitionElement * element = partition->getElement(i);
		cout << setw(10) << right << i << " : ";
		cout << element->getName() << endl;
		cout << setw(10) << " " << "Best model: " << element->getBestModel()->getModel()->getName() << endl;
	}
	cout << "************************************" << endl;
//	ModelSet * models = unique_partition.getModelset();
//	for (int i=0; i<models->getNumberOfModels(); i++) {
//		models->getModel(i)->print();
//	}

	delete searchAlgo;
	delete partitionMap;
	delete options;

	return EX_OK;
}

//} /* namespace partest */
