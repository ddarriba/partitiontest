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
#include <time.h>
#include "util/PrintMeta.h"
#include "util/Utilities.h"
#include "util/GlobalDefs.h"
#include "util/ParTestFactory.h"
#include "options/ParTestOptions.h"
#include "parser/ArgumentParser.h"
#include "indata/PartitionMap.h"
#include "selection/ModelSelector.h"
#include "selection/PartitionSelector.h"
#include "observer/ConsoleObserver.h"

using namespace std;
using namespace partest;

int main(int argc, char *argv[]) {

	time_t iniTime, endTime;

	time(&iniTime);

	PrintMeta::print_header(cout);

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

	PartitionMap * partitionMap = new PartitionMap(options->getConfigFile(),
			options->getAlignment(), options->getRateVariation(),
			options->getDataType());

	SearchAlgorithm * searchAlgo = ParTestFactory::createSearchAlgorithm(
			options, partitionMap);

#ifdef DEBUG
	cout << "[TRACE] Starting search of the best partitioning scheme" << endl;
#endif

	PartitioningScheme * partitioningScheme = searchAlgo->start();

#ifdef DEBUG
	cout << "[TRACE] End of search of the best partitioning scheme" << endl;
#endif

	cout << "Search done... it took " << time(NULL) - iniTime << " seconds."
			<< endl;
	cout << "Best Scheme:" << endl << partitioningScheme->getCode() << endl
			<< partitioningScheme->getName() << endl;

#ifdef FAST_DNA

#ifdef DEBUG
	cout << "[TRACE] Optimizing best scheme: " << partitioningScheme->getCode() << endl;
#endif

	switch(options->getOptimizeMode()) {
		case OPT_SEARCH:
		partitioningScheme->buildCompleteModelSet(false, iotui);
		ModelOptimize * mo = ParTestFactory::createModelOptimize(options);
		ConsoleObserver * observer = new ConsoleObserver();
		mo->attach(observer);
		mo->optimizePartitioningScheme(partitioningScheme, false, 1, 1);
		PartitionSelector partSelector(&partitioningScheme, 1, options);
		break;
		case OPT_GTR:
		break;
	}
#ifdef DEBUG
	cout << "[TRACE] Optimizing best scheme at once" << endl;
#endif

	static_cast<PLLModelOptimize* >(mo)->optimizePartitioningSchemeAtOnce(partitioningScheme);

	delete mo;
	delete observer;
#endif

#ifdef DEBUG
	cout << "[TRACE] Done: " << partitioningScheme->getCode() << endl;
#endif

	cout << "Optimization done... it took " << time(NULL) - iniTime
			<< " seconds." << endl;

	if (!partitioningScheme->isOptimized()) {
		cerr << endl
				<< "[ERROR] Logic error. Search algorithm returned a non-optimized partition"
				<< endl;
		Utilities::exit_partest(EX_SOFTWARE);
	}

	cout << "***** RESULTS *****" << endl;
	PrintMeta::print_options(cout, *options);
	cout << "Number of elements:  " << partitioningScheme->getNumberOfElements()
			<< endl;
	cout << "Partitioning scheme: " << endl << partitioningScheme->getCode()
			<< endl;
	for (int i = 0; i < partitioningScheme->getNumberOfElements(); i++) {
		PartitionElement * element = partitioningScheme->getElement(i);
		cout << setw(10) << right << i << " : ";
		cout << element->getName() << endl;
		cout << setw(10) << " " << "Best model: "
				<< element->getBestModel()->getModel()->getName() << endl;
	}
	cout << "************************************" << endl;

	ofstream * rout = options->getResultsOutputStream();
	PrintMeta::print_options(*rout, *options);
	*rout << "Number of elements:  "
			<< partitioningScheme->getNumberOfElements() << endl;
	*rout << "Partitioning scheme: " << endl << partitioningScheme->getCode()
			<< endl;
	for (int i = 0; i < partitioningScheme->getNumberOfElements(); i++) {
		PartitionElement * element = partitioningScheme->getElement(i);
		*rout << setw(10) << right << i << " : ";
		*rout << element->getName() << endl;
		*rout << setw(10) << " " << "Best model: "
				<< element->getBestModel()->getModel()->getName() << endl;
	}

	*rout << endl;

	for (int i = 0; i < partitioningScheme->getNumberOfElements(); i++) {
		PartitionElement * element = partitioningScheme->getElement(i);
		*rout << setw(10) << right << i << " : ";
		*rout << element->getName() << endl;
		element->getBestModel()->getModel()->print(*rout);
	}

	if (partitioningScheme->getNumberOfElements() == 1) {
		cout << "BEST TREE: "
				<< partitioningScheme->getElement(0)->getBestModel()->getModel()->getTree()
				<< endl;
	}

	delete searchAlgo;
	delete partitionMap;
	delete options;

	return EX_OK;
}
