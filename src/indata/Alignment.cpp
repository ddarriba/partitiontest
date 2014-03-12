/*
 * Alignment.cpp
 *
 *  Created on: Jan 7, 2013
 *      Author: diego
 */

#include "Alignment.h"
#include "util/Utilities.h"
#include "util/GlobalDefs.h"
#include <fstream>
#include <iostream>

namespace partest {

Alignment::Alignment() {
	numSeqs = 0;
	numSites = 0;
	numPatterns = 0;
	numberOfFrequencies = 0;
	empiricalFrequencies = 0;
	shannonEntropy = 0.0;
	dataType = DT_NUCLEIC;
}
Alignment::Alignment(string alignmentFile, DataType dataType) :
		alignmentFile(alignmentFile), dataType(dataType), numPatterns(0),
		shannonEntropy(0), numSites(0), numSeqs(0), empiricalFrequencies(0) {

	ifstream inputFile(alignmentFile.c_str());
	if (!inputFile.good()) {
		cerr << "[ERROR] Input file \"" << alignmentFile
				<< "\" does not exist or you don't have permission on it"
				<< endl;
		Utilities::exit_partest(EX_DATAERR);
	}

	switch (dataType) {
	case DT_NUCLEIC:
		numberOfFrequencies = 4;
		break;
	case DT_PROTEIC:
		numberOfFrequencies = 20;
		break;
	}
}

Alignment::~Alignment() {
}

Alignment * Alignment::splitAlignment(int firstPosition, int lastPosition) {
#ifdef DEBUG
	cout << "[TRACE] MOCK FUNCTION Split Alignment" << endl;
#endif
	return new Alignment(alignmentFile, dataType);
}

Alignment * Alignment::splitAlignment(int * firstPosition, int * lastPosition,
		int numberOfSections) {
#ifdef DEBUG
	cout << "[TRACE] MOCK FUNCTION Split Alignment" << endl;
#endif
	return new Alignment(alignmentFile, dataType);
}

double Alignment::computeShannonEntropy(double *observedFrequencies) {
	shannonEntropy = 0.0;
	int i;
	for (i = 0; i < numberOfFrequencies; i++) {
		shannonEntropy += observedFrequencies[i]
				* Utilities::dBinaryLog(observedFrequencies[i]);
	}
	shannonEntropy *= -1;
	return shannonEntropy;
}

} /* namespace partest */

