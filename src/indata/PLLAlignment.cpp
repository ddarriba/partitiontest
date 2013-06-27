/*
 * PLLAlignment.cpp
 *
 *  Created on: Jan 21, 2013
 *      Author: diego
 */

#include "PLLAlignment.h"
#include "util/Utilities.h"
#include <algorithm>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#ifndef AXML_H
#define AXML_H
#include "axml.h"
#endif
extern "C" {
#include "parser/phylip/phylip.h"
#include "utils.h"
#include "parser/partition/part.h"
}

namespace partest {

using namespace std;

PLLAlignment::PLLAlignment(PLLAlignment * alignment, int * firstPosition,
		int * lastPosition, int numberOfSections) {
}

PLLAlignment::PLLAlignment(PLLAlignment * alignment, int firstPosition,
		int lastPosition) {
}

PLLAlignment::PLLAlignment(string alignmentFile, DataType dataType,
		string partitionsFile) :
		Alignment(alignmentFile, dataType) {

	/* pllCreateInstance: int rateHetModel, int fastScaling, int saveMemory, int useRecom, long randomNumberSeed */
	tr = pllCreateInstance(GAMMA, PLL_FALSE, PLL_FALSE, PLL_FALSE, 12345);
	phylip = pllPhylipParse(alignmentFile.c_str());

	struct pllQueue * parts = pllPartitionParse(partitionsFile.c_str());
	/* commit the partitions and build a partitions structure */
	partitions = pllPartitionsCommit(parts, phylip);
	/* destroy the  intermedia partition queue structure */
	pllQueuePartitionsDestroy(&parts);

	cout << "NTAXA = " << phylip->nTaxa << endl;
	cout << "SEQLEN = " << phylip->seqLen << endl;

	//read_phylip_msa(tr, (char *) alignmentFile.c_str(), PHYLIP_SEQUENTIAL, 0);
	numSeqs = phylip->nTaxa;
	//TODO: Temporary both are equal
	numSites = phylip->seqLen;
	numPatterns = phylip->seqLen;
}

PLLAlignment::~PLLAlignment() {
	free(tr);
}

Alignment * PLLAlignment::splitAlignment(int firstPosition, int lastPosition) {

	PLLAlignment * newAlign = new PLLAlignment(this, firstPosition,
			lastPosition);

	return newAlign;
}

Alignment * PLLAlignment::splitAlignment(int * firstPosition,
		int * lastPosition, int numberOfSections) {

	return this;
//	PLLAlignment * newAlign = new PLLAlignment(this, firstPosition,
//			lastPosition, numberOfSections);
//
//	return newAlign;

}

} /* namespace partest */
