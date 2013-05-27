/*
 * PLLAlignment.cpp
 *
 *  Created on: Jan 21, 2013
 *      Author: diego
 */

#include "PLLAlignment.h"
#include <algorithm>
#include <iostream>
#include <string.h>
#include <stdlib.h>
extern "C" {

}

namespace partest {

using namespace std;

PLLAlignment::PLLAlignment(PLLAlignment * alignment, int * firstPosition,
		int * lastPosition, int numberOfSections) {

}

PLLAlignment::PLLAlignment(PLLAlignment * alignment, int firstPosition,
		int lastPosition) {

}

PLLAlignment::PLLAlignment(string alignmentFile, DataType dataType) :
		Alignment(alignmentFile, dataType) {
	tr = (pllInstance *) malloc (sizeof (pllInstance));
	struct pllPhylip * pd;
	int *inp;
	//pd = pllPhylipParse(alignmentFile.c_str());
	//read_phylip_msa(tr, (char *) alignmentFile.c_str(), PHYLIP_SEQUENTIAL, 0);
	//numSeqs = pd->nTaxa;
	//TODO: Temporary both are equal
	//numSites = pd->seqLen;
	//numPatterns = pd->seqLen;
	//pllPhylipDump(pd);
	//pllPhylipRemoveDuplicate(pd);
	//pllPhylipDestroy(pd);
}

PLLAlignment::~PLLAlignment() {
	free (tr);
}

Alignment * PLLAlignment::splitAlignment(int firstPosition,
		int lastPosition) {

	PLLAlignment * newAlign = new PLLAlignment(this, firstPosition,
			lastPosition);

	return newAlign;
}

Alignment * PLLAlignment::splitAlignment(int * firstPosition,
		int * lastPosition, int numberOfSections) {

	PLLAlignment * newAlign = new PLLAlignment(this, firstPosition,
			lastPosition, numberOfSections);

	return newAlign;

}

} /* namespace partest */
