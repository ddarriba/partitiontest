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
}

namespace partest {

using namespace std;

PLLAlignment::PLLAlignment(PLLAlignment * alignment, int * firstPosition,
		int * lastPosition, int numberOfSections) {
	/* pllCreateInstance: int rateHetModel, int fastScaling, int saveMemory, int useRecom, long randomNumberSeed */
	tr = pllCreateInstance(GAMMA, PLL_FALSE, PLL_FALSE, PLL_FALSE, rand());

	numSeqs = alignment->phylip->nTaxa;
	numSites = 0;
	for (int cur_part = 0; cur_part < numberOfSections; cur_part++) {
		numSites += lastPosition[cur_part] - firstPosition[cur_part] + 1;
	}
	phylip = (struct pllPhylip *) malloc(sizeof(struct pllPhylip));
	phylip->nTaxa = numSeqs;
	phylip->label = (char **) calloc((numSeqs + 1), sizeof(char *));
	for (int i = 0; i < numSeqs; i++) {
		phylip->label[i + 1] = strndup(alignment->phylip->label[i + 1],
				strlen(alignment->phylip->label[i + 1]));
	}
	phylip->seqLen = numSites;
	phylip->seq = (unsigned char **) malloc(
			(numSeqs + 1) * sizeof(unsigned char *));

	phylip->seq[1] = (unsigned char *) malloc(
			(numSites + 1) * numSeqs * sizeof(unsigned char));

	phylip->weights = (int *) rax_malloc(numSites * sizeof(int));
	for (int i = 1; i <= numSeqs; i++) {
		phylip->seq[i] = phylip->seq[1] + (i - 1) * (numSites + 1);
		phylip->seq[i][numSites] = 0;
		int cur_position = 0;

		for (int cur_part = 0; cur_part < numberOfSections; cur_part++) {
			int partNumSites = lastPosition[cur_part] - firstPosition[cur_part]
					+ 1;
			for (int site = 0; site < partNumSites; site++) {
				phylip->seq[i][cur_position + site] =
						alignment->phylip->seq[i][firstPosition[cur_part] + site
								- 1];
			}
			cur_position += partNumSites;
		}
	}
	for (int site = 0; site < numSites; site++) {
		phylip->weights[site] = 1;
	}
	struct pllPartitionInfo * pi;
	struct pllPartitionRegion * region;
	pi = (struct pllPartitionInfo *) rax_calloc(1,
			sizeof(struct pllPartitionInfo));
	pllQueueInit(&parts);
	pllQueueInit(&(pi->regionList));
	pllQueueAppend(parts, (void *) pi);
	pi->partitionModel = (char *) malloc(5 * sizeof(char));
	strcpy(pi->partitionModel, "DNA");

	for (int i = 0; i < alignment->partitions->numberOfPartitions; i++) {
		if (alignment->partitions->partitionData[i]->lower
				== (firstPosition[0] - 1)) {
			pi->partitionName = strdup(
					alignment->partitions->partitionData[i]->partitionName);
			break;
		}
	}

	pi->optimizeBaseFrequencies = PLL_TRUE;
	pi->dataType = DNA_DATA;
	pi->protFreqs = PLL_FALSE;
	region = (struct pllPartitionRegion *) rax_malloc(
			sizeof(struct pllPartitionRegion));
	region->start = 1;
	region->stride = 1;
	region->end = numSites;
	pllQueueAppend(pi->regionList, (void *) region);

	/* commit the partitions and build a partitions structure */
	partitions = pllPartitionsCommit(parts, phylip);

	pllPartitionsValidate(parts, phylip);

	pllPhylipRemoveDuplicate(phylip, partitions);
	/* destroy the  intermedia partition queue structure */
	pllQueuePartitionsDestroy(&parts);

	numPatterns = phylip->seqLen;

}

PLLAlignment::PLLAlignment(PLLAlignment * alignment, int firstPosition,
		int lastPosition) {
	phylip = 0;
	tr = 0;
	parts = 0;
	partitions = 0;
	Utilities::exit_partest(EX_SOFTWARE);
}

PLLAlignment::PLLAlignment(string alignmentFile, DataType dataType,
		string partitionsFile) :
		Alignment(alignmentFile, dataType) {

	/* pllCreateInstance: int rateHetModel, int fastScaling, int saveMemory, int useRecom, long randomNumberSeed */
	tr = pllCreateInstance(GAMMA, PLL_FALSE, PLL_FALSE, PLL_FALSE, rand());
	phylip = pllPhylipParse(alignmentFile.c_str());

	/* commit the partitions and build a partitions structure */
	parts = pllPartitionParse(partitionsFile.c_str());
	partitions = pllPartitionsCommit(parts, phylip);
	pllQueuePartitionsDestroy(&parts);

	cout << "NTAXA = " << phylip->nTaxa << endl;
	cout << "SEQLEN = " << phylip->seqLen << endl;

	pllTreeInitTopologyForAlignment(tr, phylip);

	/* Connect the alignment with the tree structure */
	if (!pllLoadAlignment(tr, phylip, partitions, PLL_SHALLOW_COPY)) {
		cerr << "ERROR: Incompatible tree/alignment combination" << endl;
		Utilities::exit_partest(EX_SOFTWARE);
	}

	/* Initialize the model TODO: Put the parameters in a logical order and change the TRUE to flags */
	pllInitModel(tr, PLL_TRUE, phylip, partitions);

	numSeqs = phylip->nTaxa;
	numSites = phylip->seqLen;

//pllPhylipRemoveDuplicate(phylip, partitions);

	numPatterns = phylip->seqLen;

}

PLLAlignment::~PLLAlignment() {
	if (phylip)
		pllPhylipDestroy(phylip);
	if (tr && tr->mxtips > 0) {
		if (partitions) {
			pllPartitionsDestroy(&partitions, partitions->numberOfPartitions,
					tr->mxtips);
		}
		pllTreeDestroy(tr);
	} else if (tr) {
		pllTreeDestroy(tr);
	}
}

void PLLAlignment::destroyStructures(void) {
	if (phylip)
		pllPhylipDestroy(phylip);
	if (partitions)
		pllPartitionsDestroy(&partitions, partitions->numberOfPartitions,
				tr->mxtips);
	if (tr)
		pllTreeDestroy(tr);
	phylip = 0;
	partitions = 0;
	tr = 0;
}

Alignment * PLLAlignment::splitAlignment(int firstPosition, int lastPosition) {
	PLLAlignment * newAlign = new PLLAlignment(this, &firstPosition,
			&lastPosition, 1);
	return newAlign;
}

Alignment * PLLAlignment::splitAlignment(int * firstPosition,
		int * lastPosition, int numberOfSections) {
	PLLAlignment * newAlign = new PLLAlignment(this, firstPosition,
			lastPosition, numberOfSections);
	return newAlign;
}

} /* namespace partest */
