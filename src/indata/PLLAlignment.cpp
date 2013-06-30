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

void PLLAlignment::buildTopology(void) {

}

PLLAlignment::PLLAlignment(PLLAlignment * alignment, int * firstPosition,
		int * lastPosition, int numberOfSections) {

	/* pllCreateInstance: int rateHetModel, int fastScaling, int saveMemory, int useRecom, long randomNumberSeed */
	tr = pllCreateInstance(GAMMA, PLL_FALSE, PLL_FALSE, PLL_FALSE, 12345);

	numSeqs = alignment->phylip->nTaxa;
	numSites = 0;
	for (int cur_part = 0; cur_part < numberOfSections; cur_part++) {
		numSites += lastPosition[cur_part] - firstPosition[cur_part] + 1;

	}

	phylip = (struct pllPhylip *) malloc(sizeof(struct pllPhylip));
	phylip->nTaxa = numSeqs;
	phylip->label = alignment->phylip->label;
	phylip->seqLen = numSites;
	phylip->seq = (unsigned char **) malloc(
			(numSeqs + 1) * sizeof(unsigned char *));
	for (int i = 1; i <= numSeqs; i++) {
		phylip->seq[i] = (unsigned char *) malloc(
				numSites * sizeof(unsigned char));
	}

	int cur_position = 0;

	for (int cur_part = 0; cur_part < numberOfSections; cur_part++) {
		int partNumSites = lastPosition[cur_part] - firstPosition[cur_part] + 1;
		for (int i = 1; i <= numSeqs; i++) {
			for (int site = 0; site < partNumSites; site++) {
				phylip->seq[i][cur_position + site] =
						alignment->phylip->seq[i][firstPosition[cur_part] + site
								- 1];
			}
		}

		phylip->weights = (int *) rax_malloc((phylip->seqLen) * sizeof(int));
		for (int site = 0; site < partNumSites; site++) {
			phylip->weights[cur_position + site] = 1;
		}

		cur_position += partNumSites;
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

//	for (int cur_part = 0; cur_part < numberOfSections; cur_part++) {
	for (int i = 0; i < alignment->partitions->numberOfPartitions; i++) {
		if (alignment->partitions->partitionData[i]->lower
				== (firstPosition[0] - 1)) {
			pi->partitionName =
					alignment->partitions->partitionData[i]->partitionName;
			break;
		}
	}
//	}
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

	//read_phylip_msa(tr, (char *) alignmentFile.c_str(), PHYLIP_SEQUENTIAL, 0);

	numPatterns = phylip->seqLen;
	pllTreeInitTopologyForAlignment(tr, phylip);

	/* Connect the alignment with the tree structure */
	if (!pllLoadAlignment(tr, phylip, partitions, PLL_DEEP_COPY)) {
		cerr << "ERROR: Incompatible tree/alignment combination" << endl;
		Utilities::exit_partest(EX_SOFTWARE);
	}

	//pllTreeInitTopologyRandom(tr, phylip->nTaxa, phylip->label);

	/* Initialize the model TODO: Put the parameters in a logical order and change the TRUE to flags */
	pllInitModel(tr, PLL_TRUE, phylip, partitions);

}

PLLAlignment::PLLAlignment(PLLAlignment * alignment, int firstPosition,
		int lastPosition) {
	/* pllCreateInstance: int rateHetModel, int fastScaling, int saveMemory, int useRecom, long randomNumberSeed */
	tr = pllCreateInstance(GAMMA, PLL_FALSE, PLL_FALSE, PLL_FALSE, 12345);

	numSeqs = alignment->phylip->nTaxa;
	numSites = lastPosition - firstPosition + 1;
	/* create new phylip structure */
	phylip = (struct pllPhylip *) malloc(sizeof(struct pllPhylip));
	phylip->nTaxa = numSeqs;
	phylip->label = alignment->phylip->label;
	phylip->seqLen = numSites;
	phylip->seq = (unsigned char **) malloc(
			(numSeqs + 1) * sizeof(unsigned char *));
	for (int i = 1; i <= numSeqs; i++) {
		phylip->seq[i] = (unsigned char *) malloc(
				numSites * sizeof(unsigned char));
		for (int site = 0; site < numSites; site++) {
			phylip->seq[i][site] = alignment->phylip->seq[i][firstPosition
					+ site - 1];
		}
	}

	phylip->weights = (int *) rax_malloc((phylip->seqLen) * sizeof(int));
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
				== (firstPosition - 1)) {
			pi->partitionName =
					alignment->partitions->partitionData[i]->partitionName;
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

	//read_phylip_msa(tr, (char *) alignmentFile.c_str(), PHYLIP_SEQUENTIAL, 0);

	numPatterns = phylip->seqLen;
	pllTreeInitTopologyForAlignment(tr, phylip);

	/* Connect the alignment with the tree structure */
	if (!pllLoadAlignment(tr, phylip, partitions, PLL_DEEP_COPY)) {
		cerr << "ERROR: Incompatible tree/alignment combination" << endl;
		Utilities::exit_partest(EX_SOFTWARE);
	}

	pllTreeInitTopologyRandom(tr, phylip->nTaxa, phylip->label);

	/* Initialize the model TODO: Put the parameters in a logical order and change the TRUE to flags */
	pllInitModel(tr, PLL_TRUE, phylip, partitions);

}

PLLAlignment::PLLAlignment(string alignmentFile, DataType dataType,
		string partitionsFile) :
		Alignment(alignmentFile, dataType) {

	/* pllCreateInstance: int rateHetModel, int fastScaling, int saveMemory, int useRecom, long randomNumberSeed */
	tr = pllCreateInstance(GAMMA, PLL_FALSE, PLL_FALSE, PLL_FALSE, 12345);
	phylip = pllPhylipParse(alignmentFile.c_str());

	parts = pllPartitionParse(partitionsFile.c_str());

	/* commit the partitions and build a partitions structure */
	partitions = pllPartitionsCommit(parts, phylip);
	//pllQueuePartitionsDestroy(&parts);
	cout << "NTAXA = " << phylip->nTaxa << endl;
	cout << "SEQLEN = " << phylip->seqLen << endl;

	numSeqs = phylip->nTaxa;
	numSites = phylip->seqLen;

	//pllPhylipRemoveDuplicate(phylip, partitions);

	numPatterns = phylip->seqLen;

	pllTreeInitTopologyForAlignment(tr, phylip);

	/* Connect the alignment with the tree structure */
	if (!pllLoadAlignment(tr, phylip, partitions, PLL_SHALLOW_COPY)) {
		cerr << "ERROR: Incompatible tree/alignment combination" << endl;
		Utilities::exit_partest(EX_SOFTWARE);
	}

	/* Initialize the model TODO: Put the parameters in a logical order and change the TRUE to flags */
	pllInitModel(tr, PLL_TRUE, phylip, partitions);
}

PLLAlignment::~PLLAlignment() {
	free(tr);
}

Alignment * PLLAlignment::splitAlignment(int firstPosition, int lastPosition) {

	cout << "SPLITTING " << firstPosition << " to " << lastPosition << endl;

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
