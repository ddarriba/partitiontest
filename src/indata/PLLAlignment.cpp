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
#include "pll.h"
#endif

namespace partest {

using namespace std;

PLLAlignment::PLLAlignment(PLLAlignment * alignment, int * firstPosition,
		int * lastPosition, int numberOfSections) {

	/* pllCreateInstance: int rateHetModel, int fastScaling, int saveMemory, int useRecom, long randomNumberSeed */
	pllInstanceAttr * attr = (pllInstanceAttr *) rax_malloc(
			sizeof(pllInstanceAttr));
	attr->rateHetModel = GAMMA;
	attr->fastScaling = PLL_FALSE;
	attr->saveMemory = PLL_FALSE;
	attr->useRecom = PLL_FALSE;
	attr->randomNumberSeed = rand();
	attr->numberOfThreads = 1;

	tr = pllCreateInstance(attr);
	rax_free(attr);
	numSeqs = alignment->phylip->sequenceCount;
	numSites = 0;
	for (int cur_part = 0; cur_part < numberOfSections; cur_part++) {
		numSites += lastPosition[cur_part] - firstPosition[cur_part] + 1;
	}

	phylip = (pllAlignmentData *) malloc(sizeof(pllAlignmentData));
	phylip->sequenceCount = numSeqs;
	phylip->sequenceLabels = (char **) calloc((numSeqs + 1), sizeof(char *));
	for (int i = 0; i < numSeqs; i++) {
		phylip->sequenceLabels[i + 1] = strndup(
				alignment->phylip->sequenceLabels[i + 1],
				strlen(alignment->phylip->sequenceLabels[i + 1]));
	}
	phylip->sequenceLength = numSites;
	phylip->sequenceData = (unsigned char **) malloc(
			(numSeqs + 1) * sizeof(unsigned char *));

	phylip->sequenceData[1] = (unsigned char *) malloc(
			(numSites + 1) * numSeqs * sizeof(unsigned char));

	phylip->siteWeights = (int *) rax_malloc(numSites * sizeof(int));
	for (int i = 1; i <= numSeqs; i++) {
		phylip->sequenceData[i] = phylip->sequenceData[1]
				+ (i - 1) * (numSites + 1);
		phylip->sequenceData[i][numSites] = 0;
		int cur_position = 0;

		for (int cur_part = 0; cur_part < numberOfSections; cur_part++) {
			int partNumSites = lastPosition[cur_part] - firstPosition[cur_part]
					+ 1;
			for (int site = 0; site < partNumSites; site++) {
				phylip->sequenceData[i][cur_position + site] =
						alignment->phylip->sequenceData[i][firstPosition[cur_part]
								+ site - 1];
			}
			cur_position += partNumSites;
		}
	}

	for (int site = 0; site < numSites; site++) {
		phylip->siteWeights[site] = 1;
	}
	struct pllPartitionInfo * pi;
	struct pllPartitionRegion * region;
	pi = (struct pllPartitionInfo *) rax_calloc(1,
			sizeof(struct pllPartitionInfo));
	pllQueueInit(&pllPartitions);
	pllQueueInit(&(pi->regionList));
	pllQueueAppend(pllPartitions, (void *) pi);
	pi->partitionModel = (char *) malloc(1);

	for (int i = 0; i < alignment->partitions->numberOfPartitions; i++) {
		if (alignment->partitions->partitionData[i]->lower
				== (firstPosition[0] - 1)) {
			pi->partitionName = strdup(
					alignment->partitions->partitionData[i]->partitionName);
			break;
		}
	}

	pi->protFreqs = -1;
	pi->protModels = -1;
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
	partitions = pllPartitionsCommit(pllPartitions, phylip);
	pllPartitionsValidate(pllPartitions, phylip);
	pllPhylipRemoveDuplicate(phylip, partitions);
	/* destroy the  intermedia partition queue structure */
	pllQueuePartitionsDestroy(&pllPartitions);
	numPatterns = phylip->sequenceLength;

}

PLLAlignment::PLLAlignment(PLLAlignment * alignment, int firstPosition,
		int lastPosition) {
	phylip = 0;
	tr = 0;
	pllPartitions = 0;
	partitions = 0;
	Utilities::exit_partest(EX_SOFTWARE);
}

PLLAlignment::PLLAlignment(string alignmentFile, DataType dataType,
		struct pllQueue * pllPartitions) :
		Alignment(alignmentFile, dataType) {

	/* pllCreateInstance: int rateHetModel, int fastScaling, int saveMemory, int useRecom, long randomNumberSeed */
	pllInstanceAttr * attr = (pllInstanceAttr *) rax_malloc(
			sizeof(pllInstanceAttr));
	attr->rateHetModel = GAMMA;
	attr->fastScaling = PLL_FALSE;
	attr->saveMemory = PLL_FALSE;
	attr->useRecom = PLL_FALSE;
	attr->randomNumberSeed = rand();
	attr->numberOfThreads = 1;
	tr = pllCreateInstance(attr);
	rax_free(attr);

	phylip = pllParsePHYLIP(alignmentFile.c_str());

	/* commit the partitions and build a partitions structure */
	this->pllPartitions = pllPartitions;

	partitions = pllPartitionsCommit(pllPartitions, phylip);

	pllTreeInitTopologyForAlignment(tr, phylip);

	/* Connect the alignment with the tree structure */
	if (!pllLoadAlignment(tr, phylip, partitions, PLL_SHALLOW_COPY)) {
		cerr << "ERROR: Incompatible tree/alignment combination" << endl;
		Utilities::exit_partest(EX_SOFTWARE);
	}

	/* Initialize the model TODO: Put the parameters in a logical order and change the TRUE to flags */
	//pllInitModel(tr, PLL_TRUE, phylip, partitions);

	numSeqs = phylip->sequenceCount;
	numSites = phylip->sequenceLength;
	numPatterns = phylip->sequenceLength;
}

PLLAlignment::~PLLAlignment() {
	if (phylip)
		pllAlignmentDataDestroy(phylip);
	if (tr && tr->mxtips > 0) {
		if (partitions) {
			pllPartitionsDestroy(tr, &partitions);
		}
		pllDestroyInstance(tr);
	}
	else if (tr) {
		pllDestroyInstance(tr);
	}
}

void PLLAlignment::destroyStructures(void) {
	if (phylip)
		pllAlignmentDataDestroy(phylip);
	if (partitions)
		pllPartitionsDestroy(tr, &partitions);
	if (tr)
	pllDestroyInstance(tr);
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
