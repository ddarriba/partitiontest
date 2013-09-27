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

char convert(unsigned char c) {
	switch (c) {
	case 1:
		return 'A';
	case 2:
		return 'C';
	case 4:
		return 'G';
	case 8:
		return 'T';
	default:
		return 'X';
	}
}

PLLAlignment::PLLAlignment(PLLAlignment * alignment, int * firstPosition,
		int * lastPosition, int numberOfSections) {

	numSeqs = alignment->phylip->sequenceCount;
	numSites = 0;
	for (int i = 0; i < numberOfSections; i++) {
		numSites += lastPosition[i] - firstPosition[i] + 1;
	}
	phylip = pllInitAlignmentData(numSeqs, numSites);
	phylip->sequenceCount = numSeqs;
	phylip->sequenceLength = numSites;
	for (int i = 0; i < numSeqs; i++) {
		phylip->sequenceLabels[i + 1] = strdup(
				alignment->getPhylip()->sequenceLabels[i + 1]);
	}
	phylip->siteWeights = (int *) malloc(numSites * sizeof(int));
	int nextSite = 0;
	for (int i = 0; i < numberOfSections; i++) {
		for (int j = firstPosition[i] - 1; j < lastPosition[i]; j++) {
			phylip->siteWeights[nextSite] = 1;
			for (int k = 0; k < numSeqs; k++) {
				phylip->sequenceData[k + 1][nextSite] =
						alignment->getPhylip()->sequenceData[k + 1][j];
			}
			nextSite++;
		}
	}

//	for (int i = 0; i < numberOfSections; i++) {
//		cout << "[PD*] " << i + 1 << " " << firstPosition[i] << " "
//				<< lastPosition[i] << endl;
//		cout << "[PD*] ";
//	}
//	for (int k = 0; k < numSites; k++) {
//		cout << convert(phylip->sequenceData[1][k]);
//	}
//	cout << endl;

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

	struct pllPartitionRegion * pregion;
	struct pllPartitionInfo * pinfo;

	pllQueueInit(&pllPartitions);

	pinfo = (pllPartitionInfo *) malloc(sizeof(struct pllPartitionInfo));
	pllQueueInit(&(pinfo->regionList));
	pllQueueAppend(pllPartitions, (void *) pinfo);

	pinfo->partitionName = (char *) malloc(15);
	strcpy(pinfo->partitionName, "SPLITTED");
	pinfo->partitionModel = (char *) malloc(1);

	pinfo->protModels = -1;
	pinfo->protFreqs = -1;
	pinfo->dataType = DNA_DATA;
	pinfo->optimizeBaseFrequencies = PLL_TRUE;

	pregion = (struct pllPartitionRegion *) malloc(
			sizeof(struct pllPartitionRegion));
	pregion->start = 1;
	pregion->end = numSites;
	pregion->stride = 1;
	pllQueueAppend(pinfo->regionList, (void *) pregion);

	partitions = pllPartitionsCommit(pllPartitions, phylip);

	pllTreeInitTopologyForAlignment(tr, phylip);

	/* Connect the alignment with the tree structure */
	if (!pllLoadAlignment(tr, phylip, partitions, PLL_SHALLOW_COPY)) {
		cerr << "ERROR: Incompatible tree/alignment combination" << endl;
		Utilities::exit_partest(EX_SOFTWARE);
	}

	/* Initialize the model TODO: Put the parameters in a logical order and change the TRUE to flags */
	//pllInitModel(tr, PLL_TRUE, phylip, partitions);
	pllPhylipRemoveDuplicate (phylip, partitions);
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
	} else if (tr) {
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
