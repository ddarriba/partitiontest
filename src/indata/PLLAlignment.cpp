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
#include <sstream>
#include <stdlib.h>
#include <pll.h>
#include <parsePartition.h>

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
	dataType = alignment->dataType;
	numSeqs = alignment->phylip->sequenceCount;
	numSites = 0;
	uniqueName += "_";
	for (int i = 0; i < numberOfSections; i++) {
		numSites += lastPosition[i] - firstPosition[i] + 1;
		uniqueName += to_string(firstPosition[i]);
		uniqueName +=  "-";
		uniqueName += to_string(lastPosition[i]);
		uniqueName += "_";
	}
	phylip = pllInitAlignmentData(numSeqs, numSites);

	phylip->sequenceCount = numSeqs;
	phylip->sequenceLength = numSites;
	for (unsigned int i = 0; i < numSeqs; i++) {
		phylip->sequenceLabels[i + 1] = strdup(
				alignment->getPhylip()->sequenceLabels[i + 1]);
	}
	phylip->siteWeights = (int *) malloc(numSites * sizeof(int));
	int nextSite = 0;
	for (int i = 0; i < numberOfSections; i++) {
		for (int j = firstPosition[i] - 1; j < lastPosition[i]; j++) {
			phylip->siteWeights[nextSite] = 1;
			for (unsigned int k = 0; k < numSeqs; k++) {
				phylip->sequenceData[k + 1][nextSite] =
						alignment->getPhylip()->sequenceData[k + 1][j];
			}
			nextSite++;
		}
	}

	/* pllCreateInstance: int rateHetModel, int fastScaling, int saveMemory, int useRecom, long randomNumberSeed */
	pllInstanceAttr * attr = (pllInstanceAttr *) malloc(
			sizeof(pllInstanceAttr));

	attr->rateHetModel = PLL_GAMMA;
	attr->fastScaling = PLL_FALSE;
	attr->saveMemory = PLL_FALSE;
	attr->useRecom = PLL_FALSE;
	attr->randomNumberSeed = rand();
	attr->numberOfThreads = 1;

	tr = pllCreateInstance(attr);
	free(attr);

	pllPartitionRegion * pregion;
	pllPartitionInfo * pinfo;

	pllQueueInit(&pllPartitions);

	pinfo = (pllPartitionInfo *) malloc(sizeof(pllPartitionInfo));
	switch(dataType) {
		case DT_NUCLEIC:
			pinfo->protModels = -1;
			pinfo->dataType = PLL_DNA_DATA;
			break;
		case DT_PROTEIC:
			pinfo->protModels = PLL_JTT;
			pinfo->dataType = PLL_AA_DATA;
			break;
		case DT_DEFAULT:
			cerr << "ERRORS: Unknown data type";
			exit(EXIT_FAILURE);
		}
	pllQueueInit(&(pinfo->regionList));
	pllQueueAppend(pllPartitions, (void *) pinfo);

	pinfo->partitionName = (char *) malloc(15);
	strcpy(pinfo->partitionName, "SPLITTED");
	pinfo->partitionModel = (char *) malloc(1);

	pinfo->protFreqs = -1;

	pinfo->optimizeBaseFrequencies = PLL_FALSE;

	pregion = (pllPartitionRegion *) malloc(
			sizeof(pllPartitionRegion));
	pregion->start = 1;
	pregion->end = numSites;
	pregion->stride = 1;
	pllQueueAppend(pinfo->regionList, (void *) pregion);

	partitions = pllPartitionsCommit(pllPartitions, phylip);
	pllTreeInitTopologyForAlignment(tr, phylip);
	//pllPhylipRemoveDuplicate (phylip, partitions);
	numPatterns = phylip->sequenceLength;
}

PLLAlignment::PLLAlignment(PLLAlignment * alignment, int firstPosition,
		int lastPosition) {
	dataType = alignment->dataType;	phylip = 0;
	tr = 0;
	pllPartitions = 0;
	partitions = 0;
	Utilities::exit_partest(EX_SOFTWARE);
}

PLLAlignment::PLLAlignment(string alignmentFile, DataType dataType,
		pllQueue * pllPartitions) :
		Alignment(alignmentFile, dataType) {

	/* pllCreateInstance: int rateHetModel, int fastScaling, int saveMemory, int useRecom, long randomNumberSeed */
	pllInstanceAttr * attr = (pllInstanceAttr *) malloc(
			sizeof(pllInstanceAttr));
	attr->rateHetModel = PLL_GAMMA;
	attr->fastScaling = PLL_FALSE;
	attr->saveMemory = PLL_FALSE;
	attr->useRecom = PLL_FALSE;
	attr->randomNumberSeed = rand();
	attr->numberOfThreads = 1;
#ifdef DEBUG
	cout << "[TRACE] Creating PLL tree instance" << endl;
#endif
	tr = pllCreateInstance(attr);
	free(attr);
#ifdef DEBUG
	cout << "[TRACE] Creating phylip structure" << endl;
#endif
	phylip = pllParseAlignmentFile(PLL_FORMAT_PHYLIP, alignmentFile.c_str());
	/* commit the partitions and build a partitions structure */
	this->pllPartitions = pllPartitions;
#ifdef DEBUG
	cout << "[TRACE] Committing partitions " << phylip->sequenceLength << endl;
#endif
	uniqueName += "_";
	partitions = pllPartitionsCommit(pllPartitions, phylip);

	for (int i=0; i<partitions->numberOfPartitions;i++) {
		partitions->partitionData[i]->dataType = dataType;
		uniqueName += to_string(partitions->partitionData[i]->lower);
		uniqueName +=  "-";
		uniqueName += to_string(partitions->partitionData[i]->upper);
		uniqueName += "_";
	}

#ifdef DEBUG
	cout << "[TRACE] Init tree" << endl;
#endif
	pllTreeInitTopologyForAlignment(tr, phylip);

	numSeqs = phylip->sequenceCount;
	numSites = phylip->sequenceLength;
	numPatterns = phylip->sequenceLength;
#ifdef DEBUG
	cout << "[TRACE] Done " << numSeqs << " x " << numSites << endl;
#endif
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
