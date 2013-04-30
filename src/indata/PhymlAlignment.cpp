/*
 * PhymlAlignment.cpp
 *
 *  Created on: Jan 21, 2013
 *      Author: diego
 */

#include "PhymlAlignment.h"
#include <algorithm>
#include <iostream>
#include <string.h>
#include <stdlib.h>
extern "C" {
#include "partest_to_phyml.h"
}

namespace partest {

using namespace std;

struct __Calign * PhymlAlignment::build_cdata(PhymlAlignment * alignment,
		int * firstPosition, int * lastPosition, int numberOfSections) {
	struct __Calign * cdata = (struct __Calign *) malloc(
			sizeof(struct __Calign));
	// TODO: Recompute Base Frequencies
	int i, cur_index;

	cdata->format = alignment->cdata->format;

	cdata->init_len = 0;
	for (i = 0; i < numberOfSections; i++) {
		cdata->init_len += lastPosition[i] - firstPosition[i] + 1;
	}

	int numSeqs = alignment->cdata->n_otu;
	cdata->n_otu = alignment->cdata->n_otu;

	int numSites = cdata->init_len;
	cdata->sitepatt = (int *) malloc(numSites * sizeof(int));

	cur_index = 0;
	for (i = 0; i < numberOfSections; i++) {
#ifdef DEBUG
		cout << "[TRACE] PhyML Alignment. Adding section from " << firstPosition[i] << " to " << lastPosition[i] << endl;
#endif
		int sectionSites = lastPosition[i] - firstPosition[i] + 1;
		memcpy(&(cdata->sitepatt[cur_index]),
				&(alignment->cdata->sitepatt[firstPosition[i] - 1]),
				sectionSites * sizeof(int));
		cur_index += sectionSites;
	}

	std::sort(cdata->sitepatt, &(cdata->sitepatt[numSites]));

	/* count number of patterns */
	int lastValue = -1;

	int numPatterns = 0;
	for (int i = 0; i < numSites; i++) {
		if (cdata->sitepatt[i] != lastValue) {
			numPatterns++;
			lastValue = cdata->sitepatt[i];
		}
	}
	cdata->crunch_len = numPatterns;
	cdata->wght = (int *) malloc(numPatterns * sizeof(int));
	cdata->invar = (short int *) malloc(numPatterns * sizeof(short int));
	cdata->ambigu = (short int *) malloc(numPatterns * sizeof(short int));
	/* build indirection array for renaming and weights */
	int reindex[numPatterns];
	int currentIndex = 0;
	int currentWeight = 0;
	lastValue = -1;
	int invarCount = 0;
	for (int i = 0; i < numSites; i++) {
		currentWeight++;
		if (cdata->sitepatt[i] != lastValue) {
			reindex[currentIndex] = cdata->sitepatt[i];
			cdata->invar[currentIndex] =
					alignment->cdata->invar[cdata->sitepatt[i]];
			if (currentIndex > 0) {
				cdata->wght[currentIndex - 1] = currentWeight;
			}
			currentIndex++;
			if (i < (numSites - 1))
				currentWeight = 0;
			lastValue = cdata->sitepatt[i];
		}
	}

	cdata->wght[currentIndex - 1] = currentWeight;
	for (int i = 0; i < numPatterns; i++) {
		if (cdata->invar[i] > -1) {
			invarCount += cdata->wght[i];
		}
	}
	cdata->obs_pinvar = (double) invarCount / numSites;
	/* build sequences */
	cdata->c_seq = (struct __Align **) malloc(
			numSeqs * sizeof(struct __Align *));

	for (int i = 0; i < numSeqs; i++) {

		cdata->c_seq[i] = (struct __Align *) malloc(sizeof(struct __Align));
		cdata->c_seq[i]->state = (char *) malloc(numPatterns);
		cdata->c_seq[i]->is_ambigu = (short int *) malloc(
				numPatterns * sizeof(short int));
		cdata->c_seq[i]->name = (char *) malloc(
				strlen(alignment->cdata->c_seq[i]->name) + 1);
		strcpy(cdata->c_seq[i]->name, alignment->cdata->c_seq[i]->name);
		cdata->c_seq[i]->len = numPatterns;
		for (int j = 0; j < numPatterns; j++) {
			cdata->c_seq[i]->state[j] =
					alignment->cdata->c_seq[i]->state[reindex[j]]; //cdata->sitepatt[j]];
			cdata->c_seq[i]->is_ambigu[j] =
					alignment->cdata->c_seq[i]->is_ambigu[reindex[j]]; //cdata->sitepatt[j]];
		}
	}

	for (int j = 0; j < numPatterns; j++) {
		cdata->ambigu[j] = alignment->cdata->ambigu[cdata->sitepatt[j]];
	}

	if (dataType == DT_PROTEIC) {
		cdata->b_frq = (phydbl *) malloc(20 * sizeof(phydbl));
		get_aa_freqs(cdata);
	} else {
		cdata->b_frq = (phydbl *) malloc(4 * sizeof(phydbl));
		get_nt_freqs(cdata);
	}

	return cdata;
}

PhymlAlignment::PhymlAlignment(PhymlAlignment * alignment, int * firstPosition,
		int * lastPosition, int numberOfSections) {

	dataType = alignment->dataType;
	cdata = build_cdata(alignment, firstPosition, lastPosition,
			numberOfSections);
	numSeqs = cdata->n_otu;
	numSites = cdata->init_len;
	numPatterns = cdata->crunch_len;
}

PhymlAlignment::PhymlAlignment(PhymlAlignment * alignment, int firstPosition,
		int lastPosition) {

	dataType = alignment->dataType;
	cdata = build_cdata(alignment, &firstPosition, &lastPosition, 1);
	numSeqs = cdata->n_otu;
	numSites = cdata->init_len;
	numPatterns = cdata->crunch_len;
}

PhymlAlignment::PhymlAlignment(string alignmentFile, DataType dataType) :
		Alignment(alignmentFile, dataType) {

	cdata = read_data(alignmentFile.c_str(), dataType);

	numSeqs = cdata->n_otu;
	numSites = cdata->init_len;
	numPatterns = cdata->crunch_len;

	computeShannonEntropy(cdata->b_frq);
	// TEST!!!
//	PhymlAlignment * align = splitAlignment(1, numSites);
//	cdata = align->cdata;

}

PhymlAlignment::~PhymlAlignment() {
	free_cdata(cdata);
}

struct __Calign *PhymlAlignment::getCData() {
	return cdata;
}

Alignment * PhymlAlignment::splitAlignment(int firstPosition,
		int lastPosition) {

	PhymlAlignment * newAlign = new PhymlAlignment(this, firstPosition,
			lastPosition);

	return newAlign;
}

Alignment * PhymlAlignment::splitAlignment(int * firstPosition,
		int * lastPosition, int numberOfSections) {

	PhymlAlignment * newAlign = new PhymlAlignment(this, firstPosition,
			lastPosition, numberOfSections);

	return newAlign;

}

} /* namespace partest */
