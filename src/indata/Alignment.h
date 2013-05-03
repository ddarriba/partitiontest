/*  PartitionTest, fast selection of the best fit partitioning scheme for
 *  multi-gene data sets.
 *  Copyright May 2013 by Diego Darriba
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 3 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 *
 *  For any other enquiries send an Email to Diego Darriba
 *  ddarriba@udc.es
 */

/**
 * @file Alignment.h
 *
 * @brief Generic implementation for alignment definition
 */
#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_

#include <string>
#include "util/GlobalDefs.h"

using namespace std;

namespace partest {

class Alignment {
protected:
	double computeShannonEntropy(double *observedFrequencies);
	string alignmentFile; /** alignment input file */
	unsigned int numSites; /** number of sites in each sequence */
	unsigned int numPatterns; /** number of different site patterns */
	unsigned int numSeqs; /** number of sequences */
	double *empiricalFrequencies; /** observed frequencies of nucleotides/aminoacids */
	int numberOfFrequencies; /** 4 for nucleotides, 20 for amino-acids */
	DataType dataType; /** nucleotide of aminoacid data */
	double shannonEntropy; /** Shannon Entropy of the alignment */
public:
	Alignment() {
	}
	Alignment(string alignmentFile, DataType dataType);
	virtual ~Alignment();
	virtual Alignment * splitAlignment(int firstPosition, int lastPosition);
	virtual Alignment * splitAlignment(int * firstPosition, int * lastPosition,
			int numberOfSections);
	string getAlignmentFile() const {
		return alignmentFile;
	}
	unsigned int getNumSites() const {
		return numSites;
	}
	unsigned int getNumSeqs() const {
		return numSeqs;
	}
	unsigned int getNumPatterns() const {
		return numPatterns;
	}
	double getShannonEntropy() const {
		return shannonEntropy;
	}
};

} /* namespace partest */

#endif /* ALIGNMENT_H_ */
