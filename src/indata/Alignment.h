/*
 * Alignment.h
 *
 *  Created on: Jan 7, 2013
 *      Author: diego
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
