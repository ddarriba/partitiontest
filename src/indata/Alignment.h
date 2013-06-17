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
 * @brief Generic implementation for MSA definition
 */
#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_

#include <string>
#include "util/GlobalDefs.h"

using namespace std;

namespace partest {

/**
 * @brief Generic implementation for MSA
 *
 * This class defines the interface and implements common methods for input MSAs
 */
class Alignment {
public:
	/**
	 * @brief Void constructor. Instantiates an empty alignment.
	 *
	 * Void constructor. Initializes all members to 0.
	 */
	Alignment(void);
	/**
	 * @brief Instantiates an alignment from a file.
	 *
	 * @param alignmentFile File with the MSA data.
	 * @param datatype Determines whether the MSA is nucleic or proteic data.
	 */
	Alignment(string alignmentFile, DataType dataType);
	virtual ~Alignment();
	/**
	 * @brief Constructs an alignment with a sequential subset of the sites.
	 *
	 * Creates a new alignment by splitting the current one from the specified
	 * nucleotide/amino-acid positions.
	 *
	 * @param[in] firstPosition The position of the first site in the subset.
	 * @param[in] lastPosition The position of the last site in the subset.
	 * @return The new alignment
	 */
	virtual Alignment * splitAlignment(int firstPosition, int lastPosition);
	/**
	 * @brief Constructs an alignment from the concatenation of one or more subsets of the sites.
	 *
	 * Creates a new alignment by splitting the current one by merging one or more
	 * sets of the MSA sites.
	 *
	 * @param[in] firstPosition The start positions of each subset.
	 * @param[in] lastPosition The end positions of each subset.
	 * @param[in] numberOfSections The number of subsets and also the length of the previous arrays.
	 * @return The new alignment
	 */
	virtual Alignment * splitAlignment(int * firstPosition, int * lastPosition,
			int numberOfSections);

	/**
	 * @brief Gets the MSA filename.
	 * @return The MSA filename.
	 */
	string getAlignmentFile() const {
		return alignmentFile;
	}

	/**
	 * @brief Gets the number of sites in the MSA.
	 * @return The number of sites
	 */
	unsigned int getNumSites() const {
		return numSites;
	}

	/**
	 * @brief Gets the number of sequences in the MSA.
	 * @return The number of sequences or taxa
	 */
	unsigned int getNumSeqs() const {
		return numSeqs;
	}

	/**
	 * @brief Gets the number of unique site patterns in the MSA.
	 *
	 * Gets the number of unique site patterns in the MSA. If there is
	 * no pattern compression, this equals to the number of sites.
	 *
	 * @return The number of unique site patterns.
	 */
	unsigned int getNumPatterns() const {
		return numPatterns;
	}

	/**
	 * @brief Gets the Shannon Entropy of the alignment.
	 *
	 * @return The Shannon Entropy of the alignment.
	 */
	double getShannonEntropy() const {
		return shannonEntropy;
	}
protected:
	/**
	 * @brief Computes the Shannon Entropy of the alignment.
	 *
	 * Shannon entropy is the average unpredictability in a
	 * random variable, which is equivalent to its information content.
	 * This value is calculated as follows:
	 *   \f[
	 H(X) = \sum_{i} {P(x_i)\,I(x_i)} = -\sum_{i} {P(x_i) \log_b P(x_i)} = -\sum_{i}{\frac{n_i}{N} \log_b \frac{n_i}{N}} = \log_b N - \frac 1 N \sum_{i} {n_i \log_b n_i},
	 \f]
	 * where b is the base of the logarithm used.
	 *
	 * @param observedFrequencies The empirical state frequencies in the alignment.
	 * @return The Shannon Entropy of the alignment.
	 */
	double computeShannonEntropy(double *observedFrequencies);
	string alignmentFile; /** Alignment input file. */
	unsigned int numSites; /** Number of sites in each sequence. */
	unsigned int numPatterns; /** Number of different site patterns. */
	unsigned int numSeqs; /** Number of sequences. */
	double *empiricalFrequencies; /** Observed frequencies of nucleotides/amino-acids. */
	int numberOfFrequencies; /** 4 for nucleotides, 20 for amino-acids. */
	DataType dataType; /** Nucleotide or amino-acid data. */
	/** @brief Shannon Entropy of the alignment.
	 * Shannon entropy is the average unpredictability in a
	 * random variable, which is equivalent to its information content.
	 * */
	double shannonEntropy;
};

} /* namespace partest */

#endif /* ALIGNMENT_H_ */
