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
 *  For any other inquiries send an Email to Diego Darriba
 *  ddarriba@udc.es
 */

/**
 * @file NucleicModel.h
 *
 * @brief Nucleotide substitution model.
 */

#ifndef NUCLEICMODEL_H_
#define NUCLEICMODEL_H_

#include "Model.h"
#include "util/GlobalDefs.h"

#define NUM_NUC_FREQS 4 /** Number of states for nucleotide substitution models */
#define NUM_RATES 6 /** Number of substitution rates */

namespace partest {

/**
 * @brief Nucleotide substitution model.
 */
class NucleicModel: public Model {
public:
	/**
	 * Constructor of a nucleotide substitution model.
	 *
	 * @param matrix Nucleotide substitution scheme.
	 * @param rateVariation The rate variation and frequencies parameters (+I, +G, +F).
	 * @param numberOfTaxa Number of taxa (required for computing the free parameters.
	 */
	NucleicModel(NucMatrix matrix, bitMask rateVariation, int numberOfTaxa);
	NucMatrix getMatrix(void);
	virtual void setFrequencies(const double * frequencies);
	virtual void allocateRates(void);
	void setRates(const double * rates);
	double distanceTo(Model * other);
	virtual void print(ostream& out, const char * prefix = "");
	virtual ~NucleicModel();
private:
	NucMatrix matrix; /** Nucleotide substitution scheme */
};

} /* namespace partest */
#endif /* NUCLEICMODELL_H_ */
