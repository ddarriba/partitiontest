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
 * @file ProteicModel.h
 *
 * @brief Amino-acid replacement model.
 */

#ifndef PROTEICMODEL_H_
#define PROTEICMODEL_H_

#include "Model.h"
#include "util/GlobalDefs.h"

#define NUM_PROT_FREQS 20 /** Number of states for amino-acid replacement models */

namespace partest {

/**
 * @brief Amino-acid replacement model.
 */
class ProteicModel: public Model {
public:
	/**
	 * @brief Constructor of a amino-acid replacement model.
	 *
	 * @param matrix Amino-acid replacement matrix.
	 * @param rateVariation The rate variation and frequencies parameters (+I, +G, +F).
	 * @param numberOfTaxa Number of taxa (required for computing the free parameters.
	 */
  ProteicModel(ProtMatrix matrix, bitMask rateVariation, int numberOfTaxa);
  virtual void setFrequencies(const double * frequencies);
  virtual void allocateRates(void) { /* do nothing */ }
  void setRates(const double * rates);
  ProtMatrix getMatrix(void);
  void setMatrix(ProtMatrix  matrix);
  double distanceTo(Model * other);
  virtual void print(ostream& out, const char * prefix = "");
  virtual ~ProteicModel();
private:
  /**
   * @brief Compute euclidean distances between amino-acid replacement matrices
   */
  double getEuclideanDistance(ProtMatrix m1, ProtMatrix m2);
  ProtMatrix matrix; /** Amino-acid replacement matrix. */
};

} /* namespace partest */
#endif /* PROTEICMODEL_H_ */
