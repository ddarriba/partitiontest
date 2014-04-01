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
 * @file ModelSet.h
 */
#ifndef MODELSET_H_
#define MODELSET_H_

#include "ProteicModel.h"
#include "NucleicModel.h"
#include "util/Utilities.h"
#include "util/GlobalDefs.h"

namespace partest {

/**
 * @brief Set of substitution models for the same data.
 */
class ModelSet {
public:
	/**
	 * @brief Constructor for the model set
	 *
	 * This constructor creates a model set with all the possible
	 * combinations of the rate variation and frequencies parameters.
	 *
	 * e.g, If rateVar equals +I|+F, there will be 4 times the number of
	 * matrices with no rate variation nor frequencies parameters,
	 * +F, +I and +I+F parameters.
	 *
	 * @param[in] rateVar Set of rate variations and frequencies (+I, +G, +F).
	 * @param[in] dataType Whether the data is nucleotide or amino-acids.
	 * @param[in] numberOfTaxa The number of taxa in the alignment.
	 */
  ModelSet(bitMask rateVar, DataType dataType, int numberOfTaxa, OptimizeMode optimizeMode = OPT_SEARCH, bool forceCompleteSet = false, bool emptySet=false);

  /**
   * @brief Gets the total number of models.
   *
   * @return The total number of models in the set.
   */
  size_t getNumberOfModels() { return numberOfModels; }

  void buildCompleteModelSet(bool clearAll = false);

  /**
   * @brief Gets the model with the index.
   *
   * @param[in] index Index of the model.
   *
   * @return The model in the index position.
   */
  Model * getModel(unsigned int index);

  /**
     * @brief Sets the model in the index.
     *
     * @param[in] model The model.
     * @param[in] index Index of the model.
     */
    void setModel(Model * model, unsigned int index);

    void allocateModels(int numberOfModels);
  DataType getDataType(void) { return dataType; }
  virtual ~ModelSet();
private:

  /**
   * @brief Constructs a set of models with certain parameters.
   *
   * @param[out] models The set of models with the specified parameters.
   * @param[in] rateVar The set of parameters for constructing the models.
   * @param[in] forceCompleteSet If FAST DNA SEARCH is active, this flag forces to build the complete candidate model set
   *
   * @return 0, if there was no errors.
   */
  int buildModelSet(Model **models, bitMask rateVar, bool forceCompleteSet = false);

  size_t numberOfModels; /** Number of models in the set. */
  Model **models; /** Array of models. */
  bitMask rateVar; /** Set of rate variation and frequencies parameters. */
  DataType dataType; /** Whether the data is nucleotide or amino-acids. */
  int numberOfTaxa; /** Number of taxa in the input data. */
  OptimizeMode optimizeMode;
};

} /* namespace partest */
#endif /* MODELSET_H_ */
