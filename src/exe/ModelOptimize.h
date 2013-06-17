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
 * @file ModelOptimize.h
 *
 * @brief Generic interface for optimizing models
 */
#ifndef MODELOPTIMIZE_H_
#define MODELOPTIMIZE_H_

#include "model/Model.h"
#include "indata/PartitioningScheme.h"
#include "indata/PartitionElement.h"
#include "options/ParTestOptions.h"
#include "observer/Observable.h"

namespace partest {

/**
 * @brief Generic interface for optimizing models
 */
class ModelOptimize: public Observable {
public:
	/**
	 * @brief Creates a new interface for Model Optimizers.
	 *
	 * @param[in]	options PartitionTest Options for model optimization.
	 */
	ModelOptimize(ParTestOptions * options);

	/**
	 * @brief Optimizes the parameters of the partitions that belong to a partitioning scheme.
	 *
	 * @param[in]	scheme Partitioning scheme to be evaluated.
	 * @param[in]	forceRecomputation If true, the scheme is evaluated even if it was already evaluated before.
	 */
	virtual int optimizePartitioningScheme(PartitioningScheme * scheme,
			bool forceRecomputation = false, int current_index = 0, int max_index = 0);

	/**
	 * Optimizes the parameters of a single partition.
	 *
	 * @param[in]	partitionElement Partition to be evaluated
	 * @return		0 if the evaluation ended successfully. Error code otherwise.
	 */
	virtual int optimizePartitionElement(PartitionElement * partitionElement, int current_index = 0, int max_index = 0);

	/**
	 * @brief Optimizes the parameters of a single model.
	 *
	 * @param[in,out]	model Model to be evaluated.
	 * @param[in]		partitionElement Partition the model belongs to.
	 * @param[in]		index The index of the model in the modelset.
	 * @param[in]		groupCount The number of candidate models.
	 *
	 * @return		0 if the evaluation ended successfully. Error code otherwise.
	 */
	virtual int optimizeModel(Model * model,
			PartitionElement * partitionElement, int index, int groupCount) = 0;
	virtual ~ModelOptimize();
protected:
	ParTestOptions * options; /** PartitionTest options for model optimization. */
};

} /* namespace partest */
#endif /* MODELOPTIMIZE_H_ */
