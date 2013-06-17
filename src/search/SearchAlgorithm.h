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
 * @file SearchAlgorithm.h
 */

#ifndef SEARCHALGORITHM_H_
#define SEARCHALGORITHM_H_

#include "indata/PartitionMap.h"
#include "indata/PartitioningScheme.h"
#include "indata/PartitionElement.h"
#include "model/ModelSet.h"
#include "options/ParTestOptions.h"
#include "util/GlobalDefs.h"
#include "indata/Alignment.h"
#include "observer/Observer.h"
#include "observer/Observable.h"

namespace partest {

/**
 * @brief Interface for search algorithms.
 *
 * Generic interface for search algorithms.
 */
class SearchAlgorithm: public Observer, public Observable {
public:

	/**
	 * @brief Creates a new search algorithm.
	 *
	 * @param options Execution options.
	 * @param partitionMap Partitions map instance.
	 */
	SearchAlgorithm(ParTestOptions * options, PartitionMap * partitionMap);
	virtual ~SearchAlgorithm();

	/**
	 * @brief Starts a new search through the complete states space.
	 *
	 * @return The best-fit partitioning scheme.
	 */
	virtual PartitioningScheme * start() = 0;

	/**
	 * @brief Starts a new search starting from a fixed point.
	 *
	 * @param startingScheme The starting point for searching.
	 *
	 * @return The best-fit partitioning scheme.
	 */
	virtual PartitioningScheme * start(PartitioningScheme * startingPoint) = 0;

	virtual void update(const ObservableInfo & info,
			ParTestOptions * run_instance = NULL) = 0;
protected:

	/**
	 * @brief Gets a PartitionElement from its id.
	 *
	 * @parameter id PartitionElement's id.
	 *
	 * @return The partition element.
	 */
	PartitionElement * getPartitionElement(t_partitionElementId id);


	/**
	 * @brief Gets the number of mapped elements.
	 *
	 * @return The number of elements in the partitions map.
	 */
	unsigned int getNumberOfElements();

	ParTestOptions * options; /** Execution options */
	PartitionMap * partitionMap; /** Map of partitions */
};

} /* namespace partest */
#endif /* SEARCHALGORITHM_H_ */
