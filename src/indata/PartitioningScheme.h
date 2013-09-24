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
 * @file PartitioningScheme.h
 *
 * @brief Definition of a Partitioning Scheme
 */
#ifndef PARTITION_H_
#define PARTITION_H_

#include "PartitionElement.h"
#include "PartitionMap.h"
#include "util/GlobalDefs.h"
#include <string>

namespace partest {

/**
 * @brief Definition of a Partitioning Scheme.
 * A partitioning scheme is a set of partitions that covers the
 * whole alignment without overlapping.
 */
class PartitioningScheme {
public:

	/**
	 * @brief Creates an empty partitioning scheme.
	 *
	 * Creates an empty partitioning scheme. The number of partitions must be declared.
	 * The new scheme is not valid until all the inner partitions are declared using
	 * the PartitioningScheme::addElement member function.
	 * @param[in] numberOfElements The number of partitions of the scheme.
	 */
	PartitioningScheme(int numberOfElements);

	/**
	 * @brief Creates a new partitioning scheme.
	 *
	 * @param[in] partition Vector of partition indexes of this scheme.
	 * @param[in] partitionMap Reference to the map of partitions.
	 */
	PartitioningScheme(t_partitioningScheme * partition,
			PartitionMap * partitionMap);
	/**
	 * @brief Gets the number of partitions of this scheme.
	 * @return The number of partitions of this scheme.
	 */
	int getNumberOfElements() {
		return numberOfElements;
	}

	/**
	 * @brief Gets the number of single genes in the scheme.
	 * @return The number of single genes.
	 */
	int getNumberOfBits() {
		return numberOfBits;
	}

	/**
	 * @brief Adds a new partition to the scheme.
	 *
	 * Adds a new partition to the scheme. This method works only when the scheme
	 * definition is not complete (i.e., when the scheme was created with the
	 * PartitioningScheme::PartitioningScheme(int numberOfElements) constructor and
	 * the number of defined partitions is less than the declared number of partitions).
	 *
	 * @param[in] element The partition to be added.
	 * @return 0, if successfully added the partition.
	 */
	int addElement(PartitionElement * element);

	/**
	 * @brief Gets a locally indexed partition.
	 *
	 * @param[in] id The local id of the partition (i.e., in range [0,numberOfElements-1]
	 * @return The partition.
	 */
	PartitionElement * getElement(int id);

	/**
	 * @brief Gets whether all the PartitionElement instances were optimized or not.
	 *
	 * @return true if all the PartitionElement instances were already optimized. false otherwise.
	 */
	bool isOptimized(void);

	void buildCompleteModelSet(bool clearAll = false);

	void resetModelSet();

	/**
	 * @brief Gets the closest pair of partitions.
	 *
	 * @param[out] el1, el2 the two closest partitions
	 */
	void getClosestPartitions(t_partitionElementId & el1,
			t_partitionElementId & el2);

	/**
	 * @brief Gets a string identifier for the scheme.
	 *
	 * @return A string identifier for the scheme.
	 */
	string getCode();

	/**
	 * @brief Gets a string name for the scheme.
	 *
	 * @return A string identifier for the scheme.
	 */
	string getName();

	virtual ~PartitioningScheme();

	double getLnL();
private:
	vector<PartitionElement*> * partitions; /** Array of reference to the partitions of this scheme */
	int currentElement; /** Current element index for the step-by-step construction of the scheme */
	int numberOfElements; /** The number of partitions */
	int numberOfBits; /** The number of single genes in the scheme. */
	/**
	 * @brief String identifier of this scheme.
	 *
	 * A string identifier shows how the single genes are grouped into the partitions
	 * that belong to the scheme. For example, 121123 means that genes 1,3 and 4 belong
	 * to the same partition, as well as genes 2 and 5. Gene 6 belongs to a third partition.
	 */
	string * code;
};

} /* namespace partest */
#endif /* PARTITION_H_ */
