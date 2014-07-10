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
 * @file PartitioningScheme.h
 *
 * @brief Definition of a Partitioning Scheme
 */

#ifndef PARTITION_H_
#define PARTITION_H_

#include "PartitionElement.h"
#include "util/GlobalDefs.h"
#include <string>
#include <ostream>

#define FULL_CODE -1

namespace partest {

typedef struct {
	PartitionElement * e1;
	PartitionElement * e2;
	double distance;
} elementPair;

/**
 * @brief Definition of a Partitioning Scheme.
 * A partitioning scheme is a set of partitions that covers the
 * whole alignment without overlapping.
 */
class PartitioningScheme {
public:

	/**
	 * @brief Creates a new partitioning scheme.
	 *
	 * @param[in] partition Vector of partition indexes of this scheme.
	 */
	PartitioningScheme(t_partitioningScheme * partition);

	/**
	 * @brief Gets the number of partitions of this scheme.
	 * @return The number of partitions of this scheme.
	 */
	size_t getNumberOfElements() {
		return numberOfElements;
	}

	/**
	 * @brief Gets a locally indexed partition.
	 *
	 * @param[in] id The local id of the partition (i.e., in range [0,numberOfElements-1]
	 * @return The partition.
	 */
	PartitionElement * getElement(size_t index);

	t_partitioningScheme getId(void) { return id; }

	/**
	 * @brief Gets whether all the PartitionElement instances were optimized or not.
	 *
	 * @return true if all the PartitionElement instances were already optimized. false otherwise.
	 */
	bool isOptimized(void);

	void buildCompleteModelSet(bool clearAll = false);

	void setTree(char * tree);

	/**
	 * @brief Gets the optimized tree.
	 *
	 * @return The optimized tree.
	 */
	char * getTree(void) const;

	/**
	 * @brief Gets the closest pair of partitions.
	 *
	 * @param[out] el1, el2 the two closest partitions
	 */
	vector<elementPair *> * getElementDistances();

	/**
	 * @brief Gets the number of lines of the code
	 */
	int getCodeLines(void);

	/**
	 * @brief Gets a string identifier for the scheme.
	 *
	 * @param[in] codeLine If the code has more than 1 line (i.e, more than 10 elements), determines the line to get.
	 *
	 * @return A string identifier for the scheme.
	 */
	string getCode(int codeLine=FULL_CODE);

	/**
	 * @brief Gets a string name for the scheme.
	 *
	 * @return A string identifier for the scheme.
	 */
	string getName();

	virtual ~PartitioningScheme();

	double getLnL();
	double getIcValue();
	double getLinkedBicValue();
	unsigned int getNumberOfFreeParameters();

	void print(ostream & out);
private:
	t_partitioningScheme id;
	vector<PartitionElement*> partitions; /** Array of reference to the partitions of this scheme */
	size_t currentElement; /** Current element index for the step-by-step construction of the scheme */
	size_t numberOfElements; /** The number of partitions */
	vector<elementPair *> * eps;
	char * tree;

	/** Number of lines of the scheme code */
	size_t codeLines;

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
