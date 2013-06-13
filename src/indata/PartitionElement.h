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
 * @file PartitionElement.h
 *
 * @brief Definition of a gene partition
 */

#ifndef PARTITION_ELEMENT_H_
#define PARTITION_ELEMENT_H_

#include "model/ModelSet.h"
#include "util/GlobalDefs.h"
#include "indata/Alignment.h"
#include "selection/SelectionModel.h"

namespace partest {

/**
 * @brief Definition of a gene partition.
 *
 * A partition is a subset of the multi-gene MSA that can contain one or more full genes.
 */
class PartitionElement {
public:
	/**
	 * @brief Constructs a new partition for a single gene.
	 *
	 * Constructs a new partition for a single gene and creates a set of candidate models.
	 *
	 * @param[in] id Unique id for the partition. Single-gene partitions have a power of 2 id.
	 * @param[in] name Name of the gene if this is a single-gene partition.
	 * @param[in] alignment The alignment this partition belongs to.
	 * @param[in] start Position of the first site of the partition in the whole alignment.
	 * @param[in] end Position of the last site of the partition in the whole alignment.
	 * @param[in] stride [1,3] if this partition uses only a codon position. 0 otherwise.
	 * @param[in] rateVariation The rate variations mask to be evaluated.
	 * @param[in] dataType Whether the data is nucleic or proteic.
	 */
	PartitionElement(t_partitionElementId id, string name,
			Alignment * alignment, int start, int end, int stride,
			bitMask rateVariation, DataType dataType);

	/**
	 * @brief Constructs a new partition for a subset of genes.
	 *
	 * Constructs a new partition for a subset of genes and creates a set of candidate models.
	 *
	 * @param[in] id Unique id for the partition. Single-gene partitions have a power of 2 id.
	 * @param[in] name Name of the gene if this is a single-gene partition.
	 * @param[in] alignment The alignment this partition belongs to.
	 * @param[in] start Positions of the first site of the subpartitions in the whole alignment.
	 * @param[in] end Positions of the last site of the subpartitions in the whole alignment.
	 * @param[in] stride [1,3] if the partitions use only a codon position. 0 otherwise.
	 * @param[in] numberOfElements Number of subpartitions to be joined.
	 * @param[in] rateVariation The rate variations mask to be evaluated.
	 * @param[in] dataType Whether the data is nucleic or proteic.
	 */
	PartitionElement(t_partitionElementId id, string name,
			Alignment * alignment, int * start, int * end, int * stride,
			int numberOfSections, bitMask rateVariation, DataType dataType);

	/**
	 * @brief Gets the set of candidate models.
	 *
	 * @return The set of candidate models.
	 */
	ModelSet * getModelset(void) {
		return modelset;
	}

	/**
	 * @brief Gets the alignment for this partition.
	 *
	 * @return The alignment for this partition.
	 */
	Alignment * getAlignment(void) {
		return alignment;
	}

	/**
	 * @brief Gets a subpartition of sequential sites of this.
	 *
	 * @param[in] first Position of the first site of the partition in the whole alignment.
	 * @param[in] last Position of the last site of the partition in the whole alignment.
	 * @return The new subpartition.
	 */
	PartitionElement * splitPartition(int first, int last);

	/**
	 * @brief Gets whether the candidate models were optimized or not.
	 *
	 * @return 1 if the models were already optimized. 0 otherwise.
	 */
	bool isOptimized(void);

	void buildCompleteModelSet(bool clearAll = false);

	virtual ~PartitionElement();

	/**
	 * @brief Gets the unique identifier of the partition.
	 *
	 * @return The unique identifier of the partition.
	 */
	t_partitionElementId getId(void) {
		return id;
	}

	/**
	 * @brief Gets the higher index in the partition.
	 *
	 * @return The higher index in the partition.
	 */
	unsigned int getMaxId(void) {
		return id.at(id.size() - 1);
	}

	/**
	 * @brief Gets the lower index in the partition.
	 *
	 * @return The lower index in the partition.
	 */
	unsigned int getMinId(void) {
		return id.at(0);
	}

	/**
	 * @brief Gets the global reference to the start position of a section.
	 *
	 * @section The section of the partition. This value should be in [0, numberOfSections-1]
	 *
	 * @return The starting site of the section.
	 */
	int getStart(int section = 0);

	/**
	 * @brief Gets the global reference to the end position of a section.
	 *
	 * @section The section of the partition. This value should be in [0, numberOfSections-1]
	 *
	 * @return The ending site of the section.
	 */
	int getEnd(int section = 0);

	/**
	 * @brief Gets the stride of a section if the section is divide by a codon position.
	 *
	 * @section The section of the partition. This value should be in [0, numberOfSections-1]
	 *
	 * @return The stride or codon position [1,3]. 0, if there is no subdivision.
	 */
	int getStride(int section = 0);

	/**
	 * @brief Gets the number of sections (genes or gene codon positions) of the partition.
	 *
	 * @return The number of sections of the partition.
	 */
	int getNumberOfSections(void);

	/**
	 * @brief Gets the unique name of this partition.
	 *
	 * @return The name of the partition.
	 */
	string getName(void);

	/**
	 * @brief Gets the best model according to the criterion in the global options.
	 *
	 * After evaluating the candidate models, and selecting the best one,
	 * this member function returns the best-fit model for this partition.
	 *
	 * @return The best-fit model for this partition.
	 */
	SelectionModel * getBestModel(void);

	/**
	 * @brief Sets the best model according to the criterion in the global options.
	 *
	 * @param[in] bestModel The best-fit model for this partition.
	 */
	void setBestModel(SelectionModel * bestModel);

private:
	/**
	 * @brief Unique identifier of this partition.
	 *
	 * Unique identifier of this partition. Using this identifier the user and the
	 * program knows what genes are inside this partition. Single genes have a power
	 * of 2 id, and the partition id is the sum of all the contained genes.
	 *
	 * In base 2, this can be seen as each bit position representing whether the
	 * gene is present or not in the partition.
	 *
	 * e.g, 001101 = 13 means this partition contains genes 0,2 and 3.
	 */
	t_partitionElementId id;
	Alignment * alignment; /** MSA for this partition */
	int * start; /** The start positions of each subset. */
	int * end; /** The end positions of each subset. */
	int * stride; /** The stride or codon position [1,3]. 0, if there is no subdivision. */
	int numberOfSections; /** The number of subsets and also the length of the previous arrays. */
	string name; /** The name of the partition. */
	SelectionModel * bestModel; /** The best-fit model for this partition. */
	ModelSet * modelset; /** The set of candidate models. */
};

} /* namespace partest */
#endif /* PARTITION_ELEMENT_H_ */
