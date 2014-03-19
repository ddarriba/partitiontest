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
 * @file PartitionMap.h
 *
 * @brief Mapping of the simplest and unbreakable partitions.
 */

#ifndef PARTITIONMAP_H_
#define PARTITIONMAP_H_

#include <stdlib.h>
#include <vector>
#include "PartitionElement.h"
#include "parser/ConfigParser.h"
#include "util/GlobalDefs.h"
#include "options/ParTestOptions.h"

namespace partest {

struct partitionMappingInfo {
	partitionMappingInfo() :
			partitionId(1), partitionElement(0) {
	}
	t_partitionElementId partitionId;
	PartitionElement * partitionElement;
	~partitionMappingInfo(void) {
	} //delete partitionElement; }
};

/**
 * @brief Mapping of the simplest and unbreakable partitions.
 *
 * This class represent the mapping of the simplest and unbreakable partitions.
 * Each partition has a power of two Id, what allows to easy identify a complex partition from the id.
 * e.g., Partition 13 (1101) will be a merge of partitions 8 (1000),4 (0100) and 1 (0001).
 */
class PartitionMap {
public:

	/**
	 * @brief Constructs a new partition map.
	 *
	 * @param[in] alignment The MSA containing all the genes.
	 * @param[in] numberOfPartitions The number of minimum partitions (i.e, single genes)
	 * @param[in] rateVariation The rate variations mask to be evaluated.
	 * @param[in] dataType Whether the data is nucleic or proteic.
	 */
	PartitionMap(Alignment * alignment, unsigned int numberOfPartitions,
			bitMask rateVariation, DataType dataType, OptimizeMode optimizeMode);

	/**
	 * @brief Constructs a new partition map.
	 *
	 * @param[in] configFile Configuration file containing the description of all the partitions.
	 * @param[in] alignment The MSA containing all the genes.
	 * @param[in] rateVariation The rate variations mask to be evaluated.
	 * @param[in] dataType Whether the data is nucleic or proteic.
	 */
	PartitionMap(const char * configFile, Alignment * alignment,
			bitMask rateVariation, DataType dataType, OptimizeMode optimizeMode);

	virtual ~PartitionMap();

	/**
	 * @brief Adds a new single-gene partition to the map.
	 *
	 * Adds a new single-gene partition to the map. All this partitions will be mapped to positions
	 * with identifier power of 2. For example, the partition with partitionId 4 will be mapped to
	 * the position 16. This allows to create new partitions were each bit of the identifier means
	 * whether the single-gene partition in that position is present or not.
	 *
	 * @param[in] partitionId Sequential identifier of the partition.
	 * @param[in] iniPosition Position of the first site of the partition in the whole alignment.
	 * @param[in] endPosition Position of the last site of the partition in the whole alignment.
	 * @param[in] stride [1,3] if this partition uses only a codon position. 0 otherwise.
	 *
	 * @return true, if the partition was successfully added to the map.
	 */
	bool addPartitionElement(unsigned int partitionId, string name,
			unsigned int iniPosition, unsigned int endPosition, char stride);
	/**
	 * @brief Gets the partition in a fixed position.
	 *
	 * @param[in] partitionId Position of the partition.
	 *
	 * @return The partition.
	 */
	PartitionElement * getPartitionElement(t_partitionElementId partitionId);

	/**
	 * @brief Gets the partition in a fixed position.
	 *
	 * @param[in] partitionId Position of the partition.
	 *
	 * @return The partition.
	 */
	PartitionElement * getPartitionElement(unsigned int id);

	/**
	 * @brief Deletes a PartitionElement
	 *
	 * @parameter id PartitionElement's id.
	 */
	void deletePartitionElement(t_partitionElementId id);

	/**
	 * @brief Gets the number of elements already created in the map.
	 *
	 * Gets the number of elements already created in the map. For a set of N genes,
	 * this value belongs to the interval [N,2^N]
	 *
	 * @return The number of elements already created in the map.
	 */
	unsigned int getNumberOfElements() {
		return numberOfElements;
	}

	/**
	 * @brief Gets the number of single-gene partitions.
	 *
	 * @return The number of single-gene partitions.
	 */
	unsigned int getNumberOfPartitions() {
		return numberOfPartitions;
	}

	void purgePartitionMap(t_partitionElementId id);

private:
	Alignment * alignment; /** MSA with the information of all the genes. */
	unsigned int numberOfElements; /** Number of partitions already created in the map. */
	unsigned int numberOfPartitions; /** Number of single-gene partitions */
	vector<partitionMappingInfo> * partitions; /** Vector containing all the partitions created. */
	bitMask rateVariation; /** The rate variations mask to be evaluated. */
	DataType dataType; /** Whether the data is nucleic or proteic. */
	OptimizeMode optimizeMode;
};

} /* namespace partest */
#endif /* PARTITIONMAP_H_ */
