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
 * @file PLLAlignment.h
 *
 * @brief Implementation of the MSA for the Phylogenetic Likelihood Library.
 */
#ifndef PLLALIGNMENT_H_
#define PLLALIGNMENT_H_

#include "Alignment.h"

extern "C" {
#include "parser/phylip/phylip.h"
#include "parser/partition/part.h"
//#include "parser/phylip/phylip.h"
//#include "utils.h"
//struct pllPhylip * pllPhylipParse (const char *);
}

using namespace std;
namespace partest {

/**
 * @brief Implementation of the MSA for the Phylogenetic Likelihood Library.
 */
class PLLAlignment: public Alignment {
public:
	PLLAlignment(PLLAlignment * alignment, int firstPosition, int lastPosition);
	PLLAlignment(PLLAlignment * alignment, int * firstPosition,
			int * lastPosition, int numberOfSections);
	PLLAlignment(string alignmentFile, DataType dataType,
			string partitionsFile);
	virtual ~PLLAlignment();
	Alignment * splitAlignment(int firstPosition, int lastPosition);
	Alignment * splitAlignment(int * firstPosition, int * lastPosition,
			int numberOfSections);
	/**
	 * @brief Gets the tree data structure from PLL.
	 *
	 * @return The tree data structure.
	 */
	pllInstance * getTree() {
		return tr;
	}
	pllPhylip * getPhylip() {
		return phylip;
	}
	partitionList * getPartitions() {
		return partitions;
	}
private:
	struct pllQueue * parts; /** Single parsed partitions */
	pllInstance * tr; /** Tree structure for working with PLL. */
	partitionList * partitions; /** Partition definitions */
	struct pllPhylip * phylip; /** Alignment */
	void buildTopology(void);
};

} /* namespace partest */
#endif /* PLLAlignment_H_ */
