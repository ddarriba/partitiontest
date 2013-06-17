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
 * @file PhymlAlignment.h
 *
 * @brief Implementation of the MSA for PhyML.
 */
#ifndef PHYMLALIGNMENT_H_
#define PHYMLALIGNMENT_H_

#include "Alignment.h"
extern "C" {
typedef double phydbl;
typedef struct __Calign {
	struct __Align **c_seq; /*! compressed sequences      */
	phydbl *b_frq; /*! observed state frequencies */
	short int *invar; /*! 1 -> states are identical, 0 states vary */
	int *wght; /*! # of each site in c_align */
	short int *ambigu; /*! ambigu[i]=1 is one or more of the sequences at site
	 i display an ambiguous character */
	phydbl obs_pinvar;
	int n_otu; /*! number of taxa */
	int clean_len; /*! uncrunched sequences lenghts without gaps */
	int crunch_len; /*! crunched sequences lengths */
	int init_len; /*! length of the uncompressed sequences */
	int *sitepatt; /*! this array maps the position of the patterns in the
	 compressed alignment to the positions in the uncompressed
	 one */
	int format; /*! 0 (default): PHYLIP. 1: NEXUS. */
};

typedef struct __Align {
	char *name; /*! sequence name */
	int len; /*! sequence length */
	char *state; /*! sequence itself */
	short int *is_ambigu; /*! is_ambigu[site] = 1 if state[site] is an ambiguous character. 0 otherwise */
} align;

}
namespace partest {

/**
 * @brief Implementation of the MSA for PhyML.
 */
class PhymlAlignment: public Alignment {
public:
	PhymlAlignment(PhymlAlignment * alignment, int firstPosition,
			int lastPosition);
	PhymlAlignment(PhymlAlignment * alignment, int * firstPosition,
			int * lastPosition, int numberOfSections);
	PhymlAlignment(std::string alignmentFile, DataType dataType);
	virtual ~PhymlAlignment();
	/**
	 * @brief Gets the PhyML data structure.
	 *
	 * @return  A reference to the PhyML data structure.
	 */
	struct __Calign *getCData();
	Alignment * splitAlignment(int firstPosition, int lastPosition);
	Alignment * splitAlignment(int * firstPosition, int * lastPosition,
			int numberOfSections);
private:
	/**
	 * @brief Constructs the PhyML input data structure for the MSA or a subset.
	 *
	 * @param[in] alignment The base MSA.
	 * @param[in] firstPosition The start positions of each subset.
	 * @param[in] lastPosition The end positions of each subset.
	 * @param[in] numberOfSections The number of subsets and also the length of the previous arrays.
	 *
	 * @return  A reference to the PhyML data structure.
	 */
	struct __Calign * build_cdata(PhymlAlignment * alignment,
			int * firstPosition, int * lastPosition, int numberOfSections);
	struct __Calign * cdata; /** PhyML input data structure */
};

} /* namespace partest */
#endif /* PHYMLALIGNMENT_H_ */
