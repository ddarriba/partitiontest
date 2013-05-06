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

/*!
 * @file PLLModelOptimize.h
 *
 * @brief Implementation of model optimization using the Phylogenetic Likelihood Library
 */

#ifndef PLLMODELOPTIMIZE_H_
#define PLLMODELOPTIMIZE_H_

#include "ModelOptimize.h"
#include "model/ModelSet.h"

/* external C */
#ifndef AXML_H
#define AXML_H
#include "axml.h"
#endif
//#include "parser/phylip.h"
//extern "C" {
//void read_phylip_msa (tree * tr, const char * filename, int format, int type);
//}

namespace partest {

/**
 * @brief Implementation of model optimization using the Phylogenetic Likelihood Library
 */
class PLLModelOptimize: public ModelOptimize {
public:
	PLLModelOptimize(ParTestOptions * options);
	int optimizeModel(Model * model, PartitionElement * partitionElement,
			int index, int groupCount);
	virtual ~PLLModelOptimize();
private:
	tree * tr; /** Tree structure from PLL */
	double **empiricalFrequencies; /** Empirical frequencies of the site states. */
};

} /* namespace partest */
#endif /* PLLMODELOPTIMIZE_H_ */
