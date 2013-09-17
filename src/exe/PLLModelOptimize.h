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
#include "pll.h"
#endif

namespace partest {

/**
 * @brief Implementation of model optimization using the Phylogenetic Likelihood Library
 */
class PLLModelOptimize: public ModelOptimize {
public:
	PLLModelOptimize(ParTestOptions * options);
	static void initializeStructs(pllInstance * tree, partitionList * partitions,
			pllAlignmentData * phylip);
	double evaluateSPR(pllInstance * tr, partitionList *pr,
			bool estimateModel = true, bool estimateTopology = true);
	double evaluateNNI(pllInstance * tr, partitionList *pr,
			bool estimateModel = true);
	double optimizeParameters(pllInstance * tr,
			partitionList *partitions, bool estimateModel, bool estimateTopology);
	int optimizePartitioningScheme(PartitioningScheme * scheme,
			bool forceRecomputation = false, int current_index = 0,
			int max_index = 0);
	int optimizePartitioningSchemeAtOnce(PartitioningScheme * scheme);
	int optimizeModel(Model * model, PartitionElement * partitionElement,
			int index, int groupCount);
	virtual ~PLLModelOptimize();
private:
	PLLAlignment * alignment;
	pllInstance * tr; /** Tree structure from PLL */
};

} /* namespace partest */
#endif /* PLLMODELOPTIMIZE_H_ */
