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
 * @file ModelOptimize.h
 *
 * @brief Set of functions for optimizing models, partitions and schemes
 */

#ifndef MODELOPTIMIZE_H_
#define MODELOPTIMIZE_H_

#include "util/GlobalDefs.h"
#include "indata/PartitioningScheme.h"
#include "indata/PartitionElement.h"

#include <string>

namespace partest
{

  class ModelOptimize
  {
  public:
    ModelOptimize ();
    virtual ~ModelOptimize ();

    std::string buildStartingFixedTree (int do_ml);
    std::string buildFinalTree (PartitioningScheme * finalScheme,
                                bool reoptimizeParameters);
    std::string buildFinalTreeLinking (PartitioningScheme * finalScheme,
                                       bool reoptimizeParameters);
    int optimizePartitioningScheme (PartitioningScheme * scheme, int index = 0,
                                    int limit = 1);
    int optimizePartitionElement (PartitionElement * scheme, int index = 0,
                                  int limit = 1);
  private:
    void optimizeModel (PartitionElement * element, size_t modelIndex,
                        int limit);
    void setModelParameters (t_partitionElementId id, Model * _model,
                             pllInstance * _tree, partitionList * _partitions,
                             pllAlignmentData * _alignData, int index,
                             bool setAlphaFreqs);
  };

} /* namespace partest */

#endif /* MODELOPTIMIZE_H_ */
