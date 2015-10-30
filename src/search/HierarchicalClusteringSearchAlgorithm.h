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
 * @file HierarchicalClusteringSearch.h
 * @author Diego Darriba
 * @brief Algorithm for performing a hierarchical clustering search
 */

#ifndef HIERARCHICALCLUSTERING_H_
#define HIERARCHICALCLUSTERING_H_

#include "search/SearchAlgorithm.h"

namespace partest
{

  class HierarchicalClusteringSearchAlgorithm : public SearchAlgorithm
  {
  public:
    HierarchicalClusteringSearchAlgorithm ();
    virtual ~HierarchicalClusteringSearchAlgorithm ();
    virtual PartitioningScheme * start (PartitioningScheme * startingPoint = 0);
  };

} /* namespace partest */

#endif /* HIERARCHICALCLUSTERING_H_ */
