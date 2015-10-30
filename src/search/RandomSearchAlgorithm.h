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
 * @file RandomSearchAlgorithm.h
 * @author Diego Darriba
 * @brief Algorithm for search among random partitioning schemes
 */

#ifndef RANDOMSEARCHALGORITHM_H_
#define RANDOMSEARCHALGORITHM_H_

#include "SearchAlgorithm.h"
#include "indata/PartitioningScheme.h"

namespace partest
{

  /**
   * @brief Multi-stage search algorithm choosing random states.
   *
   * This algorithms randomly selects a set of states using the Chinese Restaurant
   * Process. On each stage reduces the states space depending on the best-fit
   * partition of the previous state.
   *
   * TODO: This algorithm can be combined with other using the best-fit partitions
   * as fixed starting points.
   */
  class RandomSearchAlgorithm : public SearchAlgorithm
  {
  public:
    RandomSearchAlgorithm ();
    virtual ~RandomSearchAlgorithm ();
    virtual PartitioningScheme * start (PartitioningScheme * startingPoint = 0);
  private:
    int getRandomPartitioningScheme (
        std::vector<PartitioningScheme *> & nextSchemes, int numberOfSchemes,
        t_partitioningScheme p0);
  };

} /* namespace partest */
#endif /* RANDOMSEARCHALGORITHM_H_ */
