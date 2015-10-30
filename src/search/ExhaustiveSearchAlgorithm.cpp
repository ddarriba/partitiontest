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
 * @file ExhaustiveSearchAlgorithm.cpp
 */

#include "ExhaustiveSearchAlgorithm.h"

#include "util/GlobalDefs.h"
#include "exe/ModelOptimize.h"
#include "exe/ModelSelector.h"
#include "exe/PartitionSelector.h"
#include "indata/PartitioningScheme.h"

#include <algorithm>    /* std::adjacent_find */
#include <iostream>
#include <vector>

namespace partest
{

  ExhaustiveSearchAlgorithm::ExhaustiveSearchAlgorithm ()
  {
  }

  ExhaustiveSearchAlgorithm::~ExhaustiveSearchAlgorithm ()
  {
  }

  PartitioningScheme * ExhaustiveSearchAlgorithm::start (
      PartitioningScheme * startingPoint)
  {

    PartitioningScheme *bestScheme = 0;
    if (startingPoint)
    {
      cerr << "[ERROR] Not implemented yet" << endl;
      exit_partest (EX_UNAVAILABLE);
    }
    else
    {
      ModelOptimize * modelOptimize = new ModelOptimize ();

      if (number_of_schemes > 0)
      {
        vector<PartitioningScheme *> candidateSchemes (number_of_schemes);
        //double bestScore;

        for (size_t currentStep = 0; currentStep < number_of_schemes;
            currentStep++)
        {
          t_partitioningScheme scheme = schemes->at (currentStep);
          candidateSchemes.at (currentStep) = new PartitioningScheme (&scheme);
          cout << timestamp () << " [EXH] Step " << currentStep + 1 << "/"
              << number_of_schemes << endl;
          modelOptimize->optimizePartitioningScheme (
              candidateSchemes.at (currentStep));
        }
        PartitionSelector ps (candidateSchemes);
        bestScheme = ps.getBestScheme ();
        //bestScore = ps.getBestScheme()->getIcValue();
      }
      else
      {
        //TODO: UNINPLEMENTED
        cerr << "TO BE IMPLEMENTED" << endl;
        exit_partest (EX_UNAVAILABLE);
      }

    }
    return bestScheme;

  }

} /* namespace partest */
