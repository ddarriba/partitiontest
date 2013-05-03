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
 * @file PhymlModelOptimize.h
 *
 * @brief Implementation of model optimization using PhyML
 */
#ifndef PHYMLMODELOPTIMIZE_H_
#define PHYMLMODELOPTIMIZE_H_

#include "ModelOptimize.h"
#include "model/Model.h"
#include "options/ParTestOptions.h"
#include "indata/PhymlAlignment.h"
#include "indata/PartitionElement.h"

namespace partest {

/**
 * @brief Implementation of model optimization using the PhyML
 */
class PhymlModelOptimize : public ModelOptimize {
public:
  PhymlModelOptimize(ParTestOptions * options);
  int optimizeModel(Model * model, PartitionElement * partitionElement, int index, int groupCount);
  virtual ~PhymlModelOptimize();
};

} /* namespace partest */
#endif /* PHYMLMODELOPTIMIZE_H_ */
