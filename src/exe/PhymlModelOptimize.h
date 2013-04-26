/*
 * PhymlModelOptimize.h
 *
 *  Created on: Jan 10, 2013
 *      Author: diego
 */

#ifndef PHYMLMODELOPTIMIZE_H_
#define PHYMLMODELOPTIMIZE_H_

#include "ModelOptimize.h"
#include "../model/Model.h"
#include "../options/ParTestOptions.h"
#include "../indata/PhymlAlignment.h"
#include "../indata/PartitionElement.h"

namespace partest {

class PhymlModelOptimize : public ModelOptimize {
public:
  PhymlModelOptimize(ParTestOptions * options);
  int optimizeModel(Model * model, PartitionElement * partitionElement, int index, int groupCount);
  virtual ~PhymlModelOptimize();
};

} /* namespace partest */
#endif /* PHYMLMODELOPTIMIZE_H_ */
