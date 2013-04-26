/*
 * ModelSet.h
 *
 *  Created on: Jan 9, 2013
 *      Author: diego
 */

#ifndef MODELSET_H_
#define MODELSET_H_

#include "Model.h"
#include "ProteicModel.h"
#include "util/Utilities.h"

namespace partest {

class ModelSet {
public:
  ModelSet(bitMask rateVar, DataType dataType, int numberOfTaxa);
  size_t getNumberOfModels() { return numberOfModels; }
  Model * getModel(unsigned int index) {
	  return models[index]; }
  virtual ~ModelSet();
private:
  int buildModelSet(Model **models, bitMask rateVar, DataType dataType);
  size_t numberOfModels;
  Model **models;
  bitMask rateVar;
  DataType dataType;
  int numberOfTaxa;
};

} /* namespace partest */
#endif /* MODELSET_H_ */
