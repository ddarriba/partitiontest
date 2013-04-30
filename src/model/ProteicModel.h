/*
 * ProteicModel.h
 *
 *  Created on: Jan 8, 2013
 *      Author: diego
 */

#ifndef PROTEICMODEL_H_
#define PROTEICMODEL_H_

#include "Model.h"
#include "../util/GlobalDefs.h"

#define NUM_PROT_FREQS 20

namespace partest {

class ProteicModel: public Model {
public:
  ProteicModel(ProtMatrix matrix, bitMask rateVariation, int numberOfTaxa);
  void setFrequencies(const double * frequencies);
  void setRates(const double * rates);
  virtual ~ProteicModel();
private:
  ProtMatrix matrix;
};

} /* namespace partest */
#endif /* PROTEICMODEL_H_ */
