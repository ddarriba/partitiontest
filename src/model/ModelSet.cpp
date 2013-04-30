/*
 * ModelSet.cpp
 *
 *  Created on: Jan 9, 2013
 *      Author: diego
 */

#include "ModelSet.h"

namespace partest {

ModelSet::ModelSet(bitMask rateVar, DataType dataType, int numberOfTaxa) :
    rateVar(rateVar), dataType(dataType), numberOfTaxa(numberOfTaxa)
{
  int numberOfParameters = Utilities::setbitsCount(rateVar>>1);
  int numberOfMatrices;
  numberOfModels = Utilities::binaryPow(numberOfParameters);
  switch (dataType) {
  case DT_NUCLEIC:
    numberOfMatrices = NUC_MATRIX_SIZE / 2;
    break;
  case DT_PROTEIC:
    numberOfMatrices = PROT_MATRIX_SIZE;
    break;
  default:
    Utilities::exit_partest(EX_OSERR);
  }
  numberOfModels *= numberOfMatrices;
  models = (Model **) malloc(numberOfModels * sizeof(Model *));

#ifdef DEBUG
  cout << "[TRACE] Creating modelset for ratevar " << rateVar
      << " (" << numberOfModels << ")" << endl;
#endif

  unsigned int current = 0;

  // Loop over the parameters
  for (bitMask rateVarLoop = 0; rateVarLoop <= rateVar>>1; rateVarLoop++) {
    // check this for avoiding duplicates
    if (!(rateVarLoop & ~(rateVar>>1))) {
#ifdef DEBUG
      cout << "[TRACE] Creating models for ratevar " << (rateVarLoop<<1)
          << " (" << current << "/" << numberOfModels << ")" << endl;
#endif
      buildModelSet(&(models[current]), rateVarLoop<<1, dataType);
      current += numberOfMatrices;
    }
  }

}

ModelSet::~ModelSet()
{
  if (models) {
    for (int i = 0; i < numberOfModels; i++) {
      delete models[i];
    }
    free(models);
  }
}

Model * ModelSet::getModel(unsigned int index) {
	  return models[index];
}

int ModelSet::buildModelSet(Model **models, bitMask rateVar,
    DataType dataType) {

  int currentIndex = 0;
  switch (dataType) {
  case DT_NUCLEIC:
	  for (int i = (rateVar & RateVarF)?1:0; i < NUC_MATRIX_SIZE; i+=2) {
	        NucMatrix nm = static_cast<NucMatrix>(i);
	        models[currentIndex++] = new NucleicModel(nm, rateVar, numberOfTaxa);
	      }
    break;
  case DT_PROTEIC:
    for (int i = 0; i < PROT_MATRIX_SIZE; i++) {
      ProtMatrix pm = static_cast<ProtMatrix>(i);
      models[i] = new ProteicModel(pm, rateVar, numberOfTaxa);
    }
    break;
  default:
    Utilities::exit_partest(EX_OSERR);
  }

  return 0;

}

} /* namespace partest */
