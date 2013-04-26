/*
 * Model.h
 *
 *  Created on: Jan 8, 2013
 *      Author: diego
 */

#ifndef MODEL_H_
#define MODEL_H_

#include "../util/GlobalDefs.h"
#include <string>

using namespace std;

namespace partest {

class Model {
public:
	/**
	 * Generic constructor of a model.
	 *
	 * @param rateVariation The rate variation (+I, +G, +F).
	 * @param numberOfTaxa Number of taxa (required for computing
	 *        the free parameters.
	 */
  Model(bitMask rateVariation, int numberOfTaxa);
  string getName(void);
  string getMatrixName(void);
  bitMask getRateVariation(void);
  string getTree(void);
  bool isPInv(void);
  bool isGamma(void);
  bool isPF(void);
  bool isOptimized(void);
  double getLnL(void) { return lnL; }
  double getAlpha(void) { return alpha; }
  double getpInv(void) { return pInv; }
  int getNumberOfFreeParameters(void) { return (modelFreeParameters + treeFreeParameters); }
  int getModelFreeParameters(void) { return modelFreeParameters; }
  int getTreeFreeParameters(void) { return treeFreeParameters; }
  void setLnL(double lnL) { this->lnL = lnL; }
  void setAlpha(double alpha);
  void setpInv(double pInv);
  void setTree(string tree);
  void setTree(char * tree);
  virtual void setFrequencies(const double * frequencies);
  virtual void setRates(const double * rates);
  virtual void print();
  virtual ~Model();
protected:
  /** Rate variation (+I, +G, +F) */
  bitMask rateVariation;
  /** Likelihood */
  double lnL;
  /** Alpha value of gamma distribution */
  double alpha;
  /** Proportion of invariable sites */
  double pInv;
  /** Substitution rates */
  double *rates;
  /** State frequencies */
  double *frequencies;
  /** Full name of the model. */
  string name;
  /** Name of the model matrix. */
  string matrixName;
  /** Number of free parameters of the model. */
  int modelFreeParameters;
  /** Number of free parameters of the tree. */
  int treeFreeParameters;
  /** Most likely tree */
  string tree;
};

} /* namespace partest */

#endif /* MODEL_H_ */
