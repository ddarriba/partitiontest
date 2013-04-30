/*
 * NucleicModel.h
 *
 *  Created on: Jan 8, 2013
 *      Author: diego
 */

#ifndef NUCLEICMODEL_H_
#define NUCLEICMODEL_H_

#include "Model.h"
#include "../util/GlobalDefs.h"

#define NUM_NUC_FREQS 4

namespace partest {

class NucleicModel: public Model {
public:
	NucleicModel(NucMatrix matrix, bitMask rateVariation, int numberOfTaxa);
	void setFrequencies(const double * frequencies);
	void setRates(const double * rates);
	virtual ~NucleicModel();
private:
	NucMatrix matrix;
};

} /* namespace partest */
#endif /* NUCLEICMODELL_H_ */
