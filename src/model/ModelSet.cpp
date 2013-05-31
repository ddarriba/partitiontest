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

#include "ModelSet.h"

namespace partest {

ModelSet::ModelSet(bitMask rateVar, DataType dataType, int numberOfTaxa) :
		rateVar(rateVar), dataType(dataType), numberOfTaxa(numberOfTaxa) {
	int numberOfParameters = Utilities::setbitsCount(rateVar >> 1);
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
	numberOfMatrices = 1;
	numberOfModels *= numberOfMatrices;
	models = (Model **) malloc(numberOfModels * sizeof(Model *));

#ifdef DEBUG
	cout << "[TRACE] Creating modelset for ratevar " << rateVar
	<< " (" << numberOfModels << ")" << endl;
#endif

	unsigned int current = 0;

	// Loop over the parameters
	for (bitMask rateVarLoop = 0; rateVarLoop <= rateVar >> 1; rateVarLoop++) {
		// check this for avoiding duplicates
		if (!(rateVarLoop & ~(rateVar >> 1))) {
#ifdef DEBUG
			cout << "[TRACE] Creating models for ratevar " << (rateVarLoop<<1)
			<< " (" << current << "/" << numberOfModels << ")" << endl;
#endif
			buildModelSet(&(models[current]), rateVarLoop << 1);
			current += numberOfMatrices;
		}
	}

}

ModelSet::~ModelSet() {
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

int ModelSet::buildModelSet(Model **models, bitMask rateVar) {

	int currentIndex = 0;
	NucMatrix nm;
	switch (dataType) {
	case DT_NUCLEIC:
//		for (int i = (rateVar & RateVarF) ? 1 : 0; i < NUC_MATRIX_SIZE; i +=
//				2) {
//			NucMatrix nm = static_cast<NucMatrix>(i);
//			models[currentIndex++] = new NucleicModel(nm, rateVar,
//					numberOfTaxa);
//		}
		nm = (rateVar & RateVarF) ? NUC_MATRIX_GTR : NUC_MATRIX_SYM;
		models[currentIndex++] = new NucleicModel(nm, rateVar,
							numberOfTaxa);
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
