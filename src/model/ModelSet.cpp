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
#include <string.h>

namespace partest {

ModelSet::ModelSet(bitMask rateVar, DataType dataType, int numberOfTaxa,
		OptimizeMode optimizeMode, bool forceCompleteSet) :
		rateVar(rateVar), dataType(dataType), numberOfTaxa(numberOfTaxa), optimizeMode(
				optimizeMode) {
	int numberOfParameters = Utilities::setbitsCount(rateVar >> 1);
	int numberOfMatrices;
	numberOfModels = Utilities::binaryPow(numberOfParameters);

	if (forceCompleteSet)
		optimizeMode = OPT_SEARCH;

	if (optimizeMode == OPT_SEARCH) {
		switch (dataType) {
		case DT_NUCLEIC:
			numberOfMatrices = NUC_MATRIX_SIZE / 2;
			break;
		case DT_PROTEIC:
			numberOfMatrices = PROT_MATRIX_SIZE;
			break;
		default:
			Utilities::exit_partest(EX_OSERR);
			break;
		}
	} else {
		numberOfMatrices = 1;
	}

	if (optimizeMode == OPT_SEARCH) {
		numberOfModels *= numberOfMatrices;
	} else {
		numberOfModels = 1;
	}

	models = (Model **) malloc(numberOfModels * sizeof(Model *));

#ifdef DEBUG
	cout << "[TRACE] Creating modelset for ratevar " << rateVar
	<< " (" << numberOfModels << ")" << endl;
#endif

	unsigned int current = 0;

	if (optimizeMode == OPT_GTR) {
		buildModelSet(&(models[0]), 0, false);
	} else {
		// Loop over the parameters
		for (bitMask rateVarLoop = 0; rateVarLoop <= rateVar >> 1;
				rateVarLoop++) {
			// check this for avoiding duplicates
			if (!(rateVarLoop & ~(rateVar >> 1))) {
#ifdef DEBUG
				cout << "[TRACE] Creating models for ratevar " << (rateVarLoop<<1)
				<< " (" << current << "/" << numberOfModels << ")" << endl;
#endif
				buildModelSet(&(models[current]), rateVarLoop << 1,
						forceCompleteSet);
				current += numberOfMatrices;
			}
		}
	}
}

ModelSet::~ModelSet() {
	if (models) {
		for (unsigned int i = 0; i < numberOfModels; i++) {
			delete models[i];
		}
		free(models);
	}
}

void ModelSet::buildCompleteModelSet(bool clearAll) {
	int numberOfParameters = Utilities::setbitsCount(rateVar >> 1);
	int numberOfMatrices = 0;
	int newNumberOfModels = Utilities::binaryPow(numberOfParameters);

	switch (dataType) {
	case DT_NUCLEIC:
		numberOfMatrices = NUC_MATRIX_SIZE / 2;
		break;
	case DT_PROTEIC:
		numberOfMatrices = PROT_MATRIX_SIZE;
		break;
	default:
		Utilities::exit_partest(EX_OSERR);
		break;
	}

	newNumberOfModels *= numberOfMatrices;
	Model ** newModels = (Model **) malloc(newNumberOfModels * sizeof(Model *));

	unsigned int current = 0;
	// Loop over the parameters
	for (bitMask rateVarLoop = 0; rateVarLoop <= rateVar >> 1; rateVarLoop++) {
		// check this for avoiding duplicates
		if (!(rateVarLoop & ~(rateVar >> 1))) {
#ifdef DEBUG
			cout << "[TRACE] Creating models for ratevar " << (rateVarLoop<<1)
			<< " (" << current << "/" << numberOfModels << ")" << endl;
#endif
			buildModelSet(&(newModels[current]), rateVarLoop << 1, true);
			current += numberOfMatrices;
		}
	}

	if (clearAll) {
		for (unsigned int i = 0; i < numberOfModels; i++) {
			delete models[i];
		}
	} else {
		for (unsigned int i = 0; i < numberOfModels; i++) {
			for (int j = 0; j < newNumberOfModels; j++) {
				if (!strcmp(newModels[j]->getName().c_str(),
						models[i]->getName().c_str())) {
					delete newModels[j];
					newModels[j] = models[i];
					break;
				}
			}
		}
	}
	free(models);
	models = newModels;
	numberOfModels = newNumberOfModels;
}

Model * ModelSet::getModel(unsigned int index) {
	return models[index];
}

int ModelSet::buildModelSet(Model **models, bitMask rateVar,
		bool forceCompleteSet) {

#ifdef _PLL
	rateVar |= RateVarG;
#endif

	int currentIndex = 0;
	NucMatrix nm;
	switch (dataType) {
	case DT_NUCLEIC:
		if (optimizeMode == OPT_GTR) {
			nm = NUC_MATRIX_GTR;
#ifdef _PLL
			models[currentIndex++] = new NucleicModel(nm, RateVarG | RateVarF, numberOfTaxa);
#else
			models[currentIndex++] = new NucleicModel(nm,
					RateVarG | RateVarF | RateVarI, numberOfTaxa);
#endif
		} else {
			for (int i = (rateVar & RateVarF) ? 1 : 0; i < NUC_MATRIX_SIZE; i +=
					2) {
				nm = static_cast<NucMatrix>(i);
				models[currentIndex++] = new NucleicModel(nm, rateVar,
						numberOfTaxa);
			}
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
		break;
	}

	return 0;

}

} /* namespace partest */
