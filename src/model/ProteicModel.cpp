/*
 * ProteicModel.cpp
 *
 *  Created on: Jan 8, 2013
 *      Author: diego
 */

#include "ProteicModel.h"
#include <stdlib.h>
#include <iostream>

namespace partest {

ProteicModel::ProteicModel(ProtMatrix matrix, bitMask rateVariation,
		int numberOfTaxa) :
		Model(rateVariation, numberOfTaxa), matrix(matrix) {
	/* treeFreeParameters is already initialized to the number of branches */
	this->frequencies = (double *) malloc(20 * sizeof(double));

	switch (matrix) {
	case PROT_MATRIX_DAYHOFF:
		matrixName = "Dayhoff";
		break;
	case PROT_MATRIX_DCMUT:
		matrixName = "DCMut";
		break;
	case PROT_MATRIX_JTT:
		matrixName = "JTT";
		break;
	case PROT_MATRIX_MTREV:
		matrixName = "MtREV";
		break;
	case PROT_MATRIX_WAG:
		matrixName = "WAG";
		break;
	case PROT_MATRIX_RTREV:
		matrixName = "RtREV";
		break;
	case PROT_MATRIX_CPREV:
		matrixName = "CpREV";
		break;
	case PROT_MATRIX_VT:
		matrixName = "VT";
		break;
	case PROT_MATRIX_BLOSUM62:
		matrixName = "Blosum62";
		break;
	case PROT_MATRIX_MTMAM:
		matrixName = "MtMam";
		break;
	case PROT_MATRIX_LG:
		matrixName = "LG";
		break;
	case PROT_MATRIX_MTART:
		matrixName = "MtArt";
		break;
	case PROT_MATRIX_HIVB:
		matrixName = "HIVb";
		break;
	case PROT_MATRIX_HIVW:
		matrixName = "HIVw";
		break;
#ifdef _PLL
		case PROT_MATRIX_MTZOA:
		matrixName = "MtZoa";
		break;
		case PROT_MATRIX_PMB:
		matrixName = "PMB";
		break;
		case PROT_MATRIX_JTTDCMUT:
		matrixName = "JTTDCMut";
		break;
		case PROT_MATRIX_FLU:
		matrixName = "FLU";
		break;
		case PROT_MATRIX_AUTO:
		matrixName = "Auto";
		break;
		case PROT_MATRIX_GTR:
		matrixName = "GTR";
		/* 189 substitution rates have to be optimized! (190-1) */
		modelFreeParameters += 189;
		break;
#endif
	}

	name = matrixName;
	if (rateVariation & RateVarI) {
		name += "+I";
		/* proportion of invariable sites free parameter */
		modelFreeParameters++;
	}
	if (rateVariation & RateVarG) {
		name += "+G";
		/* alpha free parameter */
		modelFreeParameters++;
	}
	if (rateVariation & RateVarF) {
		name += "+F";
		/* 19 frequencies free parameters (20-1) */
		modelFreeParameters += 19;
	}

}

ProteicModel::~ProteicModel() {
	// NOTHING
}

void ProteicModel::setFrequencies(const double * frequencies) {
	for (int i = 0; i < 20; i++) {
		this->frequencies[i] = frequencies[i];
	}
}

void ProteicModel::setRates(const double * rates) {
	// Ignore
}

} /* namespace partest */
