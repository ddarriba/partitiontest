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

#include "NucleicModel.h"
#include <stdlib.h>
#include <iostream>
#include <math.h>

namespace partest {

NucleicModel::NucleicModel(NucMatrix matrix, bitMask rateVariation,
		int numberOfTaxa) :
		Model(rateVariation, numberOfTaxa), matrix(matrix) {
	/* treeFreeParameters is already initialized to the number of branches */
	this->numberOfFrequencies = NUM_NUC_FREQS;
	this->frequencies = (double *) malloc(NUM_NUC_FREQS * sizeof(double));
	this->rates = (double *) malloc(NUM_RATES * sizeof(double));
	switch (matrix) {
	case NUC_MATRIX_JC:
		name.append("JC");
		matrixName.append("000000");
		break;
	case NUC_MATRIX_F81:
		name = "F81";
		matrixName.append("000000");
		break;
	case NUC_MATRIX_K80:
		name.append("K80");
		matrixName.append("010010");
		modelFreeParameters+=1;
		break;
	case NUC_MATRIX_HKY:
		name.append("HKY");
		matrixName.append("010010");
		modelFreeParameters+=1;
		break;
	case NUC_MATRIX_TrNef:
		name.append("TrNef");
		matrixName.append("010020");
		modelFreeParameters+=2;
		break;
	case NUC_MATRIX_TrN:
		name.append("TrN");
		matrixName.append("010020");
		modelFreeParameters+=2;
		break;
	case NUC_MATRIX_TPM1:
		name.append("TPM1");
		matrixName.append("012210");
		modelFreeParameters+=2;
		break;
	case NUC_MATRIX_TPM1uf:
		name.append("TPM1uf");
		matrixName.append("012210");
		modelFreeParameters+=2;
		break;
	case NUC_MATRIX_TPM2:
		name.append("TPM2");
		matrixName.append("010212");
		modelFreeParameters+=2;
		break;
	case NUC_MATRIX_TPM2uf:
		name.append("TPM2uf");
		matrixName.append("010212");
		modelFreeParameters+=2;
		break;
	case NUC_MATRIX_TPM3:
		name.append("TPM3");
		matrixName.append("012012");
		modelFreeParameters+=2;
		break;
	case NUC_MATRIX_TPM3uf:
		name.append("TPM3uf");
		matrixName.append("012012");
		modelFreeParameters+=2;
		break;
	case NUC_MATRIX_TIM1ef:
		name.append("TIM1ef");
		matrixName.append("012230");
		modelFreeParameters+=3;
		break;
	case NUC_MATRIX_TIM1:
		name.append("TIM1");
		matrixName.append("012230");
		modelFreeParameters+=3;
		break;
	case NUC_MATRIX_TIM2ef:
		name.append("TIM2ef");
		matrixName.append("010232");
		modelFreeParameters+=3;
		break;
	case NUC_MATRIX_TIM2:
		name.append("TIM2");
		matrixName.append("010232");
		modelFreeParameters+=3;
		break;
	case NUC_MATRIX_TIM3ef:
		name.append("TIM3ef");
		matrixName.append("012032");
		modelFreeParameters+=3;
		break;
	case NUC_MATRIX_TIM3:
		name.append("TIM3");
		matrixName.append("012032");
		modelFreeParameters+=3;
		break;
	case NUC_MATRIX_TVMef:
		name.append("TVMef");
		matrixName.append("012314");
		modelFreeParameters+=4;
		break;
	case NUC_MATRIX_TVM:
		name.append("TVM");
		matrixName.append("012314");
		modelFreeParameters+=4;
		break;
	case NUC_MATRIX_SYM:
		name.append("SYM");
		matrixName.append("012345");
		modelFreeParameters+=5;
		break;
	case NUC_MATRIX_GTR:
		name.append("GTR");
		matrixName.append("012345");
		modelFreeParameters+=5;
		break;
	default:
		cerr << "ERROR: Unknown nucleic matrix" << endl;
		exit(-1);
#ifdef _PLL

#endif
	}

	if (rateVariation & RateVarF) {
		modelFreeParameters+=3;
	}

	if (rateVariation & RateVarI) {
		name .append("+I");
		/* proportion of invariable sites free parameter */
		modelFreeParameters++;
	}

	if (rateVariation & RateVarG) {
		name.append("+G");
		/* alpha free parameter */
		modelFreeParameters++;
	}
}

NucleicModel::~NucleicModel() {
	// NOTHING
}

void NucleicModel::setFrequencies(const double * frequencies) {
	for (int i = 0; i < 4; i++) {
		this->frequencies[i] = frequencies[i];
	}
}

void NucleicModel::setRates(const double * rates) {
	for (int i=0;i<NUM_RATES; i++)
		this->rates[i] = rates[i];
}

double NucleicModel::distanceTo(Model * otherModel) {
	NucleicModel * other = static_cast<NucleicModel *>(otherModel);
	double invDistance = pInv - other->pInv;
	double shapeDistance = alpha - other->alpha;
	double matrixDistance = 0.0;
	for (int i = 0; i < 6; i++) {
			matrixDistance += (rates[i] - other->rates[i])
					* (rates[i] - other->rates[i]);
		}
	double freqsDistance = 0.0;
	for (int i = 0; i < numberOfFrequencies; i++) {
		freqsDistance += (frequencies[i] - other->frequencies[i])
				* (frequencies[i] - other->frequencies[i]);
	}
	freqsDistance = sqrt(freqsDistance);

	double distance = matrixDistance + invDistance + shapeDistance
			+ freqsDistance;

	return distance;
}

} /* namespace partest */
