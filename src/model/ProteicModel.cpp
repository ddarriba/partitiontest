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

#include "ProteicModel.h"
#include <stdlib.h>
#include <iostream>
#include <math.h>

namespace partest {

ProteicModel::ProteicModel(ProtMatrix matrix, bitMask rateVariation,
		int numberOfTaxa) :
		Model(rateVariation, numberOfTaxa), matrix(matrix) {
	/* treeFreeParameters is already initialized to the number of branches */
	this->numberOfFrequencies = NUM_PROT_FREQS;
	this->frequencies = (double *) malloc(NUM_PROT_FREQS * sizeof(double));

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
	for (int i = 0; i < NUM_PROT_FREQS; i++) {
		this->frequencies[i] = frequencies[i];
	}
}

void ProteicModel::setRates(const double * rates) {
	// Ignore
}

double ProteicModel::distanceTo(ProteicModel * other) {
	double matrixDistance = getEuclideanDistance(matrix, other->matrix);
	double invDistance = pInv - other->pInv;
	double shapeDistance = alpha - other->alpha;
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

double ProteicModel::getEuclideanDistance(ProtMatrix m1, ProtMatrix m2) {
	if (m1 == m2) {
		return 0;
	} else {
		ProtMatrix lowMatrix = min(m1, m2);
		ProtMatrix highMatrix = max(m1, m2);
		int index = (highMatrix * (highMatrix - 1) / 2) + lowMatrix;
		double distances[91] = {
				554.40, 745.31, 308.82, 681.72, 369.32, 253.26, 738.50, 263.33,
				207.04, 362.29, 732.96, 259.32, 207.44, 360.79, 5.99, 775.30,
				415.07, 515.03, 478.51, 525.15, 523.88, 3209.58, 2690.17,
				2541.35, 2614.11, 2548.45, 2553.82, 2553.20, 5.51, 549.44,
				740.21, 677.20, 733.42, 727.89, 771.69, 3205.23, 1.80, 554.76,
				745.81, 682.29, 738.85, 733.31, 776.16, 3210.21, 5.76, 867.34,
				579.11, 336.95, 307.70, 523.09, 523.51, 690.97, 2571.97, 862.86,
				867.91, 859.28, 474.02, 486.63, 458.27, 520.60, 521.06, 372.84,
				2406.16, 855.75, 860.15, 623.86, 32.24, 535.48, 720.23, 658.53,
				718.21, 712.70, 761.99, 3190.63, 28.48, 32.78, 838.70, 842.68,
				15.20, 549.32, 735.75, 672.36, 731.36, 725.83, 772.55, 3204.47,
				12.71, 15.65, 855.21, 856.09, 20.71 };
		return distances[index];
	}
}
} /* namespace partest */
