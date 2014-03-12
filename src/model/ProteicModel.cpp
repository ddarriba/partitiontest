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

double distances[153] = {
				0.055819, 6.28983, 6.28287, 8.3492, 8.34977, 6.24178, 5.55446, 5.54012, 5.00539, 7.05072,
				8.62161, 8.60373, 7.95537, 9.10421, 6.24999, 8.03345, 8.0335, 6.95272, 8.01101, 6.10837,
				9.04789, 6.57179, 6.55386, 4.20234, 6.7556, 2.81282, 5.6121, 6.42513, 8.95954, 8.94053,
				8.25574, 9.17895, 5.08662, 5.82568, 8.33647, 4.86069, 10.7929, 10.7903, 7.84703, 5.74534,
				9.63371, 10.3161, 10.5329, 8.80411, 10.4791, 7.45006, 7.43539, 5.71096, 7.2435, 3.84227,
				4.45806, 7.35689, 2.74763, 4.5702, 8.85182, 9.40526, 9.39966, 8.95556, 7.78357, 8.5285,
				7.47896, 9.37941, 7.94739, 8.64823, 8.97761, 7.13019, 8.25429, 8.25228, 7.50827, 6.61555,
				7.54048, 7.29122, 8.08111, 6.8184, 8.26604, 8.16448, 6.35032, 2.62676, 8.30528, 8.28623,
				8.06045, 8.98873, 4.89016, 5.80825, 7.80521, 4.76463, 3.43454, 10.8316, 4.64089, 7.97547,
				7.67082, 9.86055, 9.85601, 6.00977, 8.49291, 8.37583, 10.4753, 8.44097, 7.76224, 10.7074,
				9.2562, 8.6585, 11.6827, 10.6534, 10.5191, 12.6687, 12.6677, 8.81821, 10.8266, 10.9484,
				13.2696, 10.4887, 10.4906, 12.6583, 10.6809, 11.5374, 14.3116, 13.2764, 12.9748, 7.50057,
				6.15949, 6.15329, 0.500757, 6.26143, 4.83252, 7.87112, 6.78713, 4.13408, 8.13932, 7.89322,
				5.6065, 8.95895, 7.51905, 7.94798, 5.9798, 8.83705, 9.90469, 9.89408, 6.59561, 9.28638,
				7.85689, 9.44792, 8.44811, 7.32588, 9.22372, 9.81774, 8.08318, 10.56, 9.76868, 8.86959,
				7.10936, 7.86896, 6.60864
		};

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
		modelFreeParameters += (NUM_AA_RATES - 1);
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

ProtMatrix ProteicModel::getMatrix(void) {
	return matrix;
}

void ProteicModel::setFrequencies(const double * frequencies) {
	for (int i = 0; i < NUM_PROT_FREQS; i++) {
		this->frequencies[i] = frequencies[i];
	}
}

void ProteicModel::setRates(const double * rates) {
	// Ignore
}

double ProteicModel::distanceTo(Model * otherModel) {
	ProteicModel * other = static_cast<ProteicModel *>(otherModel);
	//double matrixDistance = matrix!=other->matrix?getEuclideanDistance(matrix, other->matrix): 0;
	double matrixDistance = getEuclideanDistance(matrix, other->matrix);
	double invDistance = abs(pInv - other->pInv);
	double shapeDistance = abs(alpha - other->alpha);
	double freqsDistance = 0.0;
	// TODO: Store frequencies and compute euclidean distance.
//	for (int i = 0; i < numberOfFrequencies; i++) {
//		freqsDistance += (frequencies[i] - other->frequencies[i])
//				* (frequencies[i] - other->frequencies[i]);
//	}
//	freqsDistance = sqrt(freqsDistance);

	double distance = matrixDistance + invDistance + shapeDistance
			+ freqsDistance;

	// cout << " DIST " << matrix << " to " << other->matrix << " = " << matrixDistance << " + " << invDistance << " + " << shapeDistance << " + " << freqsDistance << " = " << distance << endl;

	return distance;
}

double ProteicModel::getEuclideanDistance(ProtMatrix m1, ProtMatrix m2) {
	if (m1 == m2) {
		return 0;
	} else {
		ProtMatrix lowMatrix = min(m1, m2);
		ProtMatrix highMatrix = max(m1, m2);
		int index = (highMatrix * (highMatrix - 1) / 2) + lowMatrix;

		return distances[index];
	}
}

void ProteicModel::print(ostream& cout, const char * prefix) {
	cout << prefix << "Name:  " << name << endl;
	if (isOptimized()) {
		cout << prefix << "lnL:   " << lnL << endl;
		if (isPInv()) {
			cout << prefix << "pInv:  " << pInv << endl;
		}
		if (isGamma()) {
			cout << prefix << "alpha: " << alpha << endl;
		}
		cout << prefix << "Most Likely Tree: " << tree;
	} else {
		cout << prefix << "State: Unoptimized" << endl;
	}
}

} /* namespace partest */
