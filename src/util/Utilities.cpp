/*
 * Utilities.cpp
 *
 *  Created on: Apr 8, 2014
 *      Author: diego
 */

#include "Utilities.h"

#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <iostream>

namespace partest {

char Utilities::encoding_table[] = {
		'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N',
		'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 'a', 'b',
		'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p',
		'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '0', '1', '2', '3',
		'4', '5', '6', '7', '8', '9', '-', '_' };

char Utilities::toBase64(int value) {
	if (value > 63)
		exit_partest(EX_SOFTWARE);
	return encoding_table[value];
}

unsigned long int Utilities::binaryPow(unsigned long int x) {
	unsigned long int nextId = 1;
	nextId <<= x;
	return nextId;
}

// Brian Kernighan’s Algorithm
int Utilities::setbitsCount(bitMask n) {
	size_t count = 0;
	while (n) {
		n &= (n - 1);
		count++;
	}
	return count;
}

double Utilities::mean(double series[], int n) {
	double sum = 0.0;
	for (int i = 0; i < n; i++) {
		sum += series[i];
	}
	return sum / n;
}

double Utilities::variance(double series[], int n) {

	double vmean = mean(series, n);
	double temp = 0;

	for (int i = 0; i < n; i++)
		temp += (series[i] - vmean) * (series[i] - vmean);

	return temp / (n - 1);
}

double Utilities::standardDeviation(double series[], int n) {
	return sqrt(variance(series, n));
}

double Utilities::covariance(double X[], double Y[], int n) {

	double xmean = mean(X, n);
	double ymean = mean(Y, n);
	double total = 0;

	for (int i = 0; i < n; i++)
		total += (X[i] - xmean) * (Y[i] - ymean);

	return total / n;

}

double Utilities::euclideanDistance(double X[], double Y[], int n) {

	double sum = 0.0;
	for (int i = 0; i < n; i++) {
		sum += (X[i] - Y[i]) * (X[i] - Y[i]);
	}
	return sqrt(sum);

}
double Utilities::normalizedEuclideanDistance(double X[], double Y[], int n) {
	double meanX = mean(X, n);
	double meanY = mean(Y, n);
	//double sdX = standardDeviation(X, n);
	//double sdY = standardDeviation(Y, n);
	double sum = 0.0;
	for (int i = 0; i < n; i++) {
		//sum += pow((X[i] - meanX)/sdX - (Y[i] - meanY)/sdY, 2);
		sum += pow(X[i] / meanX - Y[i] / meanY, 2);
	}
	return sqrt(sum);
}

int Utilities::numberOfBranches(int numTaxa) {
	return 2 * numTaxa - 3;
}

void Utilities::mergeIds(t_partitionElementId & dest, t_partitionElementId id1,
		t_partitionElementId id2) {
	dest.reserve(id1.size() + id2.size());
	dest.insert(dest.end(), id1.begin(), id1.end());
	dest.insert(dest.end(), id2.begin(), id2.end());
	sort(dest.begin(), dest.end());
}

bool Utilities::intersec(t_partitionElementId & e1, t_partitionElementId & e2) {
	for (size_t i = 0; i < e1.size(); i++) {
		for (size_t j = 0; j < e2.size(); j++) {
			if (e1.at(i) == e2.at(j))
				return true;
		}
	}
	return false;
}

bool Utilities::contains(t_partitionElementId vec, int num) {
	for (size_t i=0; i < vec.size(); i++) {
		int n = vec.at(i);
		if (n == num)
			return true;
	}
	return false;
}

bool Utilities::contains(t_partitioningScheme vec, t_partitionElementId id) {
	for (size_t i=0; i<vec.size(); i++) {
		t_partitionElementId eId = vec.at(i);
		if (eId == id)
			return true;
	}
	return false;
}

int Utilities::duplicateAlignmentData(pllAlignmentData ** out,
		pllAlignmentData * in) {
	if ((*out) > 0) {
		return -1;
	}
	(*out) = pllInitAlignmentData(in->sequenceCount, in->sequenceLength);
	(*out)->siteWeights = (int *) malloc((*out)->sequenceLength * sizeof(int));
	for (int seq = 1; seq <= (*out)->sequenceCount; seq++) {
		(*out)->sequenceLabels[seq] = (char *) malloc(
				strlen(in->sequenceLabels[seq]) + 1);
		strcpy((*out)->sequenceLabels[seq], in->sequenceLabels[seq]);
		memcpy((*out)->sequenceData[seq], in->sequenceData[seq],
				(*out)->sequenceLength * sizeof(unsigned char));
	}
	for (int site = 0; site < (*out)->sequenceLength; site++) {
		(*out)->siteWeights[site] = 1;
	}
	return 0;
}

std::string Utilities::getProtRaxmlName(ProtMatrix matrix) {
	std::string matrixName = getProtMatrixName(matrix);
	std::transform(matrixName.begin(), matrixName.end(), matrixName.begin(),
			::toupper);
	return matrixName;
}

std::string Utilities::getProtMatrixName(ProtMatrix matrix) {
	std::string matrixName;
	switch (matrix) {
	case PROT_MATRIX_DAYHOFF:
		matrixName = "DAYHOFF";
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
		matrixName = "Auto";
		break;
	default:
		exit_partest(EX_SOFTWARE);
	}
	return matrixName;
}

int Utilities::toLower(char * str) {
	for (int i = 0; str[i] != '\0'; i++)
		if (str[i] >= 'A' && str[i] <= 'Z')
			str[i] += 'a' - 'A';
	return 0;
}

void Utilities::printScheme(t_partitioningScheme scheme) {
	for (size_t i=0; i<scheme.size(); i++) {
		t_partitionElementId element = scheme.at(i);
		std::cout << "(";
		for (size_t j=0; j<element.size(); j++) {
			std::cout << element.at(j);
			if (j<element.size()-1) std::cout<< " ";
		}
		std::cout << ")";
	}
	std::cout << std::endl;
}

} /* namespace partest */
