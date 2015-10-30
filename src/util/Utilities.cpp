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
 *  For any other inquiries send an Email to Diego Darriba
 *  ddarriba@udc.es
 */

/**
 * @file Utilities.cpp
 * @author Diego Darriba
 */

#include "Utilities.h"
#include "indata/PartitionMap.h"
#include "indata/TreeManager.h"

#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <iostream>
#include <cassert>

namespace partest {

char Utilities::encoding_table[] = {
		'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N',
		'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 'a', 'b',
		'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p',
		'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '0', '1', '2', '3',
		'4', '5', '6', '7', '8', '9', '-', '_' };

bool Utilities::isNumeric(const char * value) {
	size_t len = strlen(value);
	bool decimalFound = false;
	for (size_t i = 0; i<len; i++) {
		if (!isdigit(value[i])) {
			if (value[i] == '.' && !decimalFound) {
				decimalFound = true;
			} else {
				return false;
			}
		}
	}
	return true;
}

bool Utilities::isInteger(const char * value) {
	size_t len = strlen(value);
	for (size_t i = 0; i<len; i++)
		if (!isdigit(value[i])) {
			return false;
		}
	return true;
}

char Utilities::toBase64(int value) {
	assert (value < 64);
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
	return (int) count;
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

double Utilities::euclideanDistance(double X[], double Y[], int n, double multiplier) {

	double sum = 0.0;
	for (int i = 0; i < n; i++) {
		sum += (X[i]*multiplier - Y[i]) * (X[i]*multiplier - Y[i]);
	}
	return sqrt(sum);

}
double Utilities::normalizedEuclideanDistance(double X[], double Y[], int n) {
	double meanX = mean(X, n);
	double meanY = mean(Y, n);
	double sdX = standardDeviation(X, n);
	double sdY = standardDeviation(Y, n);
	sdX = isfinite(sdX)?sdX:1;
	sdY = isfinite(sdY)?sdY:1;
	double sum = 0.0;
	for (int i = 0; i < n; i++) {
		sum += pow((X[i] - meanX)/sdX - (Y[i] - meanY)/sdY, 2);
		//sum += pow(X[i] / meanX - Y[i] / meanY, 2);
	}
	return sqrt(sum);
}

int Utilities::numSchemesHierarchicalClustering(int n) {
	int nSchemes = n + 1;
	int i;
	if (max_samples) {
		for (i=2; i<n; i++)
			nSchemes += max_samples<(i*(i-1)/2)?max_samples:(i*(i-1)/2);
	} else {
		int n_iter_samples;
		for (i=2; i<n; i++)
		{
			n_iter_samples = (int)(samples_percent * (i*(i-1)/2));
			nSchemes += n_iter_samples?n_iter_samples:1;
		}
	}
	return nSchemes;
}

int Utilities::numSchemesGreedy(int n) {
	return n + (n-1)*(n + n-2)/2;
}

int Utilities::numSchemesAutoSearch(int n) {
	return numSchemesHierarchicalClustering(n);
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
		int n = (int) vec.at(i);
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
	(*out)->siteWeights = (int *) malloc((size_t)(*out)->sequenceLength * sizeof(int));
	for (int seq = 1; seq <= (*out)->sequenceCount; seq++) {
		(*out)->sequenceLabels[seq] = (char *) malloc(
				strlen(in->sequenceLabels[seq]) + 1);
		strcpy((*out)->sequenceLabels[seq], in->sequenceLabels[seq]);
		memcpy((*out)->sequenceData[seq], in->sequenceData[seq],
				(size_t)(*out)->sequenceLength * sizeof(unsigned char));
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
		assert(0);
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

int Utilities::averageModelParameters(t_partitionElementId id, partitionList * partitions) {
	if (partitions->numberOfPartitions != 1)
		return -1;
	partitions->partitionData[0]->alpha = 0.0;
	for (size_t i=0; i<(data_type==DT_NUCLEIC?NUM_NUC_FREQS:NUM_PROT_FREQS); i++) {
		partitions->partitionData[0]->frequencies[i] = 0.0;
	}
	if (data_type == DT_NUCLEIC) {
		for (size_t i=0; i<NUM_DNA_RATES; i++) {
			partitions->partitionData[0]->substRates[i] = 0.0;
		}
	}
	for (size_t i=0; i<id.size(); i++) {
		t_partitionElementId nextElement(1);
		nextElement.at(0) = id.at(i);
		TreeManager * elementTree = PartitionMap::getInstance()->getPartitionElement(nextElement)->getTreeManager();
		partitions->partitionData[0]->alpha += elementTree->getAlpha();
		memcpy(partitions->partitionData[0]->frequencies, elementTree->getFrequencies(),
				(data_type==DT_NUCLEIC?NUM_NUC_FREQS:NUM_PROT_FREQS)*sizeof(double));
		if (data_type == DT_NUCLEIC) {
			memcpy(partitions->partitionData[0]->substRates, elementTree->getRates(), NUM_DNA_RATES*sizeof(double));
		}
	}
	partitions->partitionData[0]->alpha /= id.size();
	for (size_t i=0; i<(data_type==DT_NUCLEIC?NUM_NUC_FREQS:NUM_PROT_FREQS); i++) {
		partitions->partitionData[0]->frequencies[i] /= id.size();
		partitions->partitionData[0]->freqExponents[i] = 1.0;
	}
	if (data_type == DT_NUCLEIC) {
		for (size_t i=0; i<NUM_DNA_RATES; i++) {
			partitions->partitionData[0]->substRates[i] /= id.size();
		}
	}
	return 0;
}

pllNewickTree * Utilities::averageBranchLengths(t_partitionElementId id) {

	assert(pergene_branch_lengths && pergene_starting_tree);

	pllNewickTree ** nts = (pllNewickTree **) malloc(id.size() * sizeof(pllNewickTree *));
	pllStack ** infoStack = (pllStack **) malloc (id.size() * sizeof(pllStack *));

	if (!(nts && infoStack)) {
		std::cerr << "[ERROR] Cannot allocate memory." << std::endl;
		exit_partest(EX_IOERR);
	}

	pllNewickNodeInfo * ninfo, *updateInfo;
	for (size_t i=0; i<id.size(); i++) {
		nts[i] = pllNewickParseString(pergene_starting_tree[id.at(0)]);
		infoStack[i] = nts[i]->tree;
	}
	pllNewickTree * nt = nts[0];

	while(infoStack[0]) {
		/* sequentially update branch lengths on target structure */
		double sumBranches = 0.0;
		updateInfo = (pllNewickNodeInfo *) infoStack[0]->item;
		for (size_t i=0; i<id.size(); i++) {
			ninfo = (pllNewickNodeInfo *) infoStack[i]->item;
			sumBranches += atof(ninfo->branch);
			infoStack[i] = infoStack[i]->next;
		}
		sprintf(updateInfo->branch, "%lf", sumBranches/id.size());
	}

	ninfo = (pllNewickNodeInfo *) nt->tree->item;
	for (size_t i=1; i<id.size(); i++) {
		/* remove newick structures but the one to return */
		pllNewickParseDestroy(&nts[i]);
	}

	free ( nts );
	free ( infoStack );

	return nt;
}

void Utilities::smoothFrequencies(double *frequencies, int numberOfFrequencies) {
	int countScale = 0, l, loopCounter = 0;

	for (l = 0; l < numberOfFrequencies; l++)
		if (frequencies[l] < FREQ_MIN)
			countScale++;

	if (countScale > 0) {
		while (countScale > 0) {
			double correction = 0.0, factor = 1.0;

			for (l = 0; l < numberOfFrequencies; l++) {
				if (frequencies[l] == 0.0)
					correction += FREQ_MIN;
				else if (frequencies[l] < FREQ_MIN) {
					correction += (FREQ_MIN - frequencies[l]);
					factor -= (FREQ_MIN - frequencies[l]);
				}
			}

			countScale = 0;

			for (l = 0; l < numberOfFrequencies; l++) {
				if (frequencies[l] >= FREQ_MIN)
					frequencies[l] = frequencies[l] - (frequencies[l] * correction * factor);
				else
					frequencies[l] = FREQ_MIN;

				if (frequencies[l] < FREQ_MIN)
					countScale++;
			}
			assert(loopCounter < 100);
			loopCounter++;
		}
	}
}

/******************************************************************************/
/* BRENT'S OPTIMIZATION */
/******************************************************************************/

#define ITMAX 100
#define CGOLD 0.3819660
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

double Utilities::brent_opt(double ax, double bx, double cx, double tol,
		double *foptx, double *f2optx, double fax, double fbx, double fcx,
		void * params, double (*target_funk)(void *, double)) {
	int iter;
	double a, b, d = 0, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x,
			xm;
	double xw, wv, vx;
	double e = 0.0;

	a = (ax < cx ? ax : cx);
	b = (ax > cx ? ax : cx);
	x = bx;
	fx = fbx;
	if (fax < fcx) {
		w = ax;
		fw = fax;
		v = cx;
		fv = fcx;
	} else {
		w = cx;
		fw = fcx;
		v = ax;
		fv = fax;
	}

	for (iter = 1; iter <= ITMAX; iter++) {
		xm = 0.5 * (a + b);
		tol2 = 2.0 * (tol1 = tol * fabs(x) + ZEPS);
		if (fabs(x - xm) <= (tol2 - 0.5 * (b - a))) {
			*foptx = fx;
			xw = x - w;
			wv = w - v;
			vx = v - x;
			*f2optx = 2.0 * (fv * xw + fx * wv + fw * vx)
					/ (v * v * xw + x * x * wv + w * w * vx);
			return x;
		}

		if (fabs(e) > tol1) {
			r = (x - w) * (fx - fv);
			q = (x - v) * (fx - fw);
			p = (x - v) * q - (x - w) * r;
			q = 2.0 * (q - r);
			if (q > 0.0)
				p = -p;
			q = fabs(q);
			etemp = e;
			e = d;
			if (fabs(p) >= fabs(0.5 * q * etemp) || p <= q * (a - x)
					|| p >= q * (b - x))
				d = CGOLD * (e = (x >= xm ? a - x : b - x));
			else {
				d = p / q;
				u = x + d;
				if (u - a < tol2 || b - u < tol2)
					d = SIGN(tol1, xm - x);
			}
		} else {
			d = CGOLD * (e = (x >= xm ? a - x : b - x));
		}

		u = (fabs(d) >= tol1 ? x + d : x + SIGN(tol1, d));
		fu = target_funk(params, u);
		if (fu <= fx) {
			if (u >= x)
				a = x;
			else
				b = x;

			SHFT(v, w, x, u)
			SHFT(fv, fw, fx, fu)
		} else {
			if (u < x)
				a = u;
			else
				b = u;
			if (fu <= fw || w == x) {
				v = w;
				w = u;
				fv = fw;
				fw = fu;
			} else if (fu <= fv || v == x || v == w) {
				v = u;
				fv = fu;
			}
		}
	}

	*foptx = fx;
	xw = x - w;
	wv = w - v;
	vx = v - x;
	*f2optx = 2.0 * (fv * xw + fx * wv + fw * vx)
			/ (v * v * xw + x * x * wv + w * w * vx);

	return x;
}

/* most of the code for Brent optimization taken from IQ-Tree
 * http://www.cibiv.at/software/iqtree
 * --------------------------------------------------------------------------
 * Bui Quang Minh, Minh Anh Thi Nguyen, and Arndt von Haeseler (2013)
 * Ultrafast approximation for phylogenetic bootstrap.
 * Mol. Biol. Evol., 30:1188-1195. (free reprint, DOI: 10.1093/molbev/mst024) */
double Utilities::minimize_brent(double xmin, double xguess, double xmax, double xtol,
		double *fx, double *f2x, void * params,
		double (*target_funk)(void *, double))
{
  double eps, optx, ax, bx, cx, fa, fb, fc;
  int outbounds_ax, outbounds_cx;

  /* first attempt to bracketize minimum */
  if (xguess < xmin)
    xguess = xmin;
  if (xguess > xmax)
    xguess = xmax;
  eps = xguess * xtol * 50.0;
  ax = xguess - eps;
  outbounds_ax = ax < xmin;
  if (outbounds_ax)
    ax = xmin;
  bx = xguess;
  cx = xguess + eps;
  outbounds_cx = cx > xmax;
  if (outbounds_cx)
    cx = xmax;

  /* check if this works */
  fa = target_funk (params, ax);
  fb = target_funk (params, bx);
  fc = target_funk (params, cx);

  /* if it works use these borders else be conservative */
  if ((fa < fb) || (fc < fb))
  {
    if (!outbounds_ax)
      fa = target_funk (params, xmin);
    if (!outbounds_cx)
      fc = target_funk (params, xmax);
    ax = xmin;
    cx = xmax;
  }

  optx = brent_opt (ax, bx, cx, xtol, fx, f2x, fa, fb, fc, params,
                    target_funk);
  if (*fx > fb) // if worse, return initial value
  {
    *fx = target_funk (params, bx);
    return bx;
  }

  return optx; /* return optimal x */
}

} /* namespace partest */
