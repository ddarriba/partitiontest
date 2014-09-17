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

#include <pllInternal.h>

/* includes for working with files/directories (Linux) */
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <unistd.h>

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

bool Utilities::isNumeric(const char * value) {
	int len = strlen(value);
	bool decimalFound = false;
	for (int i = 0; i<len; i++) {
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
	int len = strlen(value);
	for (int i = 0; i<len; i++)
		if (!isdigit(value[i])) {
			return false;
		}
	return true;
}

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
	sdX = sdX?sdX:1;
	sdY = sdY?sdY:1;
	double sum = 0.0;
	for (int i = 0; i < n; i++) {
		sum += pow((X[i] - meanX)/sdX - (Y[i] - meanY)/sdY, 2);
		//sum += pow(X[i] / meanX - Y[i] / meanY, 2);
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
		pInfo * nextPL = PartitionMap::getInstance()->getPartitionElement(nextElement)->getPartitions()->partitionData[0];
		partitions->partitionData[0]->alpha += nextPL->alpha;
		for (size_t j = 0;
				j < (data_type==DT_NUCLEIC?NUM_NUC_FREQS:NUM_PROT_FREQS);
				j++) {
			partitions->partitionData[0]->frequencies[j] +=
					nextPL->frequencies[j];
		}
		if (data_type == DT_NUCLEIC) {
			for (size_t j = 0; j < NUM_DNA_RATES; j++) {
				partitions->partitionData[0]->substRates[j] +=
						nextPL->substRates[j];
			}
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

	if (!pergene_branch_lengths || !pergene_starting_tree)
		return 0;

	pllNewickTree * nts[id.size()];
	pllStack * infoStack[id.size()];
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
	return nt;
}


int Utilities::path_is_directory (const char* path) {
    struct stat s_buf;

    if (stat(path, &s_buf))
        return 0;

    return S_ISDIR(s_buf.st_mode);
}

int Utilities::delete_folder_tree (const char* directory_name) {
	   DIR *d = opendir(directory_name);
	   size_t path_len = strlen(directory_name);
	   int r = -1;

	   if (d)
	   {
	      struct dirent *p;

	      r = 0;

	      while (!r && (p=readdir(d)))
	      {
	          int r2 = -1;
	          char *buf;
	          size_t len;

	          /* Skip the names "." and ".." as we don't want to recurse on them. */
	          if (!strcmp(p->d_name, ".") || !strcmp(p->d_name, ".."))
	          {
	             continue;
	          }

	          len = path_len + strlen(p->d_name) + 2;
	          buf = (char *) malloc(len);

	          if (buf)
	          {
	             struct stat statbuf;

	             snprintf(buf, len, "%s/%s", directory_name, p->d_name);

	             if (!stat(buf, &statbuf))
	             {
	                if (S_ISDIR(statbuf.st_mode))
	                {
	                   r2 = delete_folder_tree(buf);
	                }
	                else
	                {
	                   r2 = unlink(buf);
	                }
	             }

	             free(buf);
	          }

	          r = r2;
	      }

	      closedir(d);
	   }

	   if (!r)
	   {
	      r = rmdir(directory_name);
	   }

	   return r;
}

} /* namespace partest */
