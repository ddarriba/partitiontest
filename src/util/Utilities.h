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
 * @file Utilities.h
 * @author Diego Darriba
 * @brief Set of utilities for PartitionTest
 */

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <pll.h>
#include <math.h>

#include "util/GlobalDefs.h"

namespace partest {

class Utilities {
public:

	/**
	 * Binary power of x.
	 * @return 2 to the x
	 */
	static unsigned long int binaryPow(unsigned long int x);

	/**
	 * Check the existence of a file
	 *
	 * @return true, if the file exists
	 */
	static inline bool existsFile (const std::string& name) {
	    if (FILE *file = fopen(name.c_str(), "r")) {
	        fclose(file);
	        return true;
	    } else {
	        return false;
	    }
	}

	static inline size_t iBinaryLog(size_t x) {
		return ceil(log(x) / log(2));
	}
	static inline double dBinaryLog(double x) {
		return (log(x) / log(2));
	}
	static inline int iDecLog(int x) {
		return x > 0 ? floor(log(x) / log(10)) : 0;
	}
	static char toBase64(int value);
	static int setbitsCount(bitMask value);

	static double mean(double series[], int n);
	static double variance(double series[], int n);
	static double standardDeviation(double series[], int n);
	static double covariance(double X[], double Y[], int n);
	static double correlation(double X[], double Y[], int n);
	static double euclideanDistance(double X[], double Y[], int n);
	static double normalizedEuclideanDistance(double X[], double Y[], int n);

	/** Number of branches according to the number of taxa. */
	static int numberOfBranches(int numTaxa);

	static void mergeIds(t_partitionElementId & dest, t_partitionElementId id1,
			t_partitionElementId id2);
	static bool intersec(t_partitionElementId & e1, t_partitionElementId & e2);
	static int duplicateAlignmentData(pllAlignmentData ** out,
			pllAlignmentData * in);

	static bool contains(t_partitionElementId vec, int num);
	static bool contains(t_partitioningScheme vec, t_partitionElementId id);

	static std::string getProtMatrixName(ProtMatrix matrix);
	static std::string getProtRaxmlName(ProtMatrix matrix);

	static int toLower(char * str);

	static void printScheme(t_partitioningScheme scheme);

	static int averageModelParameters(t_partitionElementId id, partitionList * partitions);
	static pllNewickTree * averageBranchLengths(t_partitionElementId id);
private:
	static char encoding_table[];
};

} /* namespace partest */

#endif /* UTILITIES_H_ */
