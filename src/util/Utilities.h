/*
 * Utilities.h
 *
 *  Created on: Apr 8, 2014
 *      Author: diego
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
	static inline unsigned int iBinaryLog(unsigned int x) {
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

private:
	static char encoding_table[];
};

} /* namespace partest */

#endif /* UTILITIES_H_ */
