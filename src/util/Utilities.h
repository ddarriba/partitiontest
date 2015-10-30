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

#include <math.h>

#include "util/GlobalDefs.h"

#define FREQ_MIN 0.001

namespace partest
{

  class Utilities
  {
  public:

    /* ********* NUMERIC UTILITIES ********* */

    /**
     * Binary power of x.
     * @return 2 to the x
     */
    static unsigned long int binaryPow (unsigned long int x);

    /**
     * @brief Compute the unsigned integer binary logarithm
     */
    static inline size_t iBinaryLog (size_t x)
    {
      return (size_t) ceil (log (x) / log (2));
    }

    /**
     * @brief Compute the double precission binary logarithm
     */
    static inline double dBinaryLog (double x)
    {
      return (log (x) / log (2));
    }

    /**
     * @brief Compute the integer decimal logarithm
     */
    static inline int iDecLog (int x)
    {
      return (int) (x > 0 ? floor (log (x) / log (10)) : 0);
    }

    /**
     * @brief Check if a string value is numeric
     */
    static bool isNumeric (const char * value);

    /**
     * @brief Check if a string value is integer
     */
    static bool isInteger (const char * value);

    /**
     * @brief Convert an integer value to base 64 char
     * \pre { value < 64 }
     */
    static char toBase64 (int value);

    /**
     * @brief Count the number of 1s in a binary string
     */
    static int setbitsCount (bitMask value);

    /* ********* STATISTICS ********* */

    /**
     * @brief Compute the mean value of a doubles vector
     */
    static double mean (double series[], int n);

    /**
     * @brief Compute the variance of a doubles vector
     */
    static double variance (double series[], int n);

    /**
     * @brief Compute the standard deviation of a doubles vector
     */
    static double standardDeviation (double series[], int n);

    /**
     * @brief Compute the covariance of 2 doubles vectors
     */
    static double covariance (double X[], double Y[], int n);

    /**
     * @brief Compute the correlation of 2 doubles vectors
     */
    static double correlation (double X[], double Y[], int n);

    /**
     * @brief Compute the euclidean distance between 2 doubles vectors
     */
    static double euclideanDistance (double X[], double Y[], int n,
                                     double multiplier = 1.0);

    /**
     * @brief Compute the normalized euclidean distance between 2 doubles vectors
     */
    static double normalizedEuclideanDistance (double X[], double Y[], int n);

    static int numSchemesHierarchicalClustering (int numDataBlocks);

    static int numSchemesGreedy (int numDataBlocks);

    static int numSchemesAutoSearch (int numDataBlocks);

    /* ********* IDENTIFIERS MANAGEMENT ********* */

    /**
     * @brief Merge 2 partition element ids into a new one
     */
    static void mergeIds (t_partitionElementId & dest, t_partitionElementId id1,
                          t_partitionElementId id2);
    /**
     * @brief Get the intersection between 2 partition element ids
     */
    static bool intersec (t_partitionElementId & e1, t_partitionElementId & e2);

    /**
     * @brief Check if a partition element id contains a certain partition
     */
    static bool contains (t_partitionElementId vec, int num);

    /**
     * @brief Check if a partitioning scheme contains a certain partition element
     */
    static bool contains (t_partitioningScheme vec, t_partitionElementId id);

    /**
     * @brief Clone PLL alignment data
     */
    static int duplicateAlignmentData (pllAlignmentData ** out,
                                       pllAlignmentData * in);

    /**
     * @brief Print a partitioning scheme identifier
     */
    static void printScheme (t_partitioningScheme scheme);

    /* ********* MATRICES NAMING ********* */

    /**
     * @brief Get the name of the amino acid replacement matrix
     */
    static std::string getProtMatrixName (ProtMatrix matrix);
    /**
     * @brief Get the name of the amino acid replacement matrix as used in RAxML
     */
    static std::string getProtRaxmlName (ProtMatrix matrix);

    /* ********* MISCELANEOUS ********* */

    /**
     * @brief Convert a string to lower case
     */
    static int toLower (char * str);

    /**
     * @brief Compute number of branches according to the number of taxa.
     */
    static int numberOfBranches (int numTaxa);

    /**
     * @brief Compute the average of the model parameters
     */
    static int averageModelParameters (t_partitionElementId id,
                                       partitionList * partitions);

    /**
     * @brief Get a newick tree with the average branch lengths
     */
    static pllNewickTree * averageBranchLengths (t_partitionElementId id);

    /**
     * @brief Smooth the base frequencies such that there is no frequency below FREQ_MIN
     */
    static void smoothFrequencies (double *frequencies,
                                   int numberOfFrequencies);

    static double minimize_brent (double xmin, double xguess, double xmax,
                                  double xtol, double *fx, double *f2x,
                                  void * params,
                                  double (*target_funk) (void *, double));
  private:
    static char encoding_table[]; /** Table for encoding base 64 values */

    static double brent_opt (double ax, double bx, double cx, double tol,
                             double *foptx, double *f2optx, double fax,
                             double fbx, double fcx, void * params,
                             double (*target_funk) (void *, double));
  };

} /* namespace partest */

#endif /* UTILITIES_H_ */
