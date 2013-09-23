#pragma once
#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include "GlobalDefs.h"
#include <map>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

using namespace std;

namespace partest {

typedef map<const string, string> CfgMap;

class Utilities {
private:
	static CfgMap config;
public:

	static void vprint(ostream & out, t_partitionElementId & v);

	static bool intersec(t_partitionElementId & e1, t_partitionElementId & e2);

	/**
	 * Common function for exiting genomictest.
	 */
	static void exit_partest(int exit_value);
	/**
	 * Computes the factorial of a number
	 *
	 */
	static long factorial(unsigned int x);
	static long combinatorial(unsigned int a, unsigned int b);
	static int bell(int n);
	static double mean(double series[], int numberOfElements);
	static double variance(double series[], int numberOfElements);
	static double standardDeviation(double series[], int numberOfElements);
	static double covariance(double XSeries[], double YSeries[],
			int numberOfElements);
	static double correlation(double XSeries[], double YSeries[],
			int numberOfElements);
	static double euclideanDistance(double XSeries[], double YSeries[],
			int numberOfElements);
	static double normalizedEuclideanDistance(double XSeries[],
			double YSeries[], int numberOfElements);

	static bool existProperty(const string& property);
	static string getValue(const string& property);

	static double stringToDouble(const string& s);
	static int stringToInt(const string& s);
	static string timeToString(time_t time);
	static int numDigits(int number);
	static int setbitsCount(unsigned int value);
	static inline string toString(int n) {
		stringstream ss;
		ss << n;
		return ss.str();
	}
	static inline bool isPowerOfTwo(unsigned int x) {
		return ((x != 0) && !(x & (x - 1)));
	}
	static inline unsigned int binaryPow(unsigned int x) {
		unsigned int nextId = 1;
		nextId <<= x;
		return nextId;
	}
	static inline unsigned int binaryLog(unsigned int x) {
		return ceil(log(x) / log(2));
	}
	static inline double dBinaryLog(double x) {
		return (log(x) / log(2));
	}
	static inline int iDecLog(int x) {
		return floor(log(x) / log(10));
	}
	static int copyFile(string initialFilePath, string outputFilePath);
	static char * getTempFilename();
	static char * getStreamTempFilename();
	static int countWords(string ss);

	static FILE * myfopen(const char *path, const char *mode,
			bool lazy = false);
	static int myGetline(char **lineptr, int *n, FILE *stream);

	static inline int numberOfBranches(int numberOfTaxa) {
		return ((2 * numberOfTaxa) - 2);
	}

	static void mergeIds(t_partitionElementId & dest, t_partitionElementId id1);

	static void mergeIds(t_partitionElementId & dest, t_partitionElementId id1,
			t_partitionElementId id2);

	class PropertiesInitializer {
	public:
		PropertiesInitializer() {
			srand(time(NULL));
			ifstream read_file;
			read_file.open("partest.properties");
			string s;
			while (!read_file.eof()) {
				read_file >> s;
				config[s.substr(0, s.find('='))] = s.substr(s.find('=') + 1);
			}

		}
	};

	friend class Utilities::PropertiesInitializer;
};

} /* namespace partest */

#endif
