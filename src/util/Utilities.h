#pragma once
#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include "GlobalDefs.h"
#include <map>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
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
		return (1 << x);
	}
	static inline unsigned int binaryLog(unsigned int x) {
		return ceil(log(x) / log(2));
	}
	static inline double dBinaryLog(double x) {
			return (log(x) / log(2));
		}
	static int copyFile(string initialFilePath, string outputFilePath);
	static char * getTempFilename();
	static char * getStreamTempFilename();
	static int countWords(string ss);

	static FILE * myfopen(const char *path, const char *mode, bool lazy=false);
	static int myGetline(char **lineptr, int *n, FILE *stream);

	static inline int numberOfBranches(int numberOfTaxa) {
		return ((2 * numberOfTaxa) - 1);
	}

	class PropertiesInitializer {
	public:
		PropertiesInitializer() {
			srand (time(NULL));ifstream read_file;
			read_file.open("partest.properties");
			string s;
			while ( !read_file.eof() ) {
				read_file >> s;
				config[s.substr(0, s.find('='))] = s.substr(s.find('=') + 1);
			}

		}
	};

	friend class Utilities::PropertiesInitializer;
};

} /* namespace partest */

#endif
