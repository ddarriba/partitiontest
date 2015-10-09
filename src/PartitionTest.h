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
 * @file PartitionTest.h
 * @author Diego Darriba
 * @brief Main class of PartitionTest
 */

#ifndef PARTITIONTEST_H_
#define PARTITIONTEST_H_

#include "util/GlobalDefs.h"
#include <string>

namespace partest {

class PartitionTest {
public:
	PartitionTest();
	virtual ~PartitionTest();

	void setDataType(DataType dataType);
	void setDoRate(bitMask doRate);
	void setIcType(InformationCriterion icType);
	void setMaxSamples(int maxSamples);
	void setSamplesPercent(double samplesPercent);
	void setOptimize(OptimizeMode optimize);
	void setSearchAlgo(SearchAlgo searchAlgo);
	void setStartingTopology(StartTopo startingTopology);

	std::string * getConfigFile() const;
	void setConfigFile(const char * configFile);
	std::string * getInputFile() const;
	void setInputFile(const char * inputFile);
	std::string * getOutputDir() const;
	void setOutputDir(const char * outputDir);
	std::string * getUserTree() const;
	void setUserTree(const char * userTree);

	bool checkParameters(void);
	bool configure(PartitionTest * ptest, int argc, char * argv[]);
private:
	std::string buildStartingTree(void);
};

}

#endif /* PARTITIONTEST_H_ */
