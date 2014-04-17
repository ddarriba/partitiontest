/*
 * PartitionTest.h
 *
 *  Created on: Apr 8, 2014
 *      Author: diego
 */

#ifndef PARTITIONTEST_H_
#define PARTITIONTEST_H_

#include "util/GlobalDefs.h"
#include <string>

using namespace std;

namespace partest {

class PartitionTest {
public:
	PartitionTest();
	virtual ~PartitionTest();

	void setDataType(DataType dataType);
	void setDoRate(bitMask doRate);
	void setIcType(InformationCriterion icType);
	void setMaxSamples(int maxSamples);
	void setOptimize(OptimizeMode optimize);
	void setSearchAlgo(SearchAlgo searchAlgo);
	void setStartingTopology(StartTopo startingTopology);

	string * getConfigFile() const;
	void setConfigFile(const char * configFile);
	string * getInputFile() const;
	void setInputFile(const char * inputFile);
	string * getOutputDir() const;
	void setOutputDir(const char * outputDir);
	string * getUserTree() const;
	void setUserTree(const char * userTree);

	bool checkParameters(void);
	bool configure(void);
private:
	string buildStartingTree(void);
};

}

#endif /* PARTITIONTEST_H_ */
