/*
 * ParTestOptions.h
 *
 *  Created on: Jan 7, 2013
 *      Author: diego
 */

#ifndef PARTESTOPTIONS_H_
#define PARTESTOPTIONS_H_

#ifdef _PHYML
#include "../indata/PhymlAlignment.h"
#else
#include "../indata/PLLAlignment.h"
#endif
#include "../util/GlobalDefs.h"
#include <string>

namespace partest {

class ParTestOptions {
public:
	ParTestOptions();
	virtual ~ParTestOptions();
	void set(const char *inputFile, DataType dataType, bitMask rateVariation,
			const char *configFile, StartTopo startingTopology, SearchAlgo searchAlgo, InformationCriterion informationCriterion,
			SampleSize sampleSize, double sampleSizeValue = 0.0, const char *user_tree = "");

	Alignment * getAlignment(void);
	char* getConfigFile(void);
	DataType getDataType(void) const;
	string getInputFile(void) const;
	bitMask getRateVariation(void) const;
	StartTopo getStartingTopology(void) const;
	char* getTreeFile(void);
	SearchAlgo getSearchAlgorithm(void) const;
	SampleSize getSampleSize(void) const;
	InformationCriterion getInformationCriterion(void) const;
	double getSampleSizeValue(void);

	void setAlignment(Alignment *alignment);
	void setConfigFile(char* configFile);
	void setInputFile(char* inputFile);
	void setTreeFile(char* treeFile);
	void setInformationCriterion(InformationCriterion informationCriterion);
	void setSampleSize(SampleSize sampleSize);
	void setSampleSizeValue(double sampleSizeValue);

private:
	Alignment * alignment;
	bitMask rateVariation;
	DataType dataType;
	StartTopo startingTopology;
	SearchAlgo searchAlgo;
	SampleSize sampleSize;
	double sampleSizeValue;
	InformationCriterion informationCriterion;
	char inputFile[256];
	char configFile[256];
	char treeFile[256];
};

} /* namespace partest */

#endif /* PARTESTOPTIONS_H_ */
