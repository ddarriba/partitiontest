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
#include <iostream>
#include <fstream>
#include <string>

namespace partest {

class ParTestOptions {
public:
	ParTestOptions();
	virtual ~ParTestOptions();
	void set(const char *inputFile, DataType dataType, bitMask rateVariation,
			const char *configFile, StartTopo startingTopology,
			SearchAlgo searchAlgo, InformationCriterion informationCriterion,
			SampleSize sampleSize, double sampleSizeValue = 0.0,
			const char *user_tree = "");

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
	string getOutputFileModels(void) const;
	string getOutputFilePartitions(void) const;
	string getOutputFileSchemes(void) const;
	ofstream * getModelsOutputStream(void) const;
	ofstream * getPartitionsOutputStream(void) const;
	ofstream * getSchemesOutputStream(void) const;

	void setAlignment(Alignment *alignment);
	void setConfigFile(char* configFile);
	void setInputFile(char* inputFile);
	void setTreeFile(char* treeFile);
	void setInformationCriterion(InformationCriterion informationCriterion);
	void setSampleSize(SampleSize sampleSize);
	void setSampleSizeValue(double sampleSizeValue);
	void setOutputFileModels(string outputFileModels);
	void setOutputFilePartitions(string outputFilePartitions);
	void setOutputFileSchemes(string outputFileSchemes);

private:
	Alignment * alignment;  /** Input MSA */
	bitMask rateVariation;  /** Rates variation to analyze */
	DataType dataType;      /** Whether input data is nucleic or proteic */
	StartTopo startingTopology;	/** Starting topology for optimization */
	SearchAlgo searchAlgo;		/** Search algorithm */
	SampleSize sampleSize;		/** Sample size mode for model selection */
	double sampleSizeValue;		/** Sample size value for model selection */
	InformationCriterion informationCriterion;	/** Statistical criterion for model selection */
	char inputFile[256];	/** Input MSA file */
	char configFile[256];	/** Input configuration file */
	char treeFile[256];		/** Input tree file */
	string outputFileModels;		/** Output models file */
	string outputFilePartitions;	/** Output partitions file */
	string outputFileSchemes;		/** Output schemes file */
	ofstream * modelsOutputStream;		/** Output models stream */
	ofstream * partitionsOutputStream;	/** Output partitions stream */
	ofstream * schemesOutputStream;		/** Output schemes stream */
};

} /* namespace partest */

#endif /* PARTESTOPTIONS_H_ */
