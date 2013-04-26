/*
 * ParTestOptions.cpp
 *
 *  Created on: Jan 7, 2013
 *      Author: diego
 */

#include "ParTestOptions.h"

#include <string.h>
#include <iostream>
#include <stdlib.h>

using namespace std;

namespace partest {

ParTestOptions::ParTestOptions() {
	dataType = DEFAULT_DATA_TYPE;
	alignment = 0;
	informationCriterion = DEFAULT_IC_TYPE;
	startingTopology = DEFAULT_STARTING_TOPOLOGY;
	rateVariation = 0;
	searchAlgo = DEFAULT_SEARCH_ALGO;
	sampleSize = DEFAULT_SAMPLE_SIZE;
	strcpy(treeFile, "");
	sampleSizeValue = 0.0;
}

ParTestOptions::~ParTestOptions() {
	delete alignment;
}

void ParTestOptions::set(const char *inputFile, DataType dataType,
		bitMask doRateVariation, const char *configFile,
		StartTopo startingTopology, SearchAlgo searchAlgo, InformationCriterion informationCriterion,
		SampleSize sampleSize, double sampleSizeValue, const char *userTree) {

	this->rateVariation = doRateVariation;
#ifdef _PHYML
	alignment = new PhymlAlignment(inputFile, dataType);
#else
	alignment = new PLLAlignment(inputFile, dataType);
#endif

	this->dataType = dataType;
	this->startingTopology = startingTopology;
	this->informationCriterion = informationCriterion;
	this->sampleSize = sampleSize;
	this->searchAlgo = searchAlgo;

	if (sampleSize == SS_CUSTOM) {
		this->sampleSizeValue = sampleSizeValue;
	} else {
		this->sampleSizeValue = 0.0;
	}

	strcpy(this->inputFile, inputFile);

	if (configFile) {
		strcpy(this->configFile, configFile);
	}
	if (userTree && strlen(userTree) > 0) {
		strcpy(this->treeFile, userTree);
	}

}

Alignment * ParTestOptions::getAlignment() {
	return alignment;
}

char* ParTestOptions::getConfigFile() {
	return configFile;
}

string ParTestOptions::getInputFile() const {
	return alignment->getAlignmentFile();
}

SearchAlgo ParTestOptions::getSearchAlgorithm() const {
	return searchAlgo;
}

bitMask ParTestOptions::getRateVariation() const {
	return rateVariation;
}

StartTopo ParTestOptions::getStartingTopology() const {
	return startingTopology;
}

char* ParTestOptions::getTreeFile() {
	return treeFile;
}

DataType ParTestOptions::getDataType() const {
	return dataType;
}

void ParTestOptions::setAlignment(Alignment *alignment) {
	this->alignment = alignment;
}

void ParTestOptions::setConfigFile(char* configFile) {
	strcpy(this->configFile, configFile);
}

void ParTestOptions::setInputFile(char* inputFile) {
	strcpy(this->inputFile, inputFile);
}

void ParTestOptions::setTreeFile(char* treeFile) {
	strcpy(this->treeFile, treeFile);
}

SampleSize ParTestOptions::getSampleSize(void) const {
	return sampleSize;
}

InformationCriterion ParTestOptions::getInformationCriterion(void) const {
	return informationCriterion;
}

double ParTestOptions::getSampleSizeValue(void) {
	return sampleSizeValue;
}

} /* namespace partest */
