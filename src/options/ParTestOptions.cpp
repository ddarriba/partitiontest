/*
 * ParTestOptions.cpp
 *
 *  Created on: Jan 7, 2013
 *      Author: diego
 */

#include "ParTestOptions.h"
#include "parser/ConfigParser.h"
#include "util/PrintMeta.h"
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
	modelsOutputStream = 0;
	partitionsOutputStream = 0;
	schemesOutputStream = 0;
}

ParTestOptions::~ParTestOptions() {
	if (modelsOutputStream && modelsOutputStream->is_open())
		modelsOutputStream->close();
	if (partitionsOutputStream && partitionsOutputStream->is_open())
		partitionsOutputStream->close();
	if (schemesOutputStream && schemesOutputStream->is_open())
		schemesOutputStream->close();
	delete alignment;
}

void ParTestOptions::set(const char *inputFile, DataType dataType,
		bitMask doRateVariation, const char *configFile,
		StartTopo startingTopology, SearchAlgo searchAlgo,
		InformationCriterion informationCriterion, SampleSize sampleSize,
		double sampleSizeValue, const char *userTree) {

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

	ConfigParser parser(configFile);
	this->outputFileModels = parser.getOutputFileModels();
	this->outputFilePartitions = parser.getOutputFilePartitions();
	this->outputFileSchemes = parser.getOutputFileSchemes();
	modelsOutputStream = new ofstream(outputFileModels.c_str());
	partitionsOutputStream = new ofstream(outputFilePartitions.c_str());
	schemesOutputStream = new ofstream(outputFileSchemes.c_str());
	PrintMeta::print_header(*modelsOutputStream);
	PrintMeta::print_header(*partitionsOutputStream);
	PrintMeta::print_header(*schemesOutputStream);
	PrintMeta::print_options(*modelsOutputStream, *this);
	PrintMeta::print_options(*partitionsOutputStream, *this);
	PrintMeta::print_options(*schemesOutputStream, *this);
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

string ParTestOptions::getOutputFileModels(void) const {
	return outputFileModels;
}

string ParTestOptions::getOutputFilePartitions(void) const {
	return outputFilePartitions;
}

string ParTestOptions::getOutputFileSchemes(void) const {
	return outputFileSchemes;
}

ofstream * ParTestOptions::getModelsOutputStream(void) const {
	return modelsOutputStream;
}

ofstream * ParTestOptions::getPartitionsOutputStream(void) const {
	return partitionsOutputStream;
}

ofstream * ParTestOptions::getSchemesOutputStream(void) const {
	return schemesOutputStream;
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

void ParTestOptions::setOutputFileModels(string outputFileModels) {
	this->outputFileModels = outputFileModels;
}

void ParTestOptions::setOutputFilePartitions(string outputFilePartitions) {
	this->outputFilePartitions = outputFilePartitions;
}

void ParTestOptions::setOutputFileSchemes(string outputFileSchemes) {
	this->outputFileSchemes = outputFileSchemes;
}

} /* namespace partest */
