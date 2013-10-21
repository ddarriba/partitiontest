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
 * @file ParTestOptions.cpp
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
	optimize = DEFAULT_OPTIMIZE;
	rateVariation = 0;
	searchAlgo = DEFAULT_SEARCH_ALGO;
	sampleSize = DEFAULT_SAMPLE_SIZE;
	strcpy(treeFile, "");
	sampleSizeValue = 0.0;
	resultsOutputStream = 0;
	modelsOutputStream = 0;
	partitionsOutputStream = 0;
	schemesOutputStream = 0;
	treeString = 0;
}

ParTestOptions::~ParTestOptions() {
	if (resultsOutputStream) {
		if (resultsOutputStream->is_open())
			resultsOutputStream->close();
		delete resultsOutputStream;
	}
	if (modelsOutputStream) {
		if (modelsOutputStream->is_open())
			modelsOutputStream->close();
		delete modelsOutputStream;
	}
	if (partitionsOutputStream) {
		if (partitionsOutputStream->is_open())
			partitionsOutputStream->close();
		delete partitionsOutputStream;
	}
	if (schemesOutputStream) {
		if (schemesOutputStream->is_open())
			schemesOutputStream->close();
		delete schemesOutputStream;
	}
	if (startingTopology == StartTopoFIXED) {
		delete alignment;
	}
}

void ParTestOptions::set(const char *inputFile, DataType dataType,
		bitMask doRateVariation, const char *configFile,
		StartTopo startingTopology, SearchAlgo searchAlgo, OptimizeMode optimize,
		InformationCriterion informationCriterion, SampleSize sampleSize,
		double sampleSizeValue, const char *userTree) {

	this->rateVariation = doRateVariation;

	this->dataType = dataType;
	this->startingTopology = startingTopology;
	this->informationCriterion = informationCriterion;
	this->sampleSize = sampleSize;
	this->searchAlgo = searchAlgo;
	this->optimize = optimize;

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
	this->outputFileModels = parser.getOutputFileResults();
	this->outputFileModels = parser.getOutputFileModels();
	this->outputFilePartitions = parser.getOutputFilePartitions();
	this->outputFileSchemes = parser.getOutputFileSchemes();
#ifdef _PLL
	this->pllPartitions = parser.getPllPartitions();
#endif
	resultsOutputStream = new ofstream(outputFileResults.c_str());
	modelsOutputStream = new ofstream(outputFileModels.c_str());
	partitionsOutputStream = new ofstream(outputFilePartitions.c_str());
	schemesOutputStream = new ofstream(outputFileSchemes.c_str());
#ifdef _PLL
	alignment = new PLLAlignment(inputFile, dataType, pllPartitions);
#else
	alignment = new PhymlAlignment(inputFile, dataType);
#endif
	PrintMeta::print_header(*resultsOutputStream);
	PrintMeta::print_header(*modelsOutputStream);
	PrintMeta::print_header(*partitionsOutputStream);
	PrintMeta::print_header(*schemesOutputStream);
	PrintMeta::print_options(*resultsOutputStream, *this);
	PrintMeta::print_options(*modelsOutputStream, *this);
	PrintMeta::print_options(*partitionsOutputStream, *this);
	PrintMeta::print_options(*schemesOutputStream, *this);
}

Alignment * ParTestOptions::getAlignment(void) {
	return alignment;
}

char* ParTestOptions::getConfigFile(void) {
	return configFile;
}

string ParTestOptions::getInputFile(void) const {
	return alignment->getAlignmentFile();
}

SearchAlgo ParTestOptions::getSearchAlgorithm(void) const {
	return searchAlgo;
}

OptimizeMode ParTestOptions::getOptimizeMode(void) const {
	return optimize;
}

bitMask ParTestOptions::getRateVariation(void) const {
	return rateVariation;
}

StartTopo ParTestOptions::getStartingTopology(void) const {
	return startingTopology;
}

char* ParTestOptions::getTreeFile(void) {
	return treeFile;
}

char* ParTestOptions::getTreeString(void) {
	return treeString;
}

DataType ParTestOptions::getDataType() const {
	return dataType;
}

string ParTestOptions::getOutputFileResults(void) const {
	return outputFileResults;
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

#ifdef _PLL
pllQueue *  ParTestOptions::getPllPartitions(void) const {
	return pllPartitions;
}
#endif

ofstream * ParTestOptions::getResultsOutputStream(void) const {
	return resultsOutputStream;
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

void ParTestOptions::setConfigFile(const char * configFile) {
	strcpy(this->configFile, configFile);
}

void ParTestOptions::setInputFile(const char * inputFile) {
	strcpy(this->inputFile, inputFile);
}

void ParTestOptions::setTreeFile(const char * treeFile) {
	strcpy(this->treeFile, treeFile);
}

void ParTestOptions::setTreeString(char * treeString, bool eager) {
	if (eager) {
		this->treeString = (char *) malloc(strlen(treeString) + 1);
		strcpy(this->treeString, treeString);
	} else {
		this->treeString = treeString;
	}
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

void ParTestOptions::setOutputFileResults(string outputFileResults) {
	this->outputFileResults = outputFileResults;
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
