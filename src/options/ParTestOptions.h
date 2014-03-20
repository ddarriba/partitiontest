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
 * @file ParTestOptions.h
 */

#ifndef PARTESTOPTIONS_H_
#define PARTESTOPTIONS_H_

#ifdef _PHYML
#include "indata/PhymlAlignment.h"
#else
#include "indata/PLLAlignment.h"
#endif
#include "util/GlobalDefs.h"
#include <iostream>
#include <fstream>
#include <string>

namespace partest {

/**
 * @brief Main PartitionTest execution options
 */
class ParTestOptions {
public:

	/**
	 * @brief Default constructor for empty options instance.
	 */
	ParTestOptions();
	virtual ~ParTestOptions();

	/**
	 * @brief Sets the PartitionTest execution parameters.
	 *
	 * @param inputFile The alignment input file.
	 * @param dataType The alignment data type.
	 * @param rateVariation The rate variation mask for building candidate models set.
	 * @param configFile The configuration file name.
	 * @param startingTopology The starting topology for each model.
	 * @param searchAlgo The search algorithm.
	 * @param maxSamples The maximum number of samples for the HCluster algorithm.
	 * @param optimize The mode for optimizing the best-fit partition.
	 * @param informationCriterion The IC to be used during the search process.
	 * @param sampleSize The sample size to be used during the model selection process.
	 * @param sampleSizeValue Custom sample size value if applicable (optional).
	 * @para userTree User tree file name (optional).
	 */
	void set(const char *inputFile, DataType dataType, bitMask rateVariation,
			const char *configFile, StartTopo startingTopology,
			SearchAlgo searchAlgo, int maxSamples, OptimizeMode optimize,
			InformationCriterion informationCriterion, SampleSize sampleSize,
			double sampleSizeValue, const char *userTree, const char *outputDir);

	/**
	 * @brief Gets the alignment.
	 *
	 * @return The alignment instance.
	 */
	Alignment * getAlignment(void);

	/**
	 * @brief Gets the configuration file name.
	 *
	 * @return The configuration file name.
	 */
	char* getConfigFile(void);

	/**
	 * @brief Gets the execution data type.
	 *
	 * @return The data type (NUCLEIC / PROTEIC).
	 */
	DataType getDataType(void) const;

	/**
	 * @brief Gets the input alignment file name.
	 *
	 * @return The input alignment file name.
	 */
	string getInputFile(void) const;

	/**
	 * @brief Gets the rate variation mask for building candidate models set.
	 *
	 * @return The rate variation mask.
	 */
	bitMask getRateVariation(void) const;

	/**
	 * @brief Gets the starting topology to be used for optimizing each model.
	 *
	 * @return The starting topology.
	 */
	StartTopo getStartingTopology(void) const;

	/**
	 * @brief Gets the tree file name to be used as starting topology.
	 *
	 * @return The input tree file name.
	 */
	char * getTreeFile(void);

	/**
	 * @brief Gets the input tree in Newick format.
	 *
	 * @return The input tree.
	 */
	char * getTreeString(void);

	/**
	 * @brief Gets the search algorithm.
	 *
	 * @return The search algorithm.
	 */
	SearchAlgo getSearchAlgorithm(void) const;

	/**
	 * @brief Gets the maximum number of samples for the HCluster algorithm.
	 *
	 * @return The maximum number of samples.
	 */
	int getMaxSamples(void) const;

	/**
	 * @brief Gets the best-fit partition optimization mode.
	 *
	 * @return The best-fit partition optimization mode.
	 */
	OptimizeMode getOptimizeMode(void) const;

	/**
	 * @brief Gets the sample size type for model selection.
	 *
	 * @return The sample size type.
	 */
	SampleSize getSampleSize(void) const;

	/**
	 * @brief Gets the criterion for model selection.
	 *
	 * @return The information criterion.
	 */
	InformationCriterion getInformationCriterion(void) const;

	/**
	 * @brief Gets the custom sample size value.
	 *
	 * @return The custom sample size value.
	 */
	double getSampleSizeValue(void);

	/**
	 * @brief Gets the output file for results.
	 *
	 * @return The output file for results.
	 */
	string getOutputFileResults(void) const;

	/**
	 * @brief Gets the output file for model selections in partitions.
	 *
	 * @return The output file for model selections in partitions.
	 */
	string getOutputFileModels(void) const;

	/**
	 * @brief Gets the output file for partition selections.
	 *
	 * @return The output file for partition selections.
	 */
	string getOutputFilePartitions(void) const;

	/**
	 * @brief Gets the output file for scheme selections.
	 *
	 * @return The output file for scheme selections.
	 */
	string getOutputFileSchemes(void) const;

#ifdef _PLL
	/**
	 * @brief Gets the partition definitions for PLL.
	 */
	pllQueue * getPllPartitions(void) const;
#endif

	/**
	 * @brief Gets the output stream for results.
	 *
	 * @return The output file for results.
	 */
	ofstream * getResultsOutputStream(void) const;

	/**
	 * @brief Gets the output stream for model selections in partitions.
	 *
	 * @return The output file for model selections in partitions.
	 */
	ofstream * getModelsOutputStream(void) const;

	/**
	 * @brief Gets the output file for partition selections.
	 *
	 * @return The output file for partition selections.
	 */
	ofstream * getPartitionsOutputStream(void) const;

	/**
	 * @brief Gets the output file for scheme selections.
	 *
	 * @return The output file for partition selections.
	 */
	ofstream * getSchemesOutputStream(void) const;

	/**
		 * @brief Gets the path for temporary files.
		 *
		 * @return The path for temporary files.
		 */
	string getOutputTmpPath(void) const;

	/**
	 * @brief Sets the input alignment.
	 *
	 * @param alignment The input alignment.
	 */
	void setAlignment(Alignment *alignment);

	/**
	 * @brief Sets the configuration file.
	 *
	 * @param configFile The configuration file name.
	 */
	void setConfigFile(const char* configFile);

	/**
	 * @brief Sets the input alignment file.
	 *
	 * @param inputFile The input alignment filename.
	 */
	void setInputFile(const char* inputFile);

	/**
	 * @brief Sets the input tree file.
	 *
	 * @param treeFile The input tree file name.
	 */
	void setTreeFile(const char* treeFile);

	/**
	 * @brief Sets the input tree in Newick format.
	 *
	 * @param treeString The input tree.
	 * @param eager If false, the variable maintains just a reference to the tree
	 */
	void setTreeString(char* treeString, bool eager = false);

	/**
	 * @brief Sets the criterion for selection.
	 *
	 * @param informationCriterion The information criterion for selection.
	 */
	void setInformationCriterion(InformationCriterion informationCriterion);

	/**
	 * @brief Sets the sample size type.
	 *
	 * @param sampleSize The sample size type.
	 */
	void setSampleSize(SampleSize sampleSize);

	/**
	 * @brief Sets the custom sample size value.
	 *
	 * @param sampleSizeValue The custom sample size value.
	 */
	void setSampleSizeValue(double sampleSizeValue);

	/**
	 * @brief Sets the file name for model selections output.
	 *
	 * @param outputFileModels The file name for model selections output.
	 */
	void setOutputFileResults(string outputFileResults);

	/**
	 * @brief Sets the file name for model selections output.
	 *
	 * @param outputFileModels The file name for model selections output.
	 */
	void setOutputFileModels(string outputFileModels);

	/**
	 * @brief Sets the file name for partition selections output.
	 *
	 * @param outputFilePartitions The file name for partition selections output.
	 */
	void setOutputFilePartitions(string outputFilePartitions);

	/**
	 * @brief Sets the file name for scheme selections output.
	 *
	 * @param outputFileSchemes The file name for scheme selections output.
	 */
	void setOutputFileSchemes(string outputFileSchemes);

private:
	Alignment * alignment; /** Input MSA */
	bitMask rateVariation; /** Rates variation to analyze */
	DataType dataType; /** Whether input data is nucleic or proteic */
	StartTopo startingTopology; /** Starting topology for optimization */
	SearchAlgo searchAlgo; /** Search algorithm */
	int maxSamples; /** Number of samples for the HCluster algorithm */
	OptimizeMode optimize; /** The best-fit partition optimize mode */
	SampleSize sampleSize; /** Sample size mode for model selection */
	double sampleSizeValue; /** Sample size value for model selection */
	InformationCriterion informationCriterion; /** Statistical criterion for model selection */
	char inputFile[256]; /** Input MSA file */
	char configFile[256]; /** Input configuration file */
	char treeFile[256]; /** Input tree file */
	char * treeString; /** Input tree in Newick format */
	char outputDir[256]; /** Directory for output files */
	string outputFileResults; /** Output results file */
	string outputFileModels; /** Output models file */
	string outputFilePartitions; /** Output partitions file */
	string outputFileSchemes; /** Output schemes file */
	string outputTmpPath; /** Output temporary files path */
#ifdef _PLL
	pllQueue * pllPartitions; /** Partitions definition for PLL */
#endif
	ofstream * resultsOutputStream; /** Output results stream */
	ofstream * modelsOutputStream; /** Output models stream */
	ofstream * partitionsOutputStream; /** Output partitions stream */
	ofstream * schemesOutputStream; /** Output schemes stream */
};

} /* namespace partest */

#endif /* PARTESTOPTIONS_H_ */
