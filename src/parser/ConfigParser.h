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
 * @file ConfigParser.h
 * @author Diego Darriba
 *
 * @brief Utility for parsing configuration files.
 */

#ifndef CONFIGPARSER_H_
#define CONFIGPARSER_H_

#define MAX_SECTIONS 20
#include "util/GlobalDefs.h"
#include <vector>

#define PARTITIONS_TAG     "partitions"
#define SCHEMES_TAG        "schemes"

#define SEARCH_TAG             "search"
#define SEARCH_ALGORITHM_TAG   "algorithm"
#define SEARCH_ALGORITHM_REPS  "replicates"

#define MODELS_TAG          "models"
#define MODELS_INCLUDE_TAG  "include"
#define MODELS_EPSILON_TAG  "epsilon"
#define MODELS_DO_F_TAG     "include-f"
#define MODELS_DO_M_TAG     "include-m"
#define MODELS_DO_I_TAG     "include-i"
#define MODELS_DO_G_TAG     "include-g"
#define MODELS_DO_IG_TAG    "include-ig"

#define INPUT_TAG              "input"
#define INPUT_MSA_TAG          "msa"
#define INPUT_TREE_TAG         "tree"
#define INPUT_STARTTOPO_TAG    "topo"
#define INPUT_KEEPBRANCHES_TAG "keep-branches"
#define INPUT_DATATYPE_TAG     "datatype"

#define OUTPUT_TAG           "output"
#define OUTPUT_TMP_PATH      "tmp"
#define OUTPUT_BASE_PATH     "path"
#define OUTPUT_RESULTS_TAG   "results"
#define OUTPUT_MODELS_TAG    "models"
#define OUTPUT_PARTS_TAG     "partitions"
#define OUTPUT_SCHEMES_TAG   "schemes"
#ifdef _WIN32
#define DEFAULT_OUTPUT_TMP_PATH   ""
#else
#define DEFAULT_OUTPUT_TMP_PATH      "/tmp"
#endif
#define DEFAULT_OUTPUT_BASE_PATH     ""
#define DEFAULT_OUTPUT_RESULTS_TAG   "results.out"
#define DEFAULT_OUTPUT_MODELS_TAG    "models.out"
#define DEFAULT_OUTPUT_PARTS_TAG     "partitions.out"
#define DEFAULT_OUTPUT_SCHEMES_TAG   "schemes.out"

namespace partest
{

  /**
   * @brief Structure with information about a single-gene partition
   */
  struct partitionInfo
  {
    t_partitionElementId partitionId; /** Partition id */
    int start[MAX_SECTIONS]; /** Starting position */
    int end[MAX_SECTIONS]; /** Ending position */
    int stride[MAX_SECTIONS]; /** Stride for codon position (0 means no codon division) */
    int numberOfSections;
    std::string name; /** Name of the gene/partition */
    ~partitionInfo (void)
    {
    }
  };

  /**
   * @brief Utility class for parsing configuration files.
   */
  class ConfigParser
  {
  public:
    /**
     * @brief Creates a new ConfigParser
     *
     * @param configFile Configuration file name.
     */
    ConfigParser (const char * configFile);

    virtual ~ConfigParser ();

    /**
     * @brief Gets the whole set of partitions.
     *
     * @return The set of partitions.
     */
    std::vector<partitionInfo> * getPartitions ();

    /**
     * @brief Gets a single partition.
     *
     * @param index The index of the partition to be retrieved
     *
     * @return The partition at index
     */
    struct partitionInfo getPartition (size_t index);

    /**
     * @brief Gets the number of partitions.
     *
     * @return The number of partitions.
     */
    size_t getNumberOfPartitions ()
    {
      return number_of_genes;
    }

    /**
     * @brief Prints format instructions for configuration files.
     */
    static void printFormat ();

    /**
     * @brief Creates a template for configuration files.
     */
    static void createTemplate ();

    /**
     * @brief Gets the file name for results output.
     */
    const std::string& getOutputFileResults () const
    {
      return outputFileResults;
    }

    /**
     * @brief Gets the file name for model selections output.
     */
    const std::string& getOutputFileModels () const
    {
      return outputFileModels;
    }

    /**
     * @brief Gets the file name for partition selections output.
     */
    const std::string& getOutputFilePartitions () const
    {
      return outputFilePartitions;
    }

    /**
     * @brief Gets the file name for scheme selections output.
     */
    const std::string& getOutputFileSchemes () const
    {
      return outputFileSchemes;
    }

    /**
     * @brief Gets the base path for output files.
     */
    const std::string& getOutputBasePath () const
    {
      return outputBasePath;
    }

    /**
     * @brief Gets the path for temporary files.
     */
    const std::string& getOutputTmpPath () const
    {
      return outputTmpPath;
    }

    const char * getInputFile (void) const
    {
      return inputFile;
    }

    const char * getUserTree (void) const
    {
      return userTree;
    }

    void createPartitions ();

  private:

    /**
     * @brief Parses a partitioning scheme from a line definition.
     *
     * @param line The line to be parsed.
     * @param[out] scheme The partitioning scheme.
     *
     * @return 0 if there was no error.
     */
    int parseScheme (char * line, t_partitioningScheme * scheme);

    /**
     * @brief Parses the complete details of a partition.
     *
     * @param line The line to be parsed.
     * @param[out] pInfo The partitionInfo structure.
     *
     * @return 0 if there was no error.
     */
    int parsePartitionDetails (char * line, struct partitionInfo * pInfo);

    void createSinglePartition ();

    const char * configFile; /** Configuration file name */
    char inputFile[256]; /** User input alignment */
    char userTree[256]; /** User input tree */
    std::vector<partitionInfo> * partitions; /** Vector of partitions */

    /** search algorithm **/
    std::string outputFileResults; /** File name for results output */
    std::string outputFileModels; /** File name for model selections output */
    std::string outputFilePartitions; /** File name for partition selections output */
    std::string outputFileSchemes; /** File name for scheme selections output */
    std::string outputBasePath; /** Base path for output files */
    std::string outputTmpPath; /** Base path for output files */
  };

} /* namespace partest */
#endif /* CONFIGPARSER_H_ */
