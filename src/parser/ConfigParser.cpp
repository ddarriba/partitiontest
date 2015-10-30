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
 * @file ConfigParser.cpp
 * @author Diego Darriba
 */

#include "ConfigParser.h"
#include "util/Utilities.h"
#include "INIReader.h"

#include <pll/parsePartition.h>
#include <string.h>
#include <assert.h>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <assert.h>

using namespace std;

namespace partest
{

  struct comparePartitionInfos
  {
    inline bool operator() (partitionInfo p1, partitionInfo p2)
    {

      if (p1.numberOfSections < 1)
        return true;
      if (p2.numberOfSections < 1)
        return false;
      if (p1.start[0] != p2.start[0])
        return p1.start[0] < p2.start[0];
      else
        return p1.stride[0] < p2.stride[0];
    }
  };

  static int checkOverlap (vector<partitionInfo> * parts, int maxSite)
  {
    /* boolean array for marking that a site was assigned a partition */
    char * used = (char *) calloc ((size_t) maxSite, sizeof(char));

    for (size_t i = 0; i < parts->size (); i++)
    {
      partitionInfo partInfo = parts->at (i);
      for (int sec = 0; sec < partInfo.numberOfSections; sec++)
      {
        for (int site = partInfo.start[sec]; site <= partInfo.end[sec]; site +=
            partInfo.stride[sec])
        {
          if (used[site])
          {
            free (used);
            cerr << "[ERROR] Detected overlapping partitions: " << endl;
            cerr << "        " << partInfo.start[sec] << "-"
                << partInfo.end[sec];
            if (partInfo.stride[sec] > 1)
            {
              cerr << "/" << partInfo.stride[sec];
            }
            cerr << endl;
            exit_partest (EX_IOERR);
          }
          used[site] = 1;
        }
      }
    }

    free (used);
    return 0;
  }

  void ConfigParser::createSinglePartition ()
  {
    if (!input_file)
    {
      cerr << "[ERROR] Input file has not been set" << endl;
      exit_partest (EX_CONFIG);
    }

    number_of_genes = 1;
    partitions = new vector<partitionInfo> (number_of_genes);
    singleGeneNames = (string **) malloc (sizeof(string*));
    singleGeneNames[0] = new string ("SinglePartition");

    phylip = pllParseAlignmentFile (PLL_FORMAT_PHYLIP, input_file->c_str ());
    if (!phylip)
    {
      cerr << "[ERROR] There was an error parsing input data \""
          << (*input_file) << "\"" << endl;
      exit_partest (EX_IOERR);
    }
    num_taxa = (size_t) phylip->sequenceCount;
    seq_len = (size_t) phylip->sequenceLength;

    pllPartitionRegion * pregion;
    pllPartitionInfo * pinfo;
    pllQueueInit (&pllPartsQueue);

    partitions->at (0).partitionId.push_back (0);
    pinfo = (pllPartitionInfo *) malloc (sizeof(pllPartitionInfo));
    pinfo->ascBias = PLL_FALSE;
    pllQueueInit (&(pinfo->regionList));
    pllQueueAppend (pllPartsQueue, (void *) pinfo);
    pinfo->partitionName = (char *) malloc (2);
    strcpy (pinfo->partitionName, "S");
    pinfo->partitionModel = (char *) malloc (1);

    switch (data_type)
      {
      case DT_NUCLEIC:
        pinfo->protModels = -1;
        pinfo->protUseEmpiricalFreqs = -1;
        pinfo->dataType = PLL_DNA_DATA;
        pinfo->optimizeBaseFrequencies = PLL_TRUE;
        break;
      case DT_PROTEIC:
        pinfo->protModels = PLL_AUTO;
        pinfo->protUseEmpiricalFreqs = PLL_FALSE;
        pinfo->dataType = PLL_AA_DATA;
        pinfo->optimizeBaseFrequencies = PLL_TRUE;
        break;
      default:
        assert(0);
        break;
      }
    pregion = (pllPartitionRegion *) malloc (sizeof(pllPartitionRegion));
    pregion->start = 1;
    pregion->end = (int) seq_len;
    pregion->stride = 1;
    pllQueueAppend (pinfo->regionList, (void *) pregion);
  }

  void ConfigParser::createPartitions ()
  {
    if (partitions)
    {
      pllPartitionRegion * pregion;
      pllPartitionInfo * pinfo;

      pllQueueInit (&pllPartsQueue);

      for (size_t i = 0; i < number_of_genes; i++)
      {
        partitions->at (i).partitionId.push_back (i);
        pinfo = (pllPartitionInfo *) malloc (sizeof(pllPartitionInfo));
        pinfo->ascBias = PLL_FALSE;
        pllQueueInit (&(pinfo->regionList));
        pllQueueAppend (pllPartsQueue, (void *) pinfo);

        pinfo->partitionName = (char *) malloc (
            (partitions->at (i).name.size () + 1) * sizeof(char));
        strcpy (pinfo->partitionName, partitions->at (i).name.c_str ());
        pinfo->partitionModel = (char *) malloc (1);

        switch (data_type)
          {
          case DT_NUCLEIC:
            pinfo->protModels = -1;
            pinfo->protUseEmpiricalFreqs = -1;
            pinfo->dataType = PLL_DNA_DATA;
            pinfo->optimizeBaseFrequencies = PLL_TRUE;
            break;
          case DT_PROTEIC:
            pinfo->protModels = PLL_AUTO;
            pinfo->protUseEmpiricalFreqs = PLL_FALSE;
            pinfo->dataType = PLL_AA_DATA;
            pinfo->optimizeBaseFrequencies = PLL_TRUE;
            break;
          default:
            assert(0);
            break;
          }

        for (int j = 0; j < partitions->at (i).numberOfSections; j++)
        {
          pregion = (pllPartitionRegion *) malloc (sizeof(pllPartitionRegion));
          pregion->start = partitions->at (i).start[j];
          pregion->end = partitions->at (i).end[j];
          pregion->stride = partitions->at (i).stride[j];
          pllQueueAppend (pinfo->regionList, (void *) pregion);
        }
      }
    }
    else
    {
      createSinglePartition ();
    }
  }

  ConfigParser::ConfigParser (const char * _configFile) :
      configFile (_configFile), partitions (0)
  {

    if (configFile != 0 && strcmp (configFile, ""))
    {

      const char * value;
      strcpy (inputFile, "");
      strcpy (userTree, "");
      //IniParser ini(configFile);
      INIReader ini (configFile);

      /** INPUT **/
      value = ini.Get (INPUT_TAG, INPUT_MSA_TAG, "").c_str ();
      if (strcmp (value, ""))
      {
        size_t prefixLen = 0;
        if (value[0] != '/')
        {
          char * basedir = (char *) malloc (strlen (configFile) + 1);
          strcpy (basedir, configFile);
          char * lastSlash = strrchr (basedir, char_separator);
          if (lastSlash)
          {
            *(lastSlash + 1) = 0;
            strcpy (inputFile, basedir);
            prefixLen = strlen (basedir);
          }
          free (basedir);
        }
        strcpy (inputFile + prefixLen, value);
      }

      value = ini.Get (INPUT_TAG, INPUT_STARTTOPO_TAG, "").c_str ();
      if (!strcmp (value, "fixed"))
      {
        starting_topology = StartTopoFIXED;
      }
      else if (!strcmp (value, "ml"))
      {
        starting_topology = StartTopoML;
      }
      else if (!strcmp (value, "mp"))
      {
        starting_topology = StartTopoMP;
      }
      else if (!strcmp (value, "user"))
      {
        starting_topology = StartTopoUSER;
      }
      else
      {
        starting_topology = DEFAULT_STARTING_TOPOLOGY;
      }

      value = ini.Get (INPUT_TAG, INPUT_TREE_TAG, "").c_str ();
      if (strcmp (value, ""))
      {
        strcpy (userTree, value);
      }

      value = ini.Get (INPUT_TAG, INPUT_DATATYPE_TAG, "nt").c_str ();
      if (!strcmp (value, "aa"))
      {
        data_type = DT_PROTEIC;
      }
      else if (!strcmp (value, "nt"))
      {
        data_type = DT_NUCLEIC;
      }
      else
      {
        cerr << "[ERROR] \"" << INPUT_DATATYPE_TAG << " = " << value
            << "\" is not a valid value. Data type should be \"nt\" or \"aa\"."
            << endl;
        exit_partest (EX_CONFIG);
      }

      value = ini.Get (INPUT_TAG, INPUT_KEEPBRANCHES_TAG, "false").c_str ();
      if (!strcmp (value, "true") && reoptimize_branch_lengths == true)
      {
        reoptimize_branch_lengths = false;
      }

      /** CANDIDATE MODELS **/
      value = ini.Get (MODELS_TAG, MODELS_INCLUDE_TAG, "all").c_str ();
      if (!strcmp (value, "all"))
      {
        optimize_mode = OPT_SEARCH;
      }
      else if (!strcmp (value, "gtr"))
      {
        optimize_mode = OPT_GTR;
      }
      else
      {
        // parse protein matrices
        protModels = 0;
        istringstream iss (value);
        string curMatrix;
        while (iss)
        {
          size_t curValue = 0;
          iss >> curMatrix;
          if (!strcmp (curMatrix.c_str (), "dayhoff"))
          {
            curValue = PROT_MATRIX_DAYHOFF;
          }
          else if (!strcmp (curMatrix.c_str (), "dcmut"))
          {
            curValue = PROT_MATRIX_DCMUT;
          }
          else if (!strcmp (curMatrix.c_str (), "jtt"))
          {
            curValue = PROT_MATRIX_JTT;
          }
          else if (!strcmp (curMatrix.c_str (), "mtrev"))
          {
            curValue = PROT_MATRIX_MTREV;
          }
          else if (!strcmp (curMatrix.c_str (), "wag"))
          {
            curValue = PROT_MATRIX_WAG;
          }
          else if (!strcmp (curMatrix.c_str (), "cprev"))
          {
            curValue = PROT_MATRIX_CPREV;
          }
          else if (!strcmp (curMatrix.c_str (), "rtrev"))
          {
            curValue = PROT_MATRIX_RTREV;
          }
          else if (!strcmp (curMatrix.c_str (), "vt"))
          {
            curValue = PROT_MATRIX_VT;
          }
          else if (!strcmp (curMatrix.c_str (), "blosum62"))
          {
            curValue = PROT_MATRIX_BLOSUM62;
          }
          else if (!strcmp (curMatrix.c_str (), "mtmam"))
          {
            curValue = PROT_MATRIX_MTMAM;
          }
          else if (!strcmp (curMatrix.c_str (), "mtart"))
          {
            curValue = PROT_MATRIX_MTART;
          }
          else if (!strcmp (curMatrix.c_str (), "hivb"))
          {
            curValue = PROT_MATRIX_HIVB;
          }
          else if (!strcmp (curMatrix.c_str (), "hivw"))
          {
            curValue = PROT_MATRIX_HIVW;
          }
          else if (!strcmp (curMatrix.c_str (), "mtzoa"))
          {
            curValue = PROT_MATRIX_MTZOA;
          }
          else if (!strcmp (curMatrix.c_str (), "pmb"))
          {
            curValue = PROT_MATRIX_PMB;
          }
          else if (!strcmp (curMatrix.c_str (), "flu"))
          {
            curValue = PROT_MATRIX_FLU;
          }
          else if (!strcmp (curMatrix.c_str (), "auto"))
          {
            curValue = PROT_MATRIX_AUTO;
          }
          protModels |= Utilities::binaryPow (curValue);
        }
        optimize_mode = OPT_CUSTOM;
      }

      /** OPTIMIZATION EPSILON THRESHOLD **/
      value = ini.Get (MODELS_TAG, MODELS_EPSILON_TAG, "n/a").c_str ();
      if (strcasecmp (value, "n/a"))
      {
        if (Utilities::isNumeric (value))
        {
          epsilon = atof (value);
          if (epsilon <= 0.0f)
          {
            cerr << "[ERROR] \"" << MODELS_EPSILON_TAG << " = " << value
                << "\" is not a valid value. Epsilon should be \"auto\", or a numeric value greater or equal than 0."
                << endl;
            exit_partest (EX_CONFIG);
          }
        }
        else
        {
          if (!strcasecmp (value, "auto"))
          {
            epsilon = AUTO_EPSILON;
          }
          else
          {
            cerr << "[ERROR] \"" << MODELS_EPSILON_TAG << " = " << value
                << "\" is not a valid value. Epsilon should be \"auto\", or a numeric value greater or equal than 0."
                << endl;
            exit_partest (EX_CONFIG);
          }
        }
      }

      /** SEARCH ALGORITHM **/
      value = ini.Get (SEARCH_TAG, SEARCH_ALGORITHM_TAG, "auto").c_str ();
      if (!strcmp (value, "k1"))
      {
        search_algo = SearchK1;
      }
      else if (!strcmp (value, "kn"))
      {
        search_algo = SearchKN;
      }
      else if (!strcmp (value, "greedy"))
      {
        search_algo = SearchGreedy;
      }
      else if (!strcmp (value, "exhaustive"))
      {
        search_algo = SearchExhaustive;
      }
      else if (!strcmp (value, "random"))
      {
        search_algo = SearchRandom;
      }
      else if (!strcmp (value, "hcluster"))
      {
        search_algo = SearchHCluster;
      }
      else if (!strcmp (value, "auto"))
      {
        search_algo = SearchAuto;
      }
      else
      {
        cerr << "Invalid search algorithm \"" << value << "\"" << endl;
        exit_partest (EX_CONFIG);
      }

      string samples_str = ini.Get (SEARCH_TAG, SEARCH_ALGORITHM_REPS, "1");
      if (Utilities::isInteger (samples_str.c_str ()))
      {
        /* apply absolute number of replicates */
        max_samples = (int) ini.GetInteger (SEARCH_TAG,
        SEARCH_ALGORITHM_REPS,
                                            1);
        samples_percent = 0.0;
      }
      else
      {
        /* apply percentage of replicates */
        max_samples = 0;
        if (samples_str[samples_str.length () - 1] == '%')
        {
          samples_str[samples_str.length () - 1] = '\0';
          if (Utilities::isInteger (samples_str.c_str ()))
          {
            samples_percent = atof (samples_str.c_str ()) / 100;
          }
          else
          {
            cerr << "Invalid samples percent " << samples_percent << endl;
            exit_partest (EX_CONFIG);
          }
        }
        else
        {
          if (Utilities::isNumeric (samples_str.c_str ()))
          {
            samples_percent = atof (samples_str.c_str ());
          }
          else
          {
            cerr << "Invalid samples percent " << samples_percent << endl;
            exit_partest (EX_CONFIG);
          }
        }
      }

      /** PARTITIONS **/
      std::vector<std::string> * keys = ini.getGenes (
      PARTITIONS_TAG);
      number_of_genes = keys->size () / 2;

      if (number_of_genes > PLL_NUM_BRANCHES)
      {
        cerr << "[ERROR] PLL is compiled for a maximum of " << PLL_NUM_BRANCHES
            << " partitions." << endl;
        cerr << "        You need to increase PLL_NUM_BRANCHES in pll.h and "
            << "recompile both PLL and PartitionTest." << endl;
        exit_partest (EX_CONFIG);
      }

      if (number_of_genes)
      {
        partitions = new vector<partitionInfo> (number_of_genes);
        singleGeneNames = (string **) calloc ((size_t) number_of_genes,
                                              sizeof(string*));
        int maxSite = 0;

        for (size_t partitionId = 0; partitionId < number_of_genes;
            partitionId++)
        {
          partitions->at (partitionId).name = keys->at (partitionId * 2);
          singleGeneNames[partitionId] = new string (
              keys->at (partitionId * 2));

          char * lineBuffer = (char *) malloc (
              keys->at (partitionId * 2 + 1).length () + 1);
          strcpy (lineBuffer, keys->at (partitionId * 2 + 1).c_str ());

          parsePartitionDetails (lineBuffer, &partitions->at (partitionId));
          for (int sec = 0; sec < partitions->at (partitionId).numberOfSections;
              sec++)
          {
            if (partitions->at (partitionId).end[sec] > maxSite)
            {
              maxSite = partitions->at (partitionId).end[sec];
            }
          }

          free (lineBuffer);
        }

        checkOverlap (partitions, maxSite);
        delete keys;

        std::sort (partitions->begin (), partitions->end (),
                   comparePartitionInfos ());
      }

      /** SCHEMES **/
      vector<string> * defSchemes = ini.getSchemes (SCHEMES_TAG);
      number_of_schemes = defSchemes->size ();
      if (number_of_schemes)
      {
        schemes = new vector<t_partitioningScheme> (number_of_schemes);

        size_t schemeId = 0;
        for (size_t i = 0; i < defSchemes->size (); i++)
        {
          string scheme = defSchemes->at (i);
          cout << "SCHEME " << scheme << endl;
          char * lineBuffer = (char *) malloc (scheme.length () + 1);
          strcpy (lineBuffer, scheme.c_str ());
          parseScheme (lineBuffer, &(schemes->at (schemeId)));
          schemeId++;
          free (lineBuffer);
          cout << "DONE " << schemeId << endl;
        }
      }
      delete defSchemes;

      /** OUTPUT **/
      value = ini.Get (OUTPUT_TAG, OUTPUT_BASE_PATH, "").c_str ();
      if (strcmp (value, ""))
      {
        if (value[strlen (value) - 1] != char_separator)
          outputBasePath = string (value) + os_separator;
        else
          outputBasePath = string (value);
      }
      value = ini.Get (OUTPUT_TAG, OUTPUT_TMP_PATH, "").c_str ();
      if (strcmp (value, ""))
      {
        if (value[strlen (value) - 1] != char_separator)
          outputTmpPath = string (value) + os_separator;
        else
          outputTmpPath = string (value);
      }
      else
      {
        outputTmpPath = outputBasePath;
      }
      value = ini.Get (OUTPUT_TAG, OUTPUT_RESULTS_TAG, "").c_str ();
      if (strcmp (value, ""))
      {
        outputFileResults = outputBasePath + string (value);
      }
      value = ini.Get (OUTPUT_TAG, OUTPUT_MODELS_TAG, "").c_str ();
      if (strcmp (value, ""))
      {
        outputFileModels = outputBasePath + string (value);
      }
      value = ini.Get (OUTPUT_TAG, OUTPUT_PARTS_TAG, "").c_str ();
      if (strcmp (value, ""))
      {
        outputFilePartitions = outputBasePath + string (value);
      }
      value = ini.Get (OUTPUT_TAG, OUTPUT_SCHEMES_TAG, "").c_str ();
      if (strcmp (value, ""))
      {
        outputFileSchemes = outputBasePath + string (value);
      }

      if ((optimize_mode == OPT_CUSTOM) & (data_type == DT_NUCLEIC))
      {
        cerr << "[ERROR] Custom model set is not available for nucleic data."
            << endl;
        cerr << "        Please use \"gtr\" or \"all\" option." << endl;
        exit_partest (EX_CONFIG);
      }
      else if ((optimize_mode == OPT_GTR) & (data_type == DT_PROTEIC))
      {
        cerr << "[ERROR] GTR option is not available for proteic data." << endl;
        cerr << "        Please use \"all\" option or define the set of models."
            << endl;
        exit_partest (EX_CONFIG);
      }

      if (strcmp (inputFile, ""))
      {
        input_file = new string (inputFile);
      }
      if (strcmp (userTree, ""))
      {
        user_tree = new string (userTree);
      }
    }
  }

  int ConfigParser::parsePartitionDetails (char * line,
                                           struct partitionInfo * pInfo)
  {

    pllQueue * partitionLine;
    pllPartitionInfo * pi;
    char * the_line = (char *) malloc (strlen (line) + 10);
    if (data_type == DT_NUCLEIC)
    {
      strcpy (the_line, "DNA, P=");
    }
    else
    {
      strcpy (the_line, "AUTO, P=");
    }
    strcat (the_line, line);
    the_line[strlen (the_line)] = '\0';
    partitionLine = pllPartitionParseString (the_line);
    if (!partitionLine)
    {
      cerr << "[ERROR] Could not read partition: \"" << line << "\"" << endl;
      exit_partest (EX_IOERR);
    }
    pi = (pllPartitionInfo *) partitionLine->head->item;
    int numberOfSections = 0;
    for (pllQueueItem * qItem = pi->regionList->head; qItem; qItem =
        qItem->next)
    {
      pllPartitionRegion * region = (pllPartitionRegion *) qItem->item;
      pInfo->start[numberOfSections] = region->start;
      pInfo->end[numberOfSections] = region->end;
      if ((region->end < region->start) || region->stride <= 0)
      {
        cerr << "[ERROR] Invalid partition: \"" << line << "\"" << endl;
        exit_partest (EX_IOERR);
      }
      pInfo->stride[numberOfSections] = region->stride;
      numberOfSections++;
    }
    pInfo->numberOfSections = numberOfSections;
    free (the_line);

    return 0;
  }

  int ConfigParser::parseScheme (char * line, t_partitioningScheme * scheme)
  {

    char * parsed = strtok (line, "(");
    char * rest = strtok (NULL, "\0");
    parsed = strtok (parsed, ")");
    while (parsed != NULL)
    {
      char * parsedPart;
      parsedPart = strtok (parsed, ",");
      t_partitionElementId nextPart;
      while (parsedPart != NULL)
      {
        Utilities::toLower (parsedPart);
        size_t nextSingleElement = number_of_genes;
        for (size_t i = 0; i < partitions->size (); i++)
        {
          partitionInfo pInfo = partitions->at (i);
          if (!pInfo.name.compare (parsedPart))
          {
            nextSingleElement = pInfo.partitionId.at (0);
            break;
          }
        }
        if (nextSingleElement >= number_of_genes)
        {
          cerr << "[ERROR] Partition " << parsedPart << " not found" << endl;
          exit_partest (EX_IOERR);
        }
        nextPart.push_back (nextSingleElement);
        parsedPart = strtok (NULL, ",");
      }
      scheme->push_back (nextPart);

      parsed = strtok (rest, "(");
      rest = strtok (NULL, "\0");
      parsed = strtok (parsed, ")");
    }

    return 0;
  }

  ConfigParser::~ConfigParser ()
  {
    if (partitions)
    {
      delete partitions;
    }
  }

  void ConfigParser::printFormat ()
  {
    cout << "Config file format:" << endl << endl;
    cout << "   ;THIS IS A COMMENT" << endl;
    cout << endl << "   ; Start of input data" << endl;
    cout << "   [" << INPUT_TAG << "]" << endl;
    cout << "   " << INPUT_MSA_TAG << " = INPUT_ALIGNMENT_FILE" << endl;
    cout << "   " << INPUT_TREE_TAG << " = INPUT_TREE_FILE" << endl;
    cout << "   " << INPUT_DATATYPE_TAG << " = {nt|aa} (default: nt)" << endl;
    cout << "   " << INPUT_KEEPBRANCHES_TAG
        << " = {true|false} (default: false)" << endl;
    cout << endl << "   ; Start of searching options" << endl;
    cout << "   [" << SEARCH_TAG << "]" << endl;
    cout << "   " << SEARCH_ALGORITHM_TAG
        << "={greedy|hcluster|random|auto} (default: auto)";
    cout << "   " << SEARCH_ALGORITHM_REPS << " = # (default: 1)" << endl;
    cout << endl << "   ; Start of candidate models description" << endl;
    cout << "   [" << MODELS_TAG << "]" << endl;
    cout << "   " << MODELS_INCLUDE_TAG
        << " = {all | gtr | [LIST]} (default:all)" << endl;
    cout << "   " << MODELS_EPSILON_TAG << " = EPSILON_THRESHOLD (default:auto)"
        << endl;
    cout << "   " << MODELS_DO_M_TAG << " = {true|false} (default:true)"
        << endl;
    cout << "   " << MODELS_DO_F_TAG << " = {true|false} (default:false)"
        << endl;
    cout << "   " << MODELS_DO_G_TAG << " = {true|false} (default:false)"
        << endl;
#ifdef _IG_MODELS
    cout << "   " << MODELS_DO_I_TAG << " = {true|false} (default:false)" << endl;
    cout << "   " << MODELS_DO_IG_TAG << " = {true|false} (default:false)" << endl;
#endif
    cout << endl << "   ; Start of partitions" << endl;
    cout << "   [" << PARTITIONS_TAG << "]" << endl;
    cout << "   PART1 = INI1-END1" << endl;
    cout << "   PART1 = INI2-END2" << endl;
    cout << "   ..." << endl;
    cout << "   PARTn = INIn-ENDn" << endl;
    cout << endl << "   ; Start of schemes" << endl;
    cout << "   [" << SCHEMES_TAG << "]" << endl;
    cout << "   S1 = (GENE1,GENE2)(GENE3)..." << endl;
    cout << "   S2 = (GENE1,GENE2,GENE3)..." << endl;
    cout << "   ..." << endl;
    cout << "   Sn = (GENE1)(GENE2,GENE3,...)" << endl;
    cout << endl << "   ; Start of output section" << endl;
    cout << "   [" << OUTPUT_TAG << "]" << endl;
    cout << "   " << OUTPUT_BASE_PATH
        << " = OUTPUT_BASE_URL (default:partitiontest_FILENAME" << endl;
    cout << "   " << OUTPUT_MODELS_TAG
        << " = OUTPUT_MODELS_FILE (default: $path/models)" << endl;
    cout << "   " << OUTPUT_SCHEMES_TAG
        << " = OUTPUT_SCHEMES_FILE (default: $path/schemes)" << endl;
    cout << "   " << OUTPUT_RESULTS_TAG
        << " = OUTPUT_RESULTS_FILE (default: $path/results)" << endl;
    cout << endl << "Example:" << endl << endl;
    cout << "   [" << PARTITIONS_TAG << "]" << endl;
    cout << "   DNA1 = 1-976" << endl;
    cout << "   DNA2 = 976-1803" << endl;
    cout << "   DNA3 = 1804-2700" << endl;
    cout << "   DNA4 = 2701-3815" << endl << endl;
    cout << INPUT_TAG << "/" << INPUT_DATATYPE_TAG << endl;
    cout << "      nt - Nucleic (DNA) data" << endl;
    cout << "      aa - Amino Acid (protein) data" << endl;
    cout << SEARCH_TAG << "/" << SEARCH_ALGORITHM_TAG << endl;
    cout << "      greedy   - Greedy search" << endl;
    cout << "      hcluster - Hierarchical Clustering search" << endl;
    cout << "      random   - Random search" << endl;
    cout << "      auto     - Auto selected search algorithm" << endl;
    cout << SEARCH_TAG << "/" << SEARCH_ALGORITHM_REPS << endl;
    cout
        << "      (int) # - Maximum number of replicates on each HCluster/Random step"
        << endl;
    cout << MODELS_TAG << "/" << MODELS_INCLUDE_TAG << endl;
    cout << "      all    - Evaluate the whole set of models" << endl;
    cout << "      gtr    - Evaluate only gtr models (Only for DNA data)"
        << endl;
    cout
        << "      [LIST] - List of matrices to evaluate (Only for protein data)"
        << endl;
    cout
        << "               { dayhoff, dcmut, jtt, mtrev, cprev, rtrev, wag, vt,"
        << endl;
    cout
        << "                 blossum62, mtmam, mtart, hivb, hivw, mtzoa, pmb, flu }"
        << endl;
    cout << "               e.g.,  " << MODELS_INCLUDE_TAG
        << "=dayhoff dcmut jtt" << endl;
    cout << SCHEMES_TAG << endl;
    cout << "      Define a set of schemes to evaluate instead of searching."
        << endl;
    cout << "               e.g., "
        << "S1=(GENE1,GENE2)(GENE3,GENE5,GENE6)(GENE4)" << endl;
    cout << OUTPUT_TAG << endl;
    cout
        << "      Define urls for the output files. set to N/A for avoid output"
        << endl;
    cout << "               e.g., " << OUTPUT_MODELS_TAG << "=N/A" << endl;
  }

  void ConfigParser::createTemplate ()
  {
    cout << ";THIS IS A COMMENT" << endl;
    cout << "; Start of input data" << endl;
    cout << "[" << INPUT_TAG << "]" << endl;
    cout << INPUT_MSA_TAG << " = input.phy" << endl;
    cout << INPUT_TREE_TAG << " = input.tree" << endl;
    cout << INPUT_DATATYPE_TAG << " = nt" << endl;
    cout << INPUT_KEEPBRANCHES_TAG << " = false" << endl << endl;
    cout << "; Start of candidate models description" << endl;
    cout << "[" << MODELS_TAG << "]" << endl;
    cout << MODELS_INCLUDE_TAG << " = all" << endl << endl;
    cout << MODELS_EPSILON_TAG << " = 1.0" << endl;
    cout << MODELS_DO_M_TAG << " = true" << endl;
    cout << MODELS_DO_F_TAG << " = true" << endl;
    cout << MODELS_DO_G_TAG << " = true" << endl;
#ifdef _IG_MODELS
    cout << MODELS_DO_I_TAG << " = false" << endl;
    cout << MODELS_DO_IG_TAG << " = false" << endl;
#endif
    cout << "; Start of output data" << endl;
    cout << "[" << OUTPUT_TAG << "]" << endl;
    cout << OUTPUT_BASE_PATH << " = /tmp" << endl;
    cout << "; Start of searching options" << endl;
    cout << "[" << SEARCH_TAG << "]" << endl;
    cout << SEARCH_ALGORITHM_TAG << " = auto" << endl;
    cout << SEARCH_ALGORITHM_REPS << " = 1" << endl;
    cout << "[" << PARTITIONS_TAG << "]" << endl;
    cout << "PART1 = INI1-END1" << endl;
    cout << "PART2 = INI2-END2" << endl;
    cout << "..." << endl;
    cout << "PARTn = INIn-ENDn" << endl;
  }

  vector<partitionInfo> * ConfigParser::getPartitions ()
  {
    return partitions;
  }

  struct partitionInfo ConfigParser::getPartition (size_t index)
  {
    if (!partitions)
    {
      cerr << "[ERROR] No partitions were defined" << endl;
      assert(0);
    }
    if (index >= number_of_genes)
    {
      cerr << "[ERROR] Requested partition does not exist" << endl;
      assert(0);
    }
    return partitions->at (index);
  }

} /* namespace partest */
