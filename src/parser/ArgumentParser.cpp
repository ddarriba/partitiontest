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
 * @file ArgumentParser.cpp
 */

#include "ArgumentParser.h"
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include "ConfigParser.h"
#include "util/Utilities.h"
#include "util/FileUtilities.h"

using namespace std;

namespace partest
{

#ifdef _IG_MODELS
#define NUM_ARGUMENTS 28
#else
#define NUM_ARGUMENTS 26
#endif

  void ArgumentParser::init ()
  {
    index = 1;
    subindex = 1;
  }

  ArgumentParser::ArgumentParser (PartitionTest * _ptest) :
      index (0), subindex (0), ptest (_ptest)
  {

    option options_list[] =
      {
        { ARG_HELP, 'h', "help", false },
        { ARG_CONFIG_FILE, 'c', "config-file", true },
        { ARG_CONFIG_HELP, 0, "config-help", false },
        { ARG_CONFIG_TEMPLATE, 0, "config-template", false },
        { ARG_DATA_TYPE, 'd', "data-type", true },
        { ARG_DISABLE_CHECKPOINT, 0, "disable-ckp", false },
        { ARG_DISABLE_OUTPUT, 0, "disable-output", false },
        { ARG_EPSILON, 'e', "epsilon", true },
        { ARG_INPUT_FORMAT, 'f', "input-format", true },
        { ARG_FORCE_OVERRIDE, 0, "force-override", false },
        { ARG_FREQUENCIES, 'F', "empirical-frequencies", false },
        { ARG_PERGENE_BL, 'g', "pergene-bl", false },
#ifdef _IG_MODELS
        { ARG_GAMMA, 'G', "gamma-rates", false},
        { ARG_INV, 'I', "invariant-sites", false},
#endif
        { ARG_INPUT_FILE, 'i', "input-file", true },
        { ARG_KEEP_BRANCH_LENGTHS, 'k', "keep-branches", false },
        { ARG_SAMPLE_SIZE, 'n', "sample-size", true },
        { ARG_NON_STOP, 'N', "non-stop", false },
        { ARG_OUTPUT, 'o', "output", true },
        { ARG_OPTIMIZE, 'O', "optimize", true },
        { ARG_NUM_PROCS, 'p', "num-procs", true },
        { ARG_HCLUSTER_REPS, 'r', "replicates", true },
        { ARG_IC_TYPE, 's', "selection-criterion", true },
        { ARG_SEARCH_ALGORITHM, 'S', "search", true },
        { ARG_TOPOLOGY, 't', "topology", true },
        { ARG_FINAL_TREE, 'T', "get-final-tree", false },
        { ARG_USER_TREE, 'u', "user-tree", true },
        { ARG_VERBOSE, 'v', "verbose", true },
        { ARG_VERSION, 'V', "version", false }
      };

    size_t size = NUM_ARGUMENTS * sizeof(option);

    arguments = (option *) malloc (size);
    memcpy (arguments, options_list, size);
    init ();
  }

  ArgumentParser::~ArgumentParser ()
  {
    free (arguments);
  }

  ArgIndex ArgumentParser::get_opt (int argc, char *argv[], char *argument,
                                    char *value)
  {
    if (index >= argc)
      return ARG_END;
    else
    {
      ArgIndex arg_index = ARG_NULL;
      char *arg = argv[index];
      strcpy (argument, argv[index]);
      if (*arg == '-')
      {
        if (*(arg + 1) == '-')
        {
          /* check long option */
          for (int j = 0; j < NUM_ARGUMENTS; j++)
          {
            if (arguments[j].long_code > 0)
            {
              if (!strcmp (arguments[j].long_code, &(arg[2])))
              {
                arg_index = arguments[j].index;
                if (arguments[j].required_value)
                {
                  index++;
                  if (index >= argc || argv[index][0] == '-')
                  {
                    cerr << "ERROR! Argument " << arg << " requires a value."
                        << endl;
                    exit_partest (EX_CONFIG);
                  }
                  strcpy (value, argv[index]);
                }
              }
            }
          }
          index++;
        }
        else
        {
          /* check short option */
          int num_arguments = (int) strlen (arg + 1);
          for (int j = 0; j < NUM_ARGUMENTS; j++)
          {
            if (arguments[j].char_code == arg[subindex])
            {
              arg_index = arguments[j].index;
              if (arguments[j].required_value)
              {
                if (num_arguments > 1)
                {
                  cerr << "ERROR! Argument " << arg[subindex]
                      << " cannot be used in a row." << endl;
                  exit_partest (EX_CONFIG);
                }
                else
                {
                  index++;
                  if (index >= argc || argv[index][0] == '-')
                  {
                    cerr << "ERROR! Argument " << arg[subindex]
                        << " requires a value." << endl;
                    exit_partest (EX_CONFIG);
                  }
                  strcpy (value, argv[index]);
                }
              }
            }
          }
          if (num_arguments > subindex)
          {
            subindex++;
          }
          else
          {
            subindex = 1;
            index++;
          }
          if (arg_index == ARG_NULL)
          {
            cerr << "[ERROR] \"" << argument << "\" option not recognized."
                << endl;
            exit_partest (EX_CONFIG);
          }
        }
      }
      else
      {
        /* not an option! */
        cerr << "[ERROR] \"" << argv[index] << "\" is not a valid option."
            << endl;
        cerr << "Arguments must start with '-' (short) or '--' (long)" << endl;
        exit_partest (EX_CONFIG);
      }

      return arg_index;
    }
  }

  bool ArgumentParser::parseConfigFile (int argc, char *argv[])
  {
    int argument_index;
    char value[256];
    char argument[256];
    bool configSet = false;

    init ();
    while ((argument_index = get_opt (argc, argv, argument, value)) != ARG_END)
    {
      if (argument_index == ARG_CONFIG_FILE)
      {
        /* input configuration file */
        if (!FileUtilities::existsFile (value))
        {
          cerr << "[ERROR] \"Configuration file " << value
              << "\" does not exist." << endl;
          exit_partest (EX_IOERR);
        }
        configSet = true;
        ptest->setConfigFile (value);
        break;
      }
    }
    return configSet;
  }

  void ArgumentParser::parse (int argc, char *argv[])
  {

    int argument_index;
    char value[256];
    char argument[256];
    char _input_file[256] = "";
    char _user_tree[256] = "";
    char _config_file[256] = "";
    char _output_dir[256] = "";

#ifdef _SELECT_SAMPLE_SIZE
    double sampleSizeValue = 0.0;
    SampleSize sampleSize = SS_DEFAULT;
#endif

    init ();
    while ((argument_index = get_opt (argc, argv, argument, value)) != ARG_END)
    {
      switch (argument_index)
        {
        case ARG_HELP:
          /* display usage */
          exit_partest (EX_USAGE);
          break;
        case ARG_INPUT_FILE:
          /* input alignment file */
          strcpy (_input_file, value);
          if (!FileUtilities::existsFile (_input_file))
          {
            cerr << "[ERROR] \"Input file " << value << "\" does not exist."
                << endl;
            exit_partest (EX_IOERR);
          }
          break;
        case ARG_USER_TREE:
          /* input user defined topology */
          strcpy (_user_tree, value);
          if (!FileUtilities::existsFile (_user_tree))
          {
            cerr << "[ERROR] \"Tree file " << value << "\" does not exist."
                << endl;
            exit_partest (EX_IOERR);
          }
          break;
        case ARG_DATA_TYPE:
          /* data type (nucleotide or amino-acid) */
          if (!strcmp (value, ARG_DT_PROTEIC))
          {
            ptest->setDataType (DT_PROTEIC);
          }
          else if (!strcmp (value, ARG_DT_NUCLEIC))
          {
            ptest->setDataType (DT_NUCLEIC);
          }
          else
          {
            cerr << "[ERROR] \"-d " << value
                << "\" is not a valid data type. Use one of the following:"
                << endl;
            cerr << "  -d " << setw (8) << left << ARG_DT_PROTEIC
                << "Protein sequence alignment" << endl;
            cerr << "  -d " << setw (8) << left << ARG_DT_NUCLEIC
                << "DNA sequence alignment (DEFAULT)" << endl;
            exit_partest (EX_CONFIG);
          }
          break;
        case ARG_DISABLE_CHECKPOINT:
          /* disable checkpointing files */
          ckpAvailable = false;
          break;
        case ARG_DISABLE_OUTPUT:
          /* disable writing output files */
          outputAvailable = false;
          ckpAvailable = false;
          break;
        case ARG_EPSILON:
          /* epsilon used for optimization algorithm */
          if (Utilities::isNumeric (value))
          {
            epsilon = atof (value);
            if (epsilon <= 0.0f)
            {
              cerr << "[ERROR] \"-e " << value
                  << "\" is not a valid value. Epsilon should be \"auto\","
                  << " or a numeric value greater or equal than 0."
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
              cerr << "[ERROR] \"-e " << value
                  << "\" is not a valid value. Epsilon should be \"auto\","
                  << " or a numeric value greater or equal than 0."
                  << endl;
              exit_partest (EX_CONFIG);
            }
          }
          break;
        case ARG_FORCE_OVERRIDE:
          /* disable writing output files */
          force_overriding = true;
          break;
        case ARG_TOPOLOGY:
          /* starting topology (Fixed, Parsimony or User-defined) */
          if (!strcmp (value, ARG_TOPO_MP))
          {
            ptest->setStartingTopology (StartTopoMP);
          }
          else if (!strcmp (value, ARG_TOPO_ML))
          {
            ptest->setStartingTopology (StartTopoML);
          }
          else if (!strcmp (value, ARG_TOPO_FIXED))
          {
            ptest->setStartingTopology (StartTopoFIXED);
          }
          else if (!strcmp (value, ARG_TOPO_USER))
          {
            ptest->setStartingTopology (StartTopoUSER);
          }
          else
          {
            cerr << "[ERROR] \"-t " << value
                << "\" is not a valid input topology. Use one of the following:"
                << endl;
            cerr << "  -t " << setw (8) << left << ARG_TOPO_MP
                << "Maximum Parsimony topology" << endl;
            cerr << "  -t " << setw (8) << left << ARG_TOPO_ML
                << "Maximum Likelihood topology" << endl;
            cerr << "  -t " << setw (8) << left << ARG_TOPO_FIXED
                << "Fixed Maximum Likelihood topology for every model" << endl;
            cerr << "  -t " << setw (8) << left << ARG_TOPO_USER
                << "User-defined topology" << endl;
            exit_partest (EX_CONFIG);
          }
          break;
        case ARG_SEARCH_ALGORITHM:
          /* search algorithm (HCluster, Greedy, Random or Exhaustive) */
          if (!strcmp (value, ARG_SEARCH_K1))
          {
            ptest->setSearchAlgo (SearchK1);
          }
          else if (!strcmp (value, ARG_SEARCH_KN))
          {
            ptest->setSearchAlgo (SearchKN);
          }
          else if (!strcmp (value, ARG_SEARCH_EXHAUSTIVE))
          {
            ptest->setSearchAlgo (SearchExhaustive);
          }
          else if (!strcmp (value, ARG_SEARCH_RANDOM))
          {
            ptest->setSearchAlgo (SearchRandom);
          }
          else if (!strcmp (value, ARG_SEARCH_GREEDY))
          {
            ptest->setSearchAlgo (SearchGreedy);
          }
          else if (!strcmp (value, ARG_SEARCH_GREEDY_EXT))
          {
            ptest->setSearchAlgo (SearchGreedyExtended);
          }
          else if (!strcmp (value, ARG_SEARCH_HIERARCHICAL))
          {
            ptest->setSearchAlgo (SearchHCluster);
          }
          else if (!strcmp (value, ARG_SEARCH_AUTO))
          {
            ptest->setSearchAlgo (SearchAuto);
          }
          else
          {
            cerr << "[ERROR] \"-S " << value
                << "\" is not a valid search algorithm."
                << " Use one of the following:"
                << endl;
            cerr << "  -S " << setw (12) << left << ARG_SEARCH_KN
                << "Evaluate only the user defined partitioning scheme" << endl;
            cerr << "  -S " << setw (12) << left << ARG_SEARCH_EXHAUSTIVE
                << "Exhaustive algorithm (horribly computationally expensive)"
                << endl;
            cerr << "  -S " << setw (12) << left << ARG_SEARCH_RANDOM
                << "Random walk algorithm (Chinese restaurant process)" << endl;
            cerr << "  -S " << setw (12) << left << ARG_SEARCH_GREEDY
                << "Greedy hill climbing algorithm" << endl;
            cerr << "  -S " << setw (12) << left << ARG_SEARCH_GREEDY_EXT
                << "Extended greedy hill climbing algorithm" << endl;
            cerr << "  -S " << setw (12) << left << ARG_SEARCH_HIERARCHICAL
                << "Hierarchical clustering algorithm" << endl;
            exit_partest (EX_CONFIG);
          }
          break;
        case ARG_HCLUSTER_REPS:
          /* number of replicates for HCluster and Random */
          if (Utilities::isInteger (value))
          {
            /* apply absolute number of replicates */
            ptest->setMaxSamples (atoi (value));
            ptest->setSamplesPercent (0.0);
          }
          else
          {
            /* apply percentage of replicates */
            ptest->setMaxSamples (0);
            if (value[strlen (value) - 1] == '%')
            {
              value[strlen (value) - 1] = '\0';
              if (Utilities::isInteger (value))
              {
                ptest->setSamplesPercent (atof (value) / 100);
              }
              else
              {
                cerr << "[ERROR] Invalid samples percent " << samples_percent
                    << endl;
                exit_partest (EX_CONFIG);
              }
            }
            else
            {
              if (Utilities::isNumeric (value))
              {
                ptest->setSamplesPercent (atof (value));
              }
              else
              {
                cerr << "[ERROR] Invalid samples percent " << samples_percent
                    << endl;
                exit_partest (EX_CONFIG);
              }
            }
          }
          break;
        case ARG_KEEP_BRANCH_LENGTHS:
          /* keep branch lengths from the initial topology */
          reoptimize_branch_lengths = false;
          break;
        case ARG_NON_STOP:
          /* continue until the end / non stop on local maxima */
          non_stop = true;
          break;
        case ARG_IC_TYPE:
          /* information criterion (aic, bic, aicc, dt) */
          if (!strcmp (value, ARG_IC_AIC))
          {
            ptest->setIcType (AIC);
          }
          else if (!strcmp (value, ARG_IC_BIC))
          {
            ptest->setIcType (BIC);
          }
          else if (!strcmp (value, ARG_IC_AICC))
          {
            ptest->setIcType (AICC);
          }
          else if (!strcmp (value, ARG_IC_DT))
          {
            ptest->setIcType (DT);
          }
          else
          {
            cerr << "[ERROR] \"-s " << value
                << "\" is not a valid criterion. Use one of the following:"
                << endl;
            cerr << "  -s " << setw (8) << left << ARG_IC_AIC
                << "Akaike Information Criterion" << endl;
            cerr << "  -s " << setw (8) << left << ARG_IC_BIC
                << "Bayesian Information Criterion" << endl;
            cerr << "  -s " << setw (8) << left << ARG_IC_AICC
                << "Corrected Akaike Information Criterion" << endl;
            cerr << "  -s " << setw (8) << left << ARG_IC_DT
                << "Decision Theory" << endl;
            exit_partest (EX_CONFIG);
          }
          break;
#ifdef _SELECT_SAMPLE_SIZE
          case ARG_SAMPLE_SIZE:
          /* sample size used for information criteria */
          sampleSizeValue = 0.0;
          if (!strcmp(value, ARG_SS_ALIGN))
          {
            sampleSize = SS_ALIGNMENT;
          }
          else if (!strcmp(value, ARG_SS_SHANNON))
          {
            sampleSize = SS_SHANNON;
          }
          else
          {
            sampleSizeValue = atof(value);
            if (!sampleSizeValue)
            {
              cerr << "[ERROR] \"-n " << value
              << "\" is not a valid sample size. Use one of the following:"
              << endl;
              cerr << "  -n " << setw(16) << left << ARG_SS_ALIGN
              << "Alignment size" << endl;
              cerr << "  -n " << setw(16) << left << ARG_SS_SHANNON
              << "Shannon entropy of the alignment" << endl;
              cerr << "  -n " << setw(16) << left << "[CUSTOM_VALUE]"
              << "Custom sample size" << endl;
              exit_partest(EX_CONFIG);
            }
            sampleSize = SS_CUSTOM;
          }
          break;
#endif
        case ARG_PERGENE_BL:
          pergene_branch_lengths = true;
          break;
        case ARG_FREQUENCIES:
          /* include empirical / unequal frequencies */
          do_rate |= (RateVarF);
          break;
        case ARG_INV:
          /* include models with a proportion of invariant sites */
          do_rate |= (RateVarI);
          break;
        case ARG_GAMMA:
          /* include models with gamma-distributed rate heterogeneity */
          do_rate |= (RateVarG);
          break;
        case ARG_OPTIMIZE:
          /* modelset to optimize */
          if (!strcmp (value, ARG_OPTIMIZE_BESTMODEL))
          {
            ptest->setOptimize (OPT_SEARCH);
          }
          else if (!strcmp (value, ARG_OPTIMIZE_GTR))
          {
            ptest->setOptimize (OPT_GTR);
          }
          else
          {
            cerr << "[ERROR] \"-n " << value
                << "\" is not a valid optimize mode. Use one of the following:"
                << endl;
            cerr << "  -O " << setw (16) << left << ARG_OPTIMIZE_BESTMODEL
                << "\t Perform a model selection on the best partition" << endl;
            cerr << "  -O " << setw (16) << left << ARG_OPTIMIZE_GTR
                << "\t Optimize only GTR models on the best partition" << endl;
            exit_partest (EX_CONFIG);
          }
          break;
        case ARG_FINAL_TREE:
          /* compute final ML tree */
          compute_final_tree = true;
          break;
        case ARG_OUTPUT:
          /* output directory */
          strcpy (_output_dir, value);
          break;
        case ARG_CONFIG_FILE:
          /* input configuration file */
          strcpy (_config_file, value);
          break;
        case ARG_CONFIG_HELP:
          /* display configuration file format */
          ConfigParser::printFormat ();
          exit_partest (EX_OK);
          break;
        case ARG_CONFIG_TEMPLATE:
          /* display configuration file template */
          ConfigParser::createTemplate ();
          exit_partest (EX_OK);
          break;
        case ARG_NUM_PROCS:
          /* set the number of threads */
#ifdef HAVE_PTHREADS
          number_of_threads = atoi (value);
#else
          cerr
          << "[ERROR] PThreads version is not available."
          << " You must recompile with PTHREADS flag."
          << endl;
          exit_partest(EX_CONFIG);
#endif
          break;
        case ARG_VERBOSE:
          /* verbosity level */
          if (Utilities::isInteger (value))
          {
            /* apply absolute number of replicates */
            verbosity = atoi (value);
          }
          else
          {
            /* error */
            cerr << "[ERROR] \"-v " << value
                << "\" is not a valid verbosity level. "
                << "It must be 0, 1 or 2." << endl;
            exit_partest (EX_CONFIG);
          }
          break;
        case ARG_VERSION:
          /* display application version */
          cout << "PartitionTest v" << PROGRAM_VERSION << " - " << PROGRAM_DATE
              << endl;
          cout
              << "Copyright (C) 2014"
              << " D.Darriba, G.L.Taboada, R.Doallo and David Posada"
              << endl;
          cout
              << "License GPLv3+: GNU GPL version 3 or later"
              << " <http://gnu.org/licenses/gpl.html>."
              << endl;
          cout
              << "This is free software:"
              << " you are free to change and redistribute it."
              << endl;
          cout << "There is NO WARRANTY, to the extent permitted by law."
              << endl << endl;
          cout << "Written by Diego Darriba (ddarriba@udc.es)." << endl;
          exit_partest (EX_OK);
          break;
        default:
          cerr << "[ERROR] \"" << argument << "\" option not recognized."
              << endl;
          exit_partest (EX_CONFIG);
          break;
        }
    }
    if (strcmp (_input_file, ""))
      ptest->setInputFile (_input_file);
    if (strcmp (_user_tree, ""))
      ptest->setUserTree (_user_tree);
    if (strcmp (_output_dir, ""))
      ptest->setOutputDir (_output_dir);
  }

} /* namespace partest */
