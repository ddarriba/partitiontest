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
 * @file PartitionTest.cpp
 * @author Diego Darriba
 */

#include "PartitionTest.h"
#include "search/SearchAlgorithm.h"
#include "search/HierarchicalClusteringSearchAlgorithm.h"
#include "search/ExhaustiveSearchAlgorithm.h"
#include "search/GreedySearchAlgorithm.h"
#include "search/RandomSearchAlgorithm.h"
#include "indata/PartitioningScheme.h"
#include "indata/PartitionMap.h"
#include "util/PrintMeta.h"
#include "util/Utilities.h"
#include "util/FileUtilities.h"
#include "parser/ArgumentParser.h"
#include "parser/ConfigParser.h"

#include <fstream>
#include <cstring>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <cassert>

using namespace std;

namespace partest
{

  bool PartitionTest::checkParameters (void)
  {
    /* check all non-default parameters have been set */
    assert(input_file && input_file->compare (""));
    return true;
  }

  PartitionTest::PartitionTest ()
  {
    phylip = 0;
    pllPartitions = 0;
    tree = 0;
    pllPartsQueue = 0;

    data_type = DEFAULT_DATA_TYPE;
    do_rate = 0;
    starting_topology = DEFAULT_STARTING_TOPOLOGY;
    search_algo = DEFAULT_SEARCH_ALGO;
    max_samples = 1;
    optimize_mode = DEFAULT_OPTIMIZE;
    ic_type = DEFAULT_IC_TYPE;
  }

  PartitionTest::~PartitionTest ()
  {
    cout << timestamp () << " Execution done." << endl;
    if (phylip)
    {
      pllAlignmentDataDestroy (phylip);
    }
    if (pllPartsQueue)
      pllQueuePartitionsDestroy (&pllPartsQueue);
    if (tree)
    {
      if (pllPartitions && pllPartitions->alphaList)
      {
//			pllPartitionsDestroy(tree, &pllPartitions);
      }
      if (tree->nameHash)
      {
        pllDestroyInstance (tree);
        tree = 0;
      }
    }

  }

  bool PartitionTest::configure (PartitionTest * ptest, int argc, char * argv[])
  {

    ArgumentParser * argParser = new ArgumentParser (ptest);
    argParser->parseConfigFile (argc, argv);

    if (!config_file)
      config_file = new string ("");

    ConfigParser parser (config_file->c_str ());

    argParser->parse (argc, argv);
    parser.createPartitions ();

    delete argParser;

    // check required arguments
    if (!input_file)
    {
      if (!strlen (parser.getInputFile ()))
      {
        cerr << "ERROR! Input File (-i) is required!" << endl;
        exit_partest (EX_CONFIG);
      }
      input_file = new string (parser.getInputFile ());
    }

    if (outputAvailable)
    {
      if (parser.getOutputBasePath ().length () > 0)
      {
        output_dir = new string (parser.getOutputBasePath ());
      }
      else
      {
        const char *baseName = basename (input_file->c_str ());
        output_dir = new string ("partest_");
        output_dir->append (baseName);
        output_dir->append (os_separator);
      }
      models_logfile = new string (*output_dir + "models");
      schemes_logfile = new string (*output_dir + "schemes");
      results_logfile = new string (*output_dir + "results");
      log_logfile = new string (*output_dir + "log");

      ckpPath = (*output_dir) + CKP_DIR;
    }

    if (I_AM_ROOT)
    {
      if (outputAvailable)
      {
        int resMkdir = mkdir (output_dir->c_str (), 0777);
        int errorType = errno;
        if (resMkdir)
        {
          if (errorType == EEXIST)
          {
            if (force_overriding)
            {
              cerr << "[WARNING] Output directory " << (*output_dir)
                  << " already exists. Output files might be overwritten."
                  << endl;
              if (remove (log_logfile->c_str ()) && errno != ENOENT)
              {
                cerr << "[ERROR] There was an error removing log file "
                    << *log_logfile
                    << ". Please remove it manually and try again." << endl;
                exit_partest (EX_IOERR);
              }
            }
            else
            {
              if (FileUtilities::existsFile (*results_logfile))
              {
                cerr << "[ERROR] Output files already exist:" << endl;
                cerr << "        " << *results_logfile << endl;
                cerr << endl
                    << "If you really want to proceed, remove result files or re-run with --force-override or --disable-output arguments."
                    << endl << endl;
                exit_partest (EX_IOERR);
              }
            }
          }
          else
          {
            cerr << "[WARNING] ***** WARNING *****" << endl;
            cerr << "[WARNING] Output directory " << (*output_dir)
                << " cannot be created. No output files will be stored."
                << endl;
            cerr << "[WARNING] ***** WARNING *****" << endl;
          }
        }
        if (!resMkdir || errorType == EEXIST)
        {
          mkdir (ckpPath.c_str (), 0777);
          /* checkpointing is set to true unless user specified the oposite */
          ckpAvailable &= true;
        }
        else
        {
          ckpAvailable = false;
        }
      }
    }

#ifdef HAVE_MPI
    int tmpInt = ckpAvailable;
    MPI_Bcast (&tmpInt, 1, MPI_INT, 0, MPI_COMM_WORLD);
    ckpAvailable = tmpInt;
#endif

    pllInstanceAttr attr;
    attr.fastScaling = PLL_FALSE;
    attr.randomNumberSeed = 0x54321;
    attr.rateHetModel = PLL_GAMMA;
    attr.saveMemory = PLL_FALSE;
    attr.useRecom = PLL_FALSE;
    attr.numberOfThreads = number_of_threads;

    switch (optimize_mode)
      {
      case OPT_SEARCH:
        number_of_models = data_type == DT_NUCLEIC ?
        NUC_MATRIX_SIZE / 2 :
                                                     PROT_MATRIX_SIZE;
        if (do_rate & RateVarF)
        {
          number_of_models *= 2;
        }
        break;
      case OPT_GTR:
        number_of_models = 1;
        break;
      case OPT_CUSTOM:
        number_of_models = Utilities::setbitsCount (protModels);
        if (do_rate & RateVarF)
        {
          number_of_models *= 2;
        }
        break;
      default:
        assert(0);
      }

    tree = pllCreateInstance (&attr);

    phylip = pllParseAlignmentFile (PLL_FORMAT_PHYLIP, input_file->c_str ());
    if (!phylip)
    {
      cerr << "[ERROR] There was an error parsing input data \""
          << (*input_file) << "\"" << endl;
      exit_partest (EX_IOERR);
    }
    pllPartitions = pllPartitionsCommit (pllPartsQueue, phylip);
    if (!pllPartitions)
    {
      cerr << "[ERROR] There was an error parsing partitions data." << endl;
      exit_partest (EX_IOERR);
    }

    /* evaluate present states */
    bool frequenciesOK = true;
    double ** freqs = pllBaseFrequenciesAlignment (phylip, pllPartitions);
    for (int i = 0; i < pllPartitions->numberOfPartitions; i++)
    {
      for (int j = 0; j < pllPartitions->partitionData[i]->states; j++)
      {
        if (freqs[i][j] < 1e-12)
        {

          cerr << "[ERROR] State ";
          if (pllPartitions->partitionData[i]->dataType == PLL_DNA_DATA)
            cerr << dnaStateNames[j];
          else
            cerr << protStateNames[j];
          cerr << " is not present in partition "
              << pllPartitions->partitionData[i]->partitionName << endl;
          frequenciesOK = false;
        }
      }
    }
    if (!frequenciesOK)
    {
      cerr << "[ERROR] There are missing states. Please fix your data" << endl;
      exit_partest (EX_IOERR);
    }
    free (freqs);

    num_taxa = (size_t) phylip->sequenceCount;
    seq_len = (size_t) phylip->sequenceLength;

    num_patterns = (size_t) phylip->sequenceLength;

    return EX_OK;
  }

  void PartitionTest::setDataType (DataType dataType)
  {
    data_type = dataType;
  }

  void PartitionTest::setDoRate (bitMask doRate)
  {
    do_rate = doRate;
  }

  void PartitionTest::setIcType (InformationCriterion icType)
  {
    ic_type = icType;
  }

  void PartitionTest::setMaxSamples (int maxSamples)
  {
    if (maxSamples < 0)
    {
      cerr << "[ERROR] Max samples must be positive: (" << maxSamples << ")"
          << endl;
      exit_partest (EX_IOERR);
    }
    max_samples = maxSamples;
  }

  void PartitionTest::setSamplesPercent (double samplesPercent)
  {
    if (samplesPercent < 0 || samplesPercent > 1)
    {
      cerr << "[ERROR] Samples percent must be between 0 and 1: ("
          << samplesPercent << ")" << endl;
      exit_partest (EX_IOERR);
    }
    samples_percent = samplesPercent;
  }

  void PartitionTest::setOptimize (OptimizeMode optimize)
  {
    optimize_mode = optimize;
  }

  void PartitionTest::setSearchAlgo (SearchAlgo searchAlgo)
  {
    search_algo = searchAlgo;
  }

  void PartitionTest::setStartingTopology (StartTopo startingTopology)
  {
    starting_topology = startingTopology;
  }

  void PartitionTest::setConfigFile (const char * configFile)
  {
    if (config_file)
    {
      cerr << "[ERROR] Config file was already set before." << endl;
      exit_partest (EX_IOERR);
    }
    config_file = new string (configFile);
  }

  void PartitionTest::setInputFile (const char * inputFile)
  {
    if (input_file)
    {
      cerr << "[ERROR] Input file was already set before. " << input_file << " "
          << inputFile << endl;
      exit_partest (EX_IOERR);
    }
    input_file = new string (inputFile);
  }

  void PartitionTest::setOutputDir (const char * outputDir)
  {
    if (output_dir)
    {
      cerr << "[ERROR] Output directory was already set before." << endl;
      exit_partest (EX_IOERR);
    }
    output_dir = new string (outputDir);
  }

  void PartitionTest::setUserTree (const char * userTree)
  {
    if (user_tree)
    {
      cerr << "[ERROR] User tree was already set before." << endl;
      exit_partest (EX_IOERR);
    }
    user_tree = new string (userTree);
  }

}

using namespace partest;
using namespace std;

int main (int argc, char * argv[])
{

#ifdef HAVE_MPI
  if (MPI_Init (&argc, &argv))
  {
    cerr << "Error initializing MPI!!" << endl;
    exit_partest (EX_PROTOCOL);
  }

  MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
  MPI_Comm_rank (MPI_COMM_WORLD, &myRank);
#endif

  PartitionTest * ptest = new PartitionTest ();

  if (I_AM_ROOT)
  {
    PrintMeta::print_header (cout);
  }

  ptest->configure (ptest, argc, argv);

  if (!config_file && !ptest->checkParameters ())
  {
    cerr << "[ERROR] Configuration parameters are missing" << endl;
    exit_partest (EX_CONFIG);
  }

  if (I_AM_ROOT)
  {
    PrintMeta::print_options (cout);
  }

  /* In case user sets a fixed subset of schemes to optimize,
   * we perform an exhaustive optimization */
  if (number_of_schemes > 0)
    search_algo = SearchExhaustive;

  /* Instantiate the search algorithm */
  SearchAlgorithm * searchAlgo = 0;
  PartitioningScheme * bestScheme;
  if (search_algo != SearchAuto)
  {
    switch (search_algo)
      {
      case SearchK1:
        {
          number_of_schemes = 1;
          if (schemes)
            delete schemes;
          schemes = new vector<t_partitioningScheme> (number_of_schemes);
          /* building K=1 scheme */
          t_partitionElementId geneId (number_of_genes);
          for (size_t gene = 0; gene < number_of_genes; gene++)
          {
            geneId.at (gene) = gene;
          }
          schemes->at (0).push_back (geneId);
          searchAlgo = new ExhaustiveSearchAlgorithm ();
          break;
        }
      case SearchKN:
        {
          number_of_schemes = 1;
          if (schemes)
            delete schemes;
          schemes = new vector<t_partitioningScheme> (number_of_schemes);
          /* building K=N scheme */
          for (size_t gene = 0; gene < number_of_genes; gene++)
          {
            t_partitionElementId geneId (1);
            geneId.at (0) = gene;
            schemes->at (0).push_back (geneId);
          }
          searchAlgo = new ExhaustiveSearchAlgorithm ();
          break;
        }
      case SearchHCluster:
        searchAlgo = new HierarchicalClusteringSearchAlgorithm ();
        break;
      case SearchGreedy:
        searchAlgo = new GreedySearchAlgorithm ();
        break;
      case SearchRandom:
        searchAlgo = new RandomSearchAlgorithm ();
        break;
      case SearchExhaustive:
        searchAlgo = new ExhaustiveSearchAlgorithm ();
        break;
      case SearchGreedyExtended:
        cerr << "[ERROR] Not implemented yet" << endl;
        exit_partest (EX_UNAVAILABLE);
        break;
      default:
        break;
      }

    bestScheme = searchAlgo->start ();
  }
  else
  {
    /* Auto search algorithm selection according to the number of partitions */
    if (number_of_genes <= 20)
    {
      searchAlgo = new GreedySearchAlgorithm ();
      cout << "Searching with greedy algorithm" << endl;
      bestScheme = searchAlgo->start ();
    }
    else
    {
      searchAlgo = new HierarchicalClusteringSearchAlgorithm ();
      cout << "Searching with hierarchical clustering" << endl;
      bestScheme = searchAlgo->start ();
      max_samples = (int) bestScheme->getNumberOfElements ();
      cout << "Searching with hierarchical clustering with " << max_samples
          << " samples" << endl;
      bestScheme = searchAlgo->start (bestScheme);
      if (bestScheme->getNumberOfElements () <= 20)
      {
        delete searchAlgo;
        searchAlgo = new GreedySearchAlgorithm ();
        cout << "Searching with greedy algorithm" << endl;
        bestScheme = searchAlgo->start (bestScheme);
      }
    }
  }

  if (I_AM_ROOT)
  {
    PrintMeta::print_results (cout, bestScheme);

    if (compute_final_tree)
    {
      ModelOptimize mo;
      mo.buildFinalTree (bestScheme, true);
    }

    if (outputAvailable && results_logfile)
    {
      ofstream ofs (results_logfile->c_str (), ios::out);
      PrintMeta::print_results_xml (ofs, bestScheme);
      ofs.close ();
    }
  }

  delete bestScheme;
  delete searchAlgo;
  if (number_of_schemes > 0)
    delete schemes;

  PartitionMap::deleteInstance ();

  delete ptest;

#ifdef HAVE_MPI
  MPI_Finalize ();
#endif

#ifdef CLEAN_ON_EXIT
  /* cleanup */
  if (ckpAvailable)
  {
    FileUtilities::deleteDirectoryRecursive(ckpPath.c_str());
  }
#endif

  exit_partest (EX_OK);
}
