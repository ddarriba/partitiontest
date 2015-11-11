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

/*
 * @file PrintMeta.cpp
 * @author Diego Darriba
 */

#include "PrintMeta.h"
#include "parser/ArgumentParser.h"
#include "util/Utilities.h"
#include "util/GlobalDefs.h"
#include <pll/parsePartition.h>

#include <iomanip>
#include <string.h>
#include <sstream>
#include <assert.h>

#define MAX_OPT_LENGTH 40
#define SHORT_OPT_LENGTH 6
#define COMPL_OPT_LENGTH MAX_OPT_LENGTH-SHORT_OPT_LENGTH

using namespace std;

namespace partest
{

  void PrintMeta::print_header (ostream& output)
  {
    output << endl;
    output << "--------------------------------------------" << endl;
    output << "|           PARTITIONTEST v" << PACKAGE_VERSION;
    output << "           |" << endl;
    output << "|                " << PROGRAM_DATE;
    output << "                |" << endl;
    output << "|          (c) Diego Darriba 2014          |" << endl;
    output << "|                                          |" << endl;
    output << "| Model selection for genomic alignments   |" << endl;
    output << "| PLL version (Stamatakis et al.)          |" << endl;
    output << "--------------------------------------------" << endl << endl;
  }

  void PrintMeta::print_options (ostream& output)
  {
    output << endl << ":: General Settings ::" << endl;
    output << setw (H_RULE_LENGTH) << setfill ('-') << "" << setfill (' ')
        << endl;
    output << setw (OPT_DESCR_LENGTH) << left << "  Input alignment:";
    if (input_file->length () > (H_RULE_LENGTH - OPT_DESCR_LENGTH))
    {
      output << endl << setw (H_RULE_LENGTH - (int) input_file->length ())
          << " ";
    }
    output << *input_file << endl;
    output << setw (OPT_DESCR_LENGTH) << left << "     Number of taxa:" << left
        << num_taxa << endl;
    output << setw (OPT_DESCR_LENGTH) << left << "     Number of sites (total):"
        << left << seq_len << endl;
    output << setw (OPT_DESCR_LENGTH) << left << "  Starting tree:";
    switch (starting_topology)
      {
      case StartTopoFIXED:
        output << "Fixed Maximum-Parsimony" << endl;
        break;
      case StartTopoFIXEDML:
        output << "Fixed Maximum-Likelihood" << endl;
        break;
      case StartTopoML:
        output << "Maximum-Likelihood" << endl;
        break;
      case StartTopoMP:
        output << "Maximum-Parsimony" << endl;
        break;
      case StartTopoUSER:
        output << "User Defined" << endl;
        break;
      default:
        cerr << "ERROR: Undefined starting topology" << endl;
        exit_partest (EX_SOFTWARE);
        break;
      }

    if (user_tree)
    {
      output << setw (OPT_DESCR_LENGTH) << left << "  Input tree:";
      if (user_tree->length () > 0)
      {
        if (user_tree->length () > (H_RULE_LENGTH - OPT_DESCR_LENGTH))
        {
          output << endl << setw (H_RULE_LENGTH - (int) user_tree->length ())
              << " ";
        }
        output << *user_tree << endl;
      }
      else
      {
        output << "N/A" << endl;
      }
    }

    output << setw (OPT_DESCR_LENGTH) << left << "  Config file:";
    if (config_file->length () > (H_RULE_LENGTH - OPT_DESCR_LENGTH))
    {
      output << endl << setw (H_RULE_LENGTH - (int) config_file->length ())
          << " ";
    }
    output << *config_file << endl;
    output << setw (OPT_DESCR_LENGTH) << left << "     Number of partitions:"
        << left << number_of_genes << endl;

    output << setw (OPT_DESCR_LENGTH) << left << "  Optimize mode:";
    switch (optimize_mode)
      {
      case OPT_GTR:
        if (data_type == DT_NUCLEIC)
        {
          output << "GTR+G" << endl;
        }
        else
        {
          output << "AUTO only" << endl;
        }
        break;
      case OPT_CUSTOM:
        output << "Custom set" << endl;
        break;
      case OPT_SEARCH:
        output << "Whole set" << endl;
        break;
      default:
        assert(0);
        break;
      }

    output << setw (OPT_DESCR_LENGTH) << left << "  Branch lengths:";
    if (reoptimize_branch_lengths)
    {
      output << "Unlinked" << endl;
    }
    else
    {
      output << "Proportional" << endl;
    }

    output << setw (OPT_DESCR_LENGTH) << left << "  Data type:";
    switch (data_type)
      {
      case DT_NUCLEIC:
        output << "Nucleic" << endl;
        if (optimize_mode != OPT_GTR)
        {
          output << setw (OPT_DESCR_LENGTH) << left
              << "  Include models with unequal frequencies:";
        }
        break;
      case DT_PROTEIC:
        output << "Proteic" << endl;
        if (optimize_mode != OPT_GTR)
        {
          output << setw (OPT_DESCR_LENGTH) << left
              << "  Include models with empirical frequencies:";
        }
        break;
      default:
        assert(0);
        break;
      }
    if (optimize_mode != OPT_GTR)
    {
      if (do_rate & RateVarF)
        output << "True " << endl;
      else
        output << "False" << endl;
#ifdef _IG_MODELS
      output << setw(OPT_DESCR_LENGTH) << left
      << "  Include models with rate variation:";
      if (do_rate & RateVarG)
      output << "True" << endl;
      else
      output << "False" << endl;
#endif
    }
    output << setw (OPT_DESCR_LENGTH) << left << "  Number of candidate models:"
        << number_of_models << endl;

    output << setw (OPT_DESCR_LENGTH) << left << "  Optimization epsilon:";
    if (epsilon == AUTO_EPSILON)
    {
      output << "Auto" << endl;
    }
    else
    {
      output << epsilon << endl;
    }

    if (number_of_schemes > 0)
    {
      output << setw (OPT_DESCR_LENGTH) << left << "  Schemes to evaluate:";
      output << left << number_of_schemes << endl;
      if (search_algo == SearchK1)
      {
        output << left << "Single partition only" << endl;
      }
      else if (search_algo == SearchKN)
      {
        output << left << "K=N scheme only" << endl;
      }
    }
    else
    {
      output << setw (OPT_DESCR_LENGTH) << left << "  Search algorithm:";
      switch (search_algo)
        {
        case SearchK1:
          output << left << "Single partition only" << endl;
          break;
        case SearchKN:
          output << left << "K=N scheme only" << endl;
          break;
        case SearchGreedy:
          output << left << "Greedy" << endl;
          break;
        case SearchGreedyExtended:
          output << left << "Greedy extended" << endl;
          break;
        case SearchHCluster:
          output << left << "Hierarchical Cluster (";
          if (max_samples)
            output << max_samples;
          else
            output << samples_percent * 100 << "%";
          output << ")" << endl;
          break;
        case SearchRandom:
          output << left << "Random" << endl;
          break;
        case SearchExhaustive:
          output << left << "Exahustive" << endl;
          break;
        case SearchAuto:
          output << left << "Auto" << endl;
          break;
        default:
          assert(0);
        }
    }
    output << setw (OPT_DESCR_LENGTH) << left << "  Max evaluations:";
    switch (search_algo)
      {
      case SearchK1:
        output << left << 1 << endl;
        break;
      case SearchKN:
        output << left << number_of_genes << endl;
        break;
      case SearchGreedy:
        output << left << Utilities::numSchemesGreedy (number_of_genes) << endl;
        break;
      case SearchGreedyExtended:
        output << left << Utilities::numSchemesGreedy (number_of_genes) << endl;
        break;
      case SearchHCluster:
        output << left
            << Utilities::numSchemesHierarchicalClustering (number_of_genes)
            << endl;
        break;
      case SearchRandom:
        output << left << max_samples << endl;
        break;
      case SearchExhaustive:
        output << left << "a lot" << endl;
        break;
      case SearchAuto:
        output << left << Utilities::numSchemesAutoSearch (number_of_genes)
            << endl;
        break;
      default:
        assert(0);
      }
    output << setw (OPT_DESCR_LENGTH) << left << "  Selection criterion:";
    switch (ic_type)
      {
      case BIC:
        output << left << "BIC" << endl;
        break;
      case AIC:
        output << left << "AIC" << endl;
        break;
      case AICC:
        output << left << "AICc" << endl;
        break;
      case DT:
        output << left << "DT" << endl;
        break;
      default:
        assert(0);
      }
    output << setw (OPT_DESCR_LENGTH) << left << "  Output path:";
    if (output_dir->length () > (H_RULE_LENGTH - OPT_DESCR_LENGTH))
    {
      output << endl << setw (H_RULE_LENGTH - (int) output_dir->length ())
          << " ";
    }
    output << *output_dir << endl;
    output << setw (H_RULE_LENGTH) << setfill ('-') << "" << setfill (' ')
        << endl << endl;
  }

  void PrintMeta::print_usage (std::ostream& out)
  {
    out << "Usage: " << PACKAGE << " -i sequenceFilename" << endl;
    out
        << "            [-c configFile] [-d nt|aa] [-F] [-h] [-N] [-O findModel|gtr]"
        << endl;
    out
        << "            [-p numberOfThreads] [-r numberOfReplicates] [-s aic|bic|aicc|dt] "
        << endl;
    out << "            [-S greedy|greedyext|hcluster|random|exhaustive]"
        << endl;
    out << "            [-t mp|fixed|user] [-u treeFile]" << endl;
    out << "            [--config-help] [--config-template]" << endl;
    out << endl;
    out << "Selects the best-fit model of amino acid or nucleotide replacement."
        << endl << endl;
    out
        << "Mandatory arguments for long options are also mandatory for short options."
        << endl;

    out << endl;
    out << " Input options:" << endl;

    out << setw (MAX_OPT_LENGTH) << left << "  -c, --config-file CONFIG_FILE"
        << "sets the input configuration file for gene partition" << endl;
    out << setw (MAX_OPT_LENGTH) << " "
        << "configuration file defines the gene bounds in the alignment"
        << endl;
    out << setw (MAX_OPT_LENGTH) << " "
        << "use argument --help-config for info about the format:" << endl;
    out << setw (SHORT_OPT_LENGTH) << " " << setw (COMPL_OPT_LENGTH)
        << "--config-help" << "shows help about configuration files" << endl;
    out << setw (SHORT_OPT_LENGTH) << " " << setw (COMPL_OPT_LENGTH)
        << "--config-template" << "generates a configuration file template"
        << endl;
    out << endl;

    out << setw (MAX_OPT_LENGTH) << left << "  -d, --data-type DATA_TYPE"
        << "sets the type of the input data" << endl;
    out << setw (SHORT_OPT_LENGTH) << " " << setw (COMPL_OPT_LENGTH)
        << "--data-type nt" << "nucleotide sequences" << endl;
    out << setw (SHORT_OPT_LENGTH) << " " << setw (COMPL_OPT_LENGTH)
        << "--data-type aa" << "amino-acid sequences" << endl;
    out << setw (MAX_OPT_LENGTH) << " "<< "default: nt" << endl;
    out << endl;

    out << setw (MAX_OPT_LENGTH) << left << "  -i, --input-file INPUT_FILE"
        << "sets the input alignment file (REQUIRED)" << endl;
    out << endl;

    out << setw (MAX_OPT_LENGTH) << left << "  -u, --user-tree TREE_FILE"
        << "sets a user-defined topology (this option ignores all" << endl;
    out << setw (MAX_OPT_LENGTH) << " "
        << "starting topologies different from \"user-defined\")" << endl;
    out << setw (MAX_OPT_LENGTH) << " " << "the tree must be in Newick format"
        << endl;
    out << endl;



    out << endl;
    out << " Model optimization options:" << endl;

    out << setw (MAX_OPT_LENGTH) << left << "  -e, --epsilon"
        << "sets model optimization epsilon (float / \"auto\")" << endl;
    out << setw (MAX_OPT_LENGTH) << " "<< "default: auto" << endl;
    out << endl;

    out << setw (MAX_OPT_LENGTH) << left << "  -F, --empirical-frequencies"
        << "includes models with empirical frequencies (+F)" << endl;
    out << endl;

    out << setw (MAX_OPT_LENGTH) << left << "  -O, --optimize OPTIMIZE_MODE"
        << "sets the model optimization for the best-fit partition" << endl;
    out << setw (SHORT_OPT_LENGTH) << " " << setw (COMPL_OPT_LENGTH)
        << "--optimize findModel"
        << "find the best-fit model for each partition" << endl;
    out << setw (SHORT_OPT_LENGTH) << " " << setw (COMPL_OPT_LENGTH)
        << "--optimize gtr"
        << "use only GTR model for each partition (nucleic data)" << endl;
    out << setw (MAX_OPT_LENGTH) << " " << "or AUTO for protein data." << endl;
    out << setw (MAX_OPT_LENGTH) << " "<< "default: findModel" << endl;
    out << endl;

    out << setw (MAX_OPT_LENGTH) << left << "  -t, --topology STARTING_TOPOLOGY"
        << "sets the starting topology for optimization" << endl;
    out << setw (SHORT_OPT_LENGTH) << " " << setw (COMPL_OPT_LENGTH)
        << "--topology mp"
        << "creates a maximum parsimony topology for each model optimization"
        << endl;
    out << setw (SHORT_OPT_LENGTH) << " " << setw (COMPL_OPT_LENGTH)
        << "--topology ml"
        << "creates a maximum likelihood topology for every model optimization"
        << endl;
    out << setw (SHORT_OPT_LENGTH) << " " << setw (COMPL_OPT_LENGTH)
        << "--topology fixed"
        << "uses a fixed MP topology for every model optimization" << endl;
    out << setw (SHORT_OPT_LENGTH) << " " << setw (COMPL_OPT_LENGTH)
        << "--topology fixedml"
        << "uses a fixed ML topology for every model optimization" << endl;
    out << setw (SHORT_OPT_LENGTH) << " " << setw (COMPL_OPT_LENGTH)
        << "--topology user"
        << "uses a user-defined topology. Requires the \"-u\" argument" << endl;
    out << setw (MAX_OPT_LENGTH) << " "
        << "however, if \"-u\" argument is used this option is automatically set"
        << endl;
    out << setw (MAX_OPT_LENGTH) << " "<< "default: mp" << endl;
    out << endl;

    out << endl;
    out << " Partitioning scheme search options:" << endl;

    out << setw (SHORT_OPT_LENGTH) << " " << setw (COMPL_OPT_LENGTH)
        << "--disable-ckp" << "disables the checkpointing" << endl;
    out << endl;

    out << setw (SHORT_OPT_LENGTH) << " " << setw (COMPL_OPT_LENGTH)
        << "--disable-output" << "disables any file-based output." << endl;
    out << setw (MAX_OPT_LENGTH) << " "
        << "this automatically disables also checkpointing" << endl;
    out << endl;

    out << setw (MAX_OPT_LENGTH) << left << "  -g, --pergene-bl"
        << "estimate per-gene branch lengths for starting topology" << endl;
    out << endl;

    out << setw (MAX_OPT_LENGTH) << left
        << "  -s, --selection-criterion CRITERION"
        << "sets the criterion for model selection" << endl;
    out << setw (MAX_OPT_LENGTH) << " "
        << "sample size for bic, aicc and dt criteria is the alignment length"
        << endl;
    out << setw (SHORT_OPT_LENGTH) << " " << setw (COMPL_OPT_LENGTH)
        << "--selection-criterion bic"
        << "Bayesian information criterion (DEFAULT)" << endl;
    out << setw (SHORT_OPT_LENGTH) << " " << setw (COMPL_OPT_LENGTH)
        << "--selection-criterion aic" << "Akaike information criterion"
        << endl;
    out << setw (SHORT_OPT_LENGTH) << " " << setw (COMPL_OPT_LENGTH)
        << "--selection-criterion aicc"
        << "corrected Akaike information criterion" << endl;
    out << setw (SHORT_OPT_LENGTH) << " " << setw (COMPL_OPT_LENGTH)
        << "--selection-criterion dt" << "decision theory" << endl;
    out << endl;

    out << setw (MAX_OPT_LENGTH) << left << "  -S, --search SEARCH_ALGORITHM"
        << "sets the search algorithm" << endl;
    out << setw (SHORT_OPT_LENGTH) << " " << setw (COMPL_OPT_LENGTH)
        << "--search k1" << "evaluate only partitioning scheme with K=1"
        << endl;
    out << setw (SHORT_OPT_LENGTH) << " " << setw (COMPL_OPT_LENGTH)
        << "--search kn" << "evaluate only partitioning scheme with K=N"
        << endl;
    out << setw (SHORT_OPT_LENGTH) << " " << setw (COMPL_OPT_LENGTH)
        << "--search greedy" << "greedy search algorithm" << endl;
    out << setw (SHORT_OPT_LENGTH) << " " << setw (COMPL_OPT_LENGTH)
        << "--search greedyext" << "extended greedy search algorithm" << endl;
    out << setw (SHORT_OPT_LENGTH) << " " << setw (COMPL_OPT_LENGTH)
        << "--search hcluster" << "hierarchical clustering algorithm" << endl;
    out << setw (SHORT_OPT_LENGTH) << " " << setw (COMPL_OPT_LENGTH)
        << "--search random" << "multiple step random sampling" << endl;
    out << setw (SHORT_OPT_LENGTH) << " " << setw (COMPL_OPT_LENGTH)
        << "--search auto" << "auto-select algorithm" << endl;
    out << setw (SHORT_OPT_LENGTH) << " " << setw (COMPL_OPT_LENGTH)
        << "--search exhaustive" << "exhaustive search" << endl;
    out << setw (MAX_OPT_LENGTH) << " "<< "default: auto" << endl;
    out << endl;

    out << setw (MAX_OPT_LENGTH) << left << "  -k, --keep-branches"
        << "keep branch lengths from the initial topology." << endl;
    out << setw (MAX_OPT_LENGTH) << " "
        << "this argument has no effect for initial topology different than fixed."
        << endl;
    out << endl;

    out << setw (SHORT_OPT_LENGTH) << "  -N" << setw (COMPL_OPT_LENGTH)
        << "--non-stop"
        << "algorithms do not stop if no improvement found at one step" << endl;
    out << endl;

    out << setw (MAX_OPT_LENGTH) << left
        << "  -p, --num-procs NUMBER_OF_THREADS"
        << "number of threads for model evaluation" << endl;
    out << endl;
    out << setw (MAX_OPT_LENGTH) << " " << "default: 1" << endl;

    out << setw (MAX_OPT_LENGTH) << left << "  -r, --replicates N"
        << "sets the number of replicates on hierarchical clustering" << endl;
    out << endl;
    out << setw (MAX_OPT_LENGTH) << " "<< "default: 1" << endl;

    out << setw (SHORT_OPT_LENGTH) << "  -T" << setw (COMPL_OPT_LENGTH)
        << "--get-final-tree" << "conduct final ML tree optimization" << endl;
    out << endl;

    out << setw (SHORT_OPT_LENGTH) << "  -w" << setw (COMPL_OPT_LENGTH)
        << "--weights W_R,W_F,W_A" << "sets the distances matrix weights (float) for hierarchical clustering" << endl;
    out << setw (MAX_OPT_LENGTH) << " "
        << "for W_R = subst.rates, W_F = base frequencies, W_A = Gamma shape (alpha)"
        << endl;
    out << setw (MAX_OPT_LENGTH) << " " << "default: 1.0,1.0,1.0" << endl;
    out << endl;


    out << endl;
    out << " Other options:" << endl;

    out << setw (SHORT_OPT_LENGTH) << " " << setw (COMPL_OPT_LENGTH)
        << "--force-override" << "existent output files will be overwritten."
        << endl;
    out << endl;

    out << setw (MAX_OPT_LENGTH) << left << "  -h, --help"
        << "displays this help message" << endl;
    out << endl;

    out << setw (MAX_OPT_LENGTH) << left << "  -v, --verbose LEVEL"
        << "sets the verbosity level (0=low, 1=medium, 2=high)" << endl;

    out << setw (MAX_OPT_LENGTH) << left << "  -V, --version"
        << "output version information and exit" << endl;
    out << endl;

  }

  void PrintMeta::print_results_xml (ostream & ofs,
                                     PartitioningScheme * bestScheme)
  {
    ofs << "<best_scheme num_elements=\"" << bestScheme->getNumberOfElements ()
        << "\" k=\"" << bestScheme->getNumberOfFreeParameters () << "\" lnL=\""
        << bestScheme->getLnL () << "\" BIC_score=\""
        << bestScheme->getBicValue () << "\" AIC_score=\""
        << bestScheme->getAicValue () << "\" AICC_score=\""
        << bestScheme->getAiccValue () << "\" linked_BIC_score=\""
        << bestScheme->getLinkedBicValue () << "\" linked_AIC_score=\""
        << bestScheme->getLinkedAicValue () << "\" linked_AICC_score=\""
        << bestScheme->getLinkedAiccValue () << "\">" << endl;
    ofs << "  <name>" << endl << "    " << bestScheme->getName () << endl
        << "  </name>" << endl;
    int numCodeLines = bestScheme->getCodeLines ();
    ofs << "  <code nlines=\"" << numCodeLines << "\">" << endl;
    for (int i = numCodeLines - 1; i >= 0; i--)
    {
      ofs << "    " << bestScheme->getCode (i) << endl;
    }
    ofs << "  </code>" << endl;
    if (bestScheme->getTree ())
    {
      ofs << "  <tree>" << endl;
      ofs << "    " << bestScheme->getTree () << endl;
      ofs << "  </tree>" << endl;
    }
    for (size_t i = 0; i < bestScheme->getNumberOfElements (); i++)
    {
      PartitionElement * element = bestScheme->getElement (i);
      ofs << "  <partition id=\"" << i + 1 << "\" num_elements=\""
          << element->getNumberOfSections () << "\" num_sites=\""
          << element->getNumberOfSites () << "\" num_patterns=\""
          << element->getNumberOfPatterns () << "\">" << endl;
      ofs << "    <name>" << endl << "      " << element->getName () << endl
          << "    </name>" << endl;
      ofs << "    <bestModel name=\""
          << element->getBestModel ()->getModel ()->getName () << "\" lnL=\""
          << element->getBestModel ()->getModel ()->getLnL () << "\" k=\""
          << element->getBestModel ()->getModel ()->getNumberOfFreeParameters ()
          << "\" BIC_score=\"" << element->getBestModel ()->getValue () << "\">"
          << endl;
      ofs << "    </bestModel>" << endl;
      ofs << "    <tree>" << endl;
      ofs << "      " << element->getBestModel ()->getModel ()->getTree ()
          << endl;
      ofs << "    </tree>" << endl;
      ofs << "  </partition>" << endl;
    }
    ofs << "  <raxml_control>" << endl;
    for (size_t i = 0; i < bestScheme->getNumberOfElements (); i++)
    {
      PartitionElement * element = bestScheme->getElement (i);
      switch (data_type)
        {
        case DT_NUCLEIC:
          ofs << "    DNA, ";
          break;
        case DT_PROTEIC:
          ofs << "    "
              << Utilities::getProtRaxmlName (
                  static_cast<ProteicModel *> (element->getBestModel ()->getModel ())->getMatrix ());
          if (element->getBestModel ()->getModel ()->isPF ())
          {
            ofs << "F";
          }
          ofs << ", ";

          break;
        default:
          assert(0);
        }

      vector<string> pllRegions (number_of_genes);
      pllQueueItem * qitem = pllPartsQueue->head;
      int j = 0;
      while (qitem)
      {
        pllPartitionInfo * pinfo = (pllPartitionInfo *) (qitem->item);
        pllQueueItem * ritem = pinfo->regionList->head;
        stringstream curpart;
        while (ritem)
        {
          pllPartitionRegion * pregion = (pllPartitionRegion *) ritem->item;
          curpart << pregion->start << "-" << pregion->end;
          if (pregion->stride > 1)
            curpart << "\\" << pregion->stride;
          ritem = ritem->next;
          if (ritem)
          {
            curpart << ",";
          }
        }
        pllRegions[j++] = curpart.str ();
        qitem = qitem->next;
      }

      ofs << "PART" << i + 1 << " = ";
      for (size_t j = 0; j < element->getNumberOfSections (); j++)
      {
        PEsection section = element->getSection (j);
        if (j > 0)
        {
          ofs << ", ";
        }
        ofs << pllRegions[section.id];
      }
      ofs << endl;
    }
    ofs << "  </raxml_control>" << endl;
    ofs << "</best_scheme>" << endl;
  }

  void PrintMeta::print_results (ostream & ofs, PartitioningScheme * bestScheme)
  {
    print_options (ofs);
    ofs << endl << ":: Selection results ::" << endl;
    ofs << setw (H_RULE_LENGTH) << setfill ('-') << "" << setfill (' ') << endl;
    bestScheme->print (ofs);
    ofs << setw (H_RULE_LENGTH) << setfill ('-') << "" << setfill (' ') << endl;
  }

} /* namespace partest */
