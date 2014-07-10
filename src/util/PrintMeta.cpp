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

#include <iomanip>
#include <string.h>
#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#define PACKAGE "partest-tool"
#endif

#define MAX_OPT_LENGTH 40
#define SHORT_OPT_LENGTH 6
#define COMPL_OPT_LENGTH MAX_OPT_LENGTH-SHORT_OPT_LENGTH

using namespace std;

namespace partest {

void PrintMeta::print_header(ostream& output) {
	output << endl;
	output << "--------------------------------------------" << endl;
	output << "|               PARTITIONTEST              |" << endl;
	output << "|          (c) Diego Darriba 2014          |" << endl;
	output << "|                                          |" << endl;
	output << "| Model selection for genomic alignments   |" << endl;
	output << "| PLL version (Stamatakis et al.)          |" << endl;
	output << "--------------------------------------------" << endl << endl;
}

void PrintMeta::print_options(ostream& output) {
	output << endl << ":: General Settings ::" << endl;
	output << setw(H_RULE_LENGTH) << setfill('-') << "" << setfill(' ') << endl;
	output << setw(OPT_DESCR_LENGTH) << left << "  Input alignment:";
	if (input_file->length() > (H_RULE_LENGTH - OPT_DESCR_LENGTH)) {
		output << endl << setw(H_RULE_LENGTH - input_file->length()) << " ";
	}
	output << *input_file << endl;
	output << setw(OPT_DESCR_LENGTH) << left << "     Number of taxa:"
			<< left << num_taxa << endl;
	output << setw(OPT_DESCR_LENGTH) << left
			<< "     Number of sites (total):" << left << seq_len
			<< endl;
//	output << setw(OPT_DESCR_LENGTH - 5) << left
//			<< "     Number of unique patterns:" << left
//			<< num_patterns << endl;

	if (user_tree) {
		output << setw(OPT_DESCR_LENGTH) << left << "  Input tree:";
		if (user_tree->length() > 0) {
			if (user_tree->length() > (H_RULE_LENGTH - OPT_DESCR_LENGTH)) {
				output << endl << setw(H_RULE_LENGTH - user_tree->length())
						<< " ";
			}
			output << *user_tree << endl;
		} else {
			output << "N/A" << endl;
		}
	}

	output << setw(OPT_DESCR_LENGTH) << left << "  Config file:";
	if (config_file->length() > (H_RULE_LENGTH - OPT_DESCR_LENGTH)) {
		output << endl << setw(H_RULE_LENGTH - config_file->length()) << " ";
	}
	output << *config_file << endl;
	output << setw(OPT_DESCR_LENGTH) << left << "     Number of partitions:"
			<< left << number_of_genes << endl;
	output << setw(OPT_DESCR_LENGTH) << left << "  Data type:";
	switch (data_type) {
	case DT_NUCLEIC:
	case DT_DEFAULT:
		output << "Nucleic" << endl;
		output << setw(OPT_DESCR_LENGTH) << left
				<< "  Include models with unequal frequencies:";
		break;
	case DT_PROTEIC:
		output << "Proteic" << endl;
		output << setw(OPT_DESCR_LENGTH) << left
				<< "  Include models with empirical frequencies:";
		break;
	}
	if (do_rate & RateVarF)
		output << "True " << endl;
	else
		output << "False" << endl;
	output << setw(OPT_DESCR_LENGTH) << left
			<< "  Include models with rate variation:";
	if (do_rate & RateVarG)
		output << "True" << endl;
	else
		output << "False" << endl;
	output << setw(OPT_DESCR_LENGTH) << left << "  Number of candidate models:"
			<< number_of_models << endl;
	if (number_of_schemes > 0) {
		output << setw(OPT_DESCR_LENGTH) << left << "  Schemes to evaluate:";
		output << left << number_of_schemes << endl;
	} else {
		output << setw(OPT_DESCR_LENGTH) << left << "  Search algorithm:";
		switch (search_algo) {
		case SearchGreedy:
			output << left << "Greedy" << endl;
			break;
		case SearchGreedyExtended:
			output << left << "Greedy extended" << endl;
			break;
		case SearchHCluster:
			output << left << "Hierarchical Cluster (" << max_samples << ")"
					<< endl;
			break;
		case SearchRandom:
			output << left << "Random" << endl;
			break;
		case SearchExhaustive:
			output << left << "Exahustive" << endl;
			break;
		case SearchDefault:
		case SearchAuto:
			output << left << "Auto" << endl;
			break;
		}
	}
	output << setw(OPT_DESCR_LENGTH) << left << "  Selection criterion:";
	switch (ic_type) {
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
	}
	output << setw(H_RULE_LENGTH) << setfill('-') << "" << setfill(' ') << endl
			<< endl;
}

void PrintMeta::print_usage(std::ostream& out) {
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

	out << setw(MAX_OPT_LENGTH) << left << "  -c, --config-file CONFIG_FILE"
			<< "Sets the input configuration file for gene partition" << endl;
	out << setw(MAX_OPT_LENGTH) << " "
			<< "Configuration file defines the gene bounds in the alignment"
			<< endl;
	out << setw(MAX_OPT_LENGTH) << " "
			<< "Use argument --help-config for info about the format:" << endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--config-help" << "Shows help about configuration files"
			<< endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--config-template" << "Generates a configuration file template"
			<< endl;
	out << endl;

	out << setw(MAX_OPT_LENGTH) << left << "  -d, --data-type DATA_TYPE"
			<< "Sets the type of the input data" << endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--data-type nt" << "Nucleotide sequences (DEFAULT)" << endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--data-type aa" << "Amino-acid sequences" << endl;
	out << endl;

	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--disable-ckp" << "Disables the checkpointing"
			<< endl;
	out << endl;

	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--disable-output" << "Disables any file-based output." << endl;
	out << setw(MAX_OPT_LENGTH) << " "
				<< "This automatically disables also checkpointing"
				<< endl;
	out << endl;

	out << setw(MAX_OPT_LENGTH) << left << "  -F, --empirical-frequencies"
			<< "Includes models with empirical frequencies (+F)" << endl;
	out << endl;

	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--force-override" << "Existent output files will be overwritten." << endl;
	out << endl;

	out << setw(MAX_OPT_LENGTH) << left << "  -h, --help"
			<< "Displays this help message" << endl;
	out << endl;

	out << setw(MAX_OPT_LENGTH) << left << "  -i, --input-file INPUT_FILE"
			<< "Sets the input alignment file (REQUIRED)" << endl;
	out << endl;

	out << setw(MAX_OPT_LENGTH) << left << "  -k, --keep-branches"
				<< "Keep branch lengths from the initial topology." << endl;
		out << setw(MAX_OPT_LENGTH) << " "
				<< "This argument has no effect for initial topology different than fixed."
				<< endl;

	out << setw(SHORT_OPT_LENGTH) << "  -N" << setw(COMPL_OPT_LENGTH)
			<< "--non-stop"
			<< "Algorithms do not stop if no improvement found at one step"
			<< endl;
	out << endl;

	out << setw(MAX_OPT_LENGTH) << left << "  -O, --optimize OPTIMIZE_MODE"
			<< "Sets the model optimization for the best-fit partition" << endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--optimize findModel"
			<< "Find the best-fit model for each partition" << endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--optimize gtr"
			<< "Use only GTR model for each partition (nucleic data)" << endl;
	out << setw(MAX_OPT_LENGTH) << " " << "or AUTO for protein data." << endl;
	out << endl;

	out << setw(MAX_OPT_LENGTH) << left << "  -p, --num-procs NUMBER_OF_THREADS"
			<< "Number of threads for model evaluation (DEFAULT: 1)" << endl;
	out << endl;

	out << setw(MAX_OPT_LENGTH) << left << "  -r, --replicates N"
			<< "Sets the number of replicates on Hierarchical Clustering"
			<< endl;
	out << endl;

	out << setw(MAX_OPT_LENGTH) << left
			<< "  -s, --selection-criterion CRITERION"
			<< "Sets the criterion for model selection" << endl;
	out << setw(MAX_OPT_LENGTH) << " "
			<< "Sample size for bic, aicc and dt criteria is the alignment length"
			<< endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--selection-criterion bic"
			<< "Bayesian Information Criterion (DEFAULT)" << endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--selection-criterion aic" << "Akaike Information Criterion"
			<< endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--selection-criterion aicc"
			<< "Corrected Akaike Information Criterion" << endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--selection-criterion dt" << "Decision Theory" << endl;
	out << endl;

	out << setw(MAX_OPT_LENGTH) << left << "  -S, --search SEARCH_ALGORITHM"
			<< "Sets the search algorithm" << endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--search greedy" << "Greedy search algorithm" << endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--search greedyext" << "Extended greedy search algorithm"
			<< endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--search hcluster" << "Hierarchical clustering algorithm"
			<< endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--search random" << "Multiple step random sampling" << endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--search auto" << "Auto-select algorithm (DEFAULT)" << endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
				<< "--search exhaustive" << "Exhaustive search" << endl;
	out << endl;

	out << setw(MAX_OPT_LENGTH) << left << "  -t, --topology STARTING_TOPOLOGY"
			<< "Sets the starting topology for optimization" << endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--topology " << ARG_TOPO_MP
			<< "(DEFAULT) Creates a maximum parsimony topology for each model optimization"
			<< endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--topology " << ARG_TOPO_ML
			<< "Creates a maximum likelihood topology for every model optimization" << endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--topology "<< ARG_TOPO_FIXED
			<< "Uses a fixed ML topology for every model optimization" << endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--topology "<< ARG_TOPO_USER
			<< "Uses a user-defined topology. Requires the \"-u\" argument"
			<< endl;
	out << setw(MAX_OPT_LENGTH) << " "
			<< "However, if \"-u\" argument is used this option is automatically set"
			<< endl;
	out << endl;

	out << setw(SHORT_OPT_LENGTH) << "  -T" << setw(COMPL_OPT_LENGTH)
			<< "--get-final-tree"
			<< "Conduct final ML tree optimization"
			<< endl;
	out << endl;

	out << setw(MAX_OPT_LENGTH) << left << "  -u, --user-tree TREE_FILE"
			<< "Sets a user-defined topology. This option ignores all" << endl;
	out << setw(MAX_OPT_LENGTH) << " "
			<< "starting topologies different from \"user-defined\"" << endl;
	out << setw(MAX_OPT_LENGTH) << " " << "The tree must be in Newick format"
			<< endl;
	out << endl;

	out << setw(MAX_OPT_LENGTH) << left << "  -v, --version"
			<< "Output version information and exit" << endl;
	out << endl;

}

void PrintMeta::print_results_xml(ostream & ofs,
		PartitioningScheme * bestScheme) {
	ofs << "<best_scheme num_elements=\"" << bestScheme->getNumberOfElements()
			<< "\" k=\"" << bestScheme->getNumberOfFreeParameters()
			<< "\" lnL=\"" << bestScheme->getLnL()
			<< "\" BIC_score=\"" << bestScheme->getBicValue()
			<< "\" AIC_score=\"" << bestScheme->getAicValue()
			<< "\" AICC_score=\"" << bestScheme->getAiccValue()
			<< "\" linked_BIC_score=\"" << bestScheme->getLinkedBicValue()
			<< "\" linked_AIC_score=\"" << bestScheme->getLinkedAicValue()
			<< "\" linked_AICC_score=\"" << bestScheme->getLinkedAiccValue()
			<< "\">" << endl;
	ofs << "  <name>" << endl << "    " << bestScheme->getName() << endl
			<< "  </name>" << endl;
	int numCodeLines = bestScheme->getCodeLines();
	ofs << "  <code nlines=\"" << numCodeLines << "\">" << endl;
	for (int i = numCodeLines - 1; i >= 0; i--) {
		ofs << "    " << bestScheme->getCode(i) << endl;
	}
	ofs << "  </code>" << endl;
	if (bestScheme->getTree()) {
		ofs << "  <tree>" << endl;
		ofs << "    " << bestScheme->getTree() << endl;
		ofs << "  </tree>" << endl;
	}
	for (size_t i = 0; i < bestScheme->getNumberOfElements(); i++) {
		PartitionElement * element = bestScheme->getElement(i);
		ofs << "  <partition id=\"" << i + 1 << "\" num_elements=\""
				<< element->getNumberOfSections() << "\" num_sites=\""
				<< element->getNumberOfSites() << "\" num_patterns=\""
				<< element->getNumberOfPatterns() << "\">" << endl;
		ofs << "    <name>" << endl << "      " << element->getName() << endl
				<< "    </name>" << endl;
		ofs << "    <bestModel name=\""	<< element->getBestModel()->getModel()->getName()
				<< "\" lnL=\"" << element->getBestModel()->getModel()->getLnL()
				<< "\" k=\"" << element->getBestModel()->getModel()->getNumberOfFreeParameters()
				<< "\" BIC_score=\"" << element->getBestModel()->getValue()
				<< "\">" << endl;
		ofs << "    </bestModel>" << endl;
		ofs << "    <tree>" << endl;
		ofs << "      " << element->getBestModel()->getModel()->getTree()
				<< endl;
		ofs << "    </tree>" << endl;
		ofs << "  </partition>" << endl;
	}
	ofs << "  <raxml_control>" << endl;
	for (size_t i = 0; i < bestScheme->getNumberOfElements(); i++) {
		PartitionElement * element = bestScheme->getElement(i);
		switch (data_type) {
		case DT_NUCLEIC:
			ofs << "    DNA, ";
			break;
		case DT_PROTEIC:
			ofs << "    "
					<< Utilities::getProtRaxmlName(
							static_cast<ProteicModel *>(element->getBestModel()->getModel())->getMatrix())
					<< ", ";

			break;
		case DT_DEFAULT:
			cerr << "Uninitialized DataType" << endl;
			exit_partest(EX_SOFTWARE);
		}
		ofs << "PART" << i + 1 << " = ";
		for (int j = 0; j < element->getNumberOfSections(); j++) {
			PEsection section = element->getSection(j);
			if (j > 0) {
				ofs << ",";
			}
			ofs << section.start << "-" << section.end;
		}
		ofs << endl;
	}
	ofs << "  </raxml_control>" << endl;
	ofs << "</best_scheme>" << endl;
}

void PrintMeta::print_results(ostream & ofs, PartitioningScheme * bestScheme) {
	print_options(ofs);
	ofs << endl << ":: Selection results ::" << endl;
	ofs << setw(H_RULE_LENGTH) << setfill('-') << "" << setfill(' ') << endl;
	bestScheme->print(ofs);
	ofs << setw(H_RULE_LENGTH) << setfill('-') << "" << setfill(' ') << endl;
}

} /* namespace partest */
