/*
 * PrintMeta.cc
 *
 *  Created on: 01/06/2012
 *      Author: diego
 */

#include "PrintMeta.h"
#include "parser/ArgumentParser.h"
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
	output << "|          (c) Diego Darriba 2012          |" << endl;
	output << "|                                          |" << endl;
	output << "| Model selection for genomic alignments   |" << endl;
#ifdef _PLL
	output << "| PLL version (Stamatakis et.al)           |" << endl;
#else
	output << "| PhyML version (Guindon and Gascuel 2011) |" << endl;
#endif
	output << "--------------------------------------------" << endl << endl;
}

void PrintMeta::print_options(ostream& output, ParTestOptions & options) {
	output << endl << ":: General Settings ::" << endl;
	output << setw(H_RULE_LENGTH) << setfill('-') << "" << setfill(' ') << endl;
	output << setw(OPT_DESCR_LENGTH) << left << "  Input alignment:";
	if (options.getInputFile().length() > (H_RULE_LENGTH - OPT_DESCR_LENGTH)) {
		output << endl << setw(H_RULE_LENGTH - options.getInputFile().length())
				<< " ";
	}
	output << options.getInputFile() << endl;
	output << setw(OPT_DESCR_LENGTH - 5) << left << "     Number of taxa:"
			<< setw(10) << right << options.getAlignment()->getNumSeqs()
			<< endl;
	output << setw(OPT_DESCR_LENGTH - 5) << left
			<< "     Number of sites (total):" << setw(10) << right
			<< options.getAlignment()->getNumSites() << endl;
	output << setw(OPT_DESCR_LENGTH - 5) << left
			<< "     Number of unique patterns:" << setw(10) << right
			<< options.getAlignment()->getNumPatterns() << endl;

	if (options.getTreeFile()) {
		output << setw(OPT_DESCR_LENGTH) << left << "  Input tree:";
		if (strlen(options.getTreeFile())
				> (H_RULE_LENGTH - OPT_DESCR_LENGTH)) {
			output << endl
					<< setw(H_RULE_LENGTH - strlen(options.getTreeFile()))
					<< " ";
		}
		output << options.getTreeFile() << endl;
	}

	output << setw(OPT_DESCR_LENGTH) << left << "  Config file:";
	if (strlen(options.getConfigFile()) > (H_RULE_LENGTH - OPT_DESCR_LENGTH)) {
		output << endl << setw(H_RULE_LENGTH - strlen(options.getConfigFile()))
				<< " ";
	}
	output << options.getConfigFile() << endl;
	output << setw(OPT_DESCR_LENGTH) << left << "  Output files:" << endl;
	output << setw(OPT_DESCR_LENGTH - 5) << left << "     Models:" << setw(10)
			<< right << options.getOutputFileModels() << endl;
	output << setw(OPT_DESCR_LENGTH - 5) << left << "     Partitions:"
			<< setw(10) << right << options.getOutputFilePartitions() << endl;
	output << setw(OPT_DESCR_LENGTH - 5) << left << "     Schemes:" << setw(10)
			<< right << options.getOutputFileSchemes() << endl;

	output << setw(OPT_DESCR_LENGTH) << left << "  Data type:";
	switch (options.getDataType()) {
	case DT_NUCLEIC:
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
	if (options.getRateVariation() & RateVarF)
		output << "True " << endl;
	else
		output << "False" << endl;
	output << setw(OPT_DESCR_LENGTH) << left
			<< "  Include models with invariant sites:";
	if (options.getRateVariation() & RateVarI)
		output << "True" << endl;
	else
		output << "False" << endl;
	output << setw(OPT_DESCR_LENGTH) << left
			<< "  Include models with rate variation:";
	if (options.getRateVariation() & RateVarG)
		output << "True" << endl;
	else
		output << "False" << endl;

//  output << setw(OPT_DESCR_LENGTH) << left << "  Candidate models:"
//      << Globals::get_number_of_models() << endl;
//  output << setw(OPT_DESCR_LENGTH) << left << "  Starting topology:";
//  switch ( run_instance.get_topology_type() ) {
//    case TT_BIONJ:
//      output << "BIONJ tree" << endl;
//      break;
//    case TT_FIXED:
//      output << "Fixed BIONJ tree" << endl;
//      break;
//    case TT_ML:
//      output << "ML tree" << endl;
//      break;
//    case TT_USER:
//      output << run_instance.get_user_tree() << endl;
//      break;
//  }
//  output << setw(OPT_DESCR_LENGTH) << left << "  Number of partitions:"
//      << run_instance.get_number_of_genes() << endl;
//  for ( unsigned int i = 0; i < run_instance.get_number_of_genes(); i++ ) {
//    output << "    (" << i + 1 << "/" << run_instance.get_number_of_genes() << ") "
//        << Globals::get_partition_name(Utilities::binary_pow(i)) << ": "
//        << run_instance.get_number_of_sites(Utilities::binary_pow(i)) << " sites" << endl;
//  }

	output << setw(H_RULE_LENGTH) << setfill('-') << "" << setfill(' ') << endl
			<< endl;
}

void PrintMeta::print_usage(std::ostream& out) {
	out << "Usage: " << PACKAGE << " -i INPUT_FILE [OPTION]..." << endl;
	out << endl;
	out << "Selects the best-fit model of amino acid or nucleotide replacement."
			<< endl << endl;
	out
			<< "Mandatory arguments for long options are also mandatory for short options."
			<< endl << endl;
	out << setw(MAX_OPT_LENGTH) << left << "  -h, --help"
			<< "Displays this help message" << endl;
	out << endl;
	out << setw(MAX_OPT_LENGTH) << left << "  -i, --input-file INPUT_FILE"
			<< "Sets the input alignment file (REQUIRED)" << endl;
	out << endl;
	out << setw(MAX_OPT_LENGTH) << left << "  -d, --data-type DATA_TYPE"
			<< "Sets the type of the input data" << endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--data-type nt" << "Nucleotide sequences (DEFAULT)" << endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--data-type aa" << "Amino-acid sequences" << endl;
	out << endl;
	out << setw(MAX_OPT_LENGTH) << left << "  -t, --topology STARTING_TOPOLOGY"
			<< "Sets the starting topology for optimization" << endl;
#ifndef _PLL
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--topology ml"
			<< "Creates a Maximum-Likelihood tree for each model optimization"
			<< endl;
	out << setw(MAX_OPT_LENGTH) << " "
			<< "(DEFAULT) Slowest but also more accurate" << endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--topology bionj"
			<< "Creates a BIONJ for each model optimization" << endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
				<< "--topology fixed"
				<< "Uses a fixed BIONJ for every model optimization" << endl;
#else
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
				<< "--topology mp"
				<< "(DEFAULT) Creates a maximum parsimony topology for each model optimization" << endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
				<< "--topology fixed"
				<< "Uses a fixed ML topology for every model optimization" << endl;
#endif
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--topology user"
			<< "Uses a user-defined topology. Requires the \"-u\" argument"
			<< endl;
	out << setw(MAX_OPT_LENGTH) << " "
			<< "However, if \"-u\" argument is used this option is automatically set"
			<< endl;
	out << endl;
	out << setw(MAX_OPT_LENGTH) << left << "  -u, --user-tree TREE_FILE"
			<< "Sets a user-defined topology. This option ignores all"
			<< endl;
	out << setw(MAX_OPT_LENGTH) << " "
			<< "starting topologies different from \"user-defined\""
			<< endl;
	out << setw(MAX_OPT_LENGTH) << " " << "The tree must be in Newick format"
			<< endl;
	out << endl;
#ifndef _PLL
	/** PLL does not support +I */
	out << setw(MAX_OPT_LENGTH) << left << "  -I, --invariant-sites"
			<< "Includes models with a proportion of invariant sites (+I)"
			<< endl;
	out << endl;
	/** PLL includes +G by default */
	out << setw(MAX_OPT_LENGTH) << left << "  -G, --gamma-rates"
			<< "Includes models with rate variation among sites (+G)" << endl;
	out << endl;
#endif
	out << setw(MAX_OPT_LENGTH) << left << "  -F, --empirical-frequencies"
			<< "Includes models with empirical frequencies (+F)" << endl;
	out << endl;
	out << setw(MAX_OPT_LENGTH) << left << "  -S, --search SEARCH_ALGORITHM"
			<< "Sets the search algorithm" << endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--search greedy" << "Greedy search algorithm (DEFAULT)" << endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--search greedyext"
			<< "Extended greedy search algorithm (DEFAULT)" << endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--search hcluster" << "Hierarchical clustering algorithm"
			<< endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--search random" << "Multiple step random sampling" << endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--search exhaustive" << "Exhaustive search" << endl;
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
	out << setw(MAX_OPT_LENGTH) << left
			<< "  -O, --optimize OPTIMIZE_MODE"
			<< "Sets the model optimization for the best-fit partition" << endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--optimize findModel"
			<< "Find the best-fit model for each partition" << endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
				<< "--optimize gtr"
				<< "Use GTR model for each partition" << endl;
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
	//  out << setw(MAX_OPT_LENGTH) << " " <<
//      << "  GENE_NAME=START,END[\\CODON_POSITION]"
//      << endl;
//  out << setw(MAX_OPT_LENGTH) << " " << "Example:" << endl;
//  out << setw(MAX_OPT_LENGTH) << " " << "  GENE1_1=1,150\\1" << endl;
//  out << setw(MAX_OPT_LENGTH) << " " << "  GENE1_2=1,150\\2" << endl;
//  out << setw(MAX_OPT_LENGTH) << " " << "  GENE1_3=1,150\\3" << endl;
//  out << setw(MAX_OPT_LENGTH) << " " << "  GENE2=151,300" << endl;
}

} /* namespace partest */
