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
	output << "|          (c) Diego Darriba 2014          |" << endl;
	output << "|                                          |" << endl;
	output << "| Model selection for genomic alignments   |" << endl;
	output << "| PLL version (Stamatakis et.al)           |" << endl;
	output << "--------------------------------------------" << endl << endl;
}

void PrintMeta::print_options(ostream& output) {
	output << endl << ":: General Settings ::" << endl;
		output << setw(H_RULE_LENGTH) << setfill('-') << "" << setfill(' ') << endl;
		output << setw(OPT_DESCR_LENGTH) << left << "  Input alignment:";
		if (input_file->length() > (H_RULE_LENGTH - OPT_DESCR_LENGTH)) {
			output << endl << setw(H_RULE_LENGTH - input_file->length())
					<< " ";
		}
		output << *input_file << endl;
		output << setw(OPT_DESCR_LENGTH - 5) << left << "     Number of taxa:"
				<< setw(10) << right << num_taxa
				<< endl;
		output << setw(OPT_DESCR_LENGTH - 5) << left
				<< "     Number of sites (total):" << setw(10) << right
				<< seq_len << endl;
		output << setw(OPT_DESCR_LENGTH - 5) << left
				<< "     Number of unique patterns:" << setw(10) << right
				<< num_patterns << endl;

		if (user_tree) {
			output << setw(OPT_DESCR_LENGTH) << left << "  Input tree:";
			if (user_tree->length() > 0) {
				if (user_tree->length()
						> (H_RULE_LENGTH - OPT_DESCR_LENGTH)) {
					output << endl
							<< setw(H_RULE_LENGTH - user_tree->length())
							<< " ";
				}
				output << *user_tree << endl;
			} else {
				output << "N/A" << endl;
			}
		}

		output << setw(OPT_DESCR_LENGTH) << left << "  Config file:";
		if (config_file->length() > (H_RULE_LENGTH - OPT_DESCR_LENGTH)) {
			output << endl << setw(H_RULE_LENGTH - config_file->length())
					<< " ";
		}
		output << *config_file << endl;

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
		output << setw(OPT_DESCR_LENGTH) << left << "  Search algorithm:";
		switch (search_algo) {
		case SearchGreedy:
		case SearchDefault:
			output << left << "Greedy" << endl;
			break;
		case SearchGreedyExtended:
			output << left << "Greedy extended" << endl;
			break;
		case SearchHCluster:
			output << left << "Hierarchical Cluster (" << max_samples
					<< ")" << endl;
			break;
		case SearchRandom:
			output << left << "Random" << endl;
			break;
		case SearchExhaustive:
			output << left << "Exahustive" << endl;
			break;
		}
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
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--topology mp"
			<< "(DEFAULT) Creates a maximum parsimony topology for each model optimization"
			<< endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--topology fixed"
			<< "Uses a fixed ML topology for every model optimization" << endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--topology user"
			<< "Uses a user-defined topology. Requires the \"-u\" argument"
			<< endl;
	out << setw(MAX_OPT_LENGTH) << " "
			<< "However, if \"-u\" argument is used this option is automatically set"
			<< endl;
	out << endl;
	out << setw(MAX_OPT_LENGTH) << left << "  -u, --user-tree TREE_FILE"
			<< "Sets a user-defined topology. This option ignores all" << endl;
	out << setw(MAX_OPT_LENGTH) << " "
			<< "starting topologies different from \"user-defined\"" << endl;
	out << setw(MAX_OPT_LENGTH) << " " << "The tree must be in Newick format"
			<< endl;
	out << endl;
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
	out << setw(MAX_OPT_LENGTH) << left << "  -r, --replicates N"
				<< "Sets the number of replicates on Hierarchical Clustering" << endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
				<< "--non-stop" << "Algorithms do not stop if no improvement found at one step"
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
	out << setw(MAX_OPT_LENGTH) << left << "  -O, --optimize OPTIMIZE_MODE"
			<< "Sets the model optimization for the best-fit partition" << endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--optimize findModel"
			<< "Find the best-fit model for each partition" << endl;
	out << setw(SHORT_OPT_LENGTH) << " " << setw(COMPL_OPT_LENGTH)
			<< "--optimize gtr" << "Use GTR model for each partition" << endl;
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
