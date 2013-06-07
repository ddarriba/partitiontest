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
#include "util/PrintMeta.h"
#include "indata/Alignment.h"

using namespace std;

namespace partest {

#define NUM_ARGUMENTS 16

void ArgumentParser::init() {
	index = 1;
	subindex = 1;
}

ArgumentParser::ArgumentParser() {
	option options_list[] = { { ARG_HELP, 'h', "help", false }, {
			ARG_INPUT_FILE, 'i', "input-file", true }, {
			ARG_INPUT_FORMAT, 'f', "input-format", true }, {
			ARG_USER_TREE, 'u', "user-tree", true }, {
			ARG_DATA_TYPE, 'd', "data-type", true }, {
			ARG_FREQUENCIES, 'F', "empirical-frequencies", false }, {
			ARG_INV, 'I', "invariant-sites", false }, {
			ARG_GAMMA, 'G', "gamma-rates", false }, {
			ARG_TOPOLOGY, 't', "topology", true }, {
			ARG_CONFIG_FILE, 'c', "config-file", true }, {
			ARG_SEARCH_ALGORITHM, 'S', "search", true }, {
			ARG_NUM_PROCS, 'p', "num-procs", true }, {
			ARG_IC_TYPE, 's', "selection-criterion", true }, {
			ARG_SAMPLE_SIZE, 'n', "sample-size", true }, {
			ARG_CONFIG_HELP, 0, "config-help", false }, {
			ARG_CONFIG_TEMPLATE, 0, "config-template", false } };
	unsigned int size = NUM_ARGUMENTS * sizeof(option);
	arguments = (option *) malloc(size);
	memcpy(arguments, options_list, size);
	init();
}

ArgumentParser::~ArgumentParser() {
	free(arguments);
}

ArgIndex ArgumentParser::get_opt(int argc, char *argv[], char *argument,
		char *value) {
	if (index >= argc)
		return ARG_END;
	else {
		ArgIndex arg_index = ARG_NULL;
		char *arg = argv[index];
		strcpy(argument, argv[index]);
		if (*arg == '-') {
			if (*(arg + 1) == '-') {
				// check long option
				for (int j = 0; j < NUM_ARGUMENTS; j++) {
					if (arguments[j].long_code > 0) {
						if (!strcmp(arguments[j].long_code, &(arg[2]))) {
							arg_index = arguments[j].index;
							if (arguments[j].required_value) {
								index++;
								if (index >= argc || argv[index][0] == '-') {
									cerr << "ERROR! Argument " << arg
											<< " requires a value." << endl;
									Utilities::exit_partest(EX_USAGE);
								}
								strcpy(value, argv[index]);
							}
						}
					}
				}
				index++;
			} else {
				// check short option
				int num_arguments = strlen(arg + 1);
				for (int j = 0; j < NUM_ARGUMENTS; j++) {
					if (arguments[j].char_code == arg[subindex]) {
						arg_index = arguments[j].index;
						if (arguments[j].required_value) {
							if (num_arguments > 1) {
								cerr << "ERROR! Argument " << arg[subindex]
										<< " cannot be used in a row." << endl;
								Utilities::exit_partest(EX_USAGE);
							} else {
								index++;
								if (index >= argc || argv[index][0] == '-') {
									cerr << "ERROR! Argument " << arg[subindex]
											<< " requires a value." << endl;
									Utilities::exit_partest(EX_USAGE);
								}
								strcpy(value, argv[index]);
							}
						}
					}
				}
				if (num_arguments > subindex) {
					subindex++;
				} else {
					subindex = 1;
					index++;
				}
				if (arg_index == ARG_NULL) {
					cerr << "[ERROR] \"" << argument
							<< "\" option not recognized." << endl;
					cerr << "Run with --help for more information." << endl;
					Utilities::exit_partest(EX_USAGE);
				}
			}
		} else {
			// NOT AN OPTION!
			cerr << "[ERROR] \"" << argv[index] << "\" is not a valid option."
					<< endl;
			cerr << "Arguments must start with '-' (short) or '--' (long)"
					<< endl;
			cerr << "Run with --help for more information." << endl;
			Utilities::exit_partest(EX_USAGE);
		}
		return arg_index;
	}
}

void ArgumentParser::fill_options(int argc, char *argv[],
		ParTestOptions *options) {

	bool requiredInputFile = false;

	init();

	int argument_index;
	char value[256];
	char argument[256];
	char input_file[256] = "";
	char user_tree[256] = "";
	char config_file[256] = "";
	bool do_f = DEFAULT_DO_F;
	bool do_i = DEFAULT_DO_I;
	bool do_g = DEFAULT_DO_G;
	double sampleSizeValue;

	DataType data_type = DEFAULT_DATA_TYPE;
	StartTopo startingTopology = DEFAULT_STARTING_TOPOLOGY;
	InformationCriterion ic_type = DEFAULT_IC_TYPE;
	SampleSize sampleSize = DEFAULT_SAMPLE_SIZE;
	SearchAlgo searchAlgo = DEFAULT_SEARCH_ALGO;
//  AlignFormat input_format = AF_PHYLIP_SEQ;

	while ((argument_index = get_opt(argc, argv, argument, value)) != ARG_END) {
		switch (argument_index) {
		case ARG_HELP:
			PrintMeta::print_usage(cout);
			Utilities::exit_partest(EX_OK);
			break;
		case ARG_INPUT_FILE:
			strcpy(input_file, value);
			requiredInputFile = true;
			break;
		case ARG_INPUT_FORMAT:
			//        if ( !strcmp(value, ARG_AF_PHYLIP_SEQ) ) {
//          input_format = AF_PHYLIP_SEQ;
//        } else if ( !strcmp(value, ARG_AF_PHYLIP_INT) ) {
//          input_format = AF_PHYLIP_INT;
//        } else if ( !strcmp(value, ARG_AF_NEXUS) ) {
//          input_format = AF_NEXUS;
//        } else if ( !strcmp(value, ARG_AF_CLUSTAL) ) {
//          input_format = AF_CLUSTAL;
//        } else if ( !strcmp(value, ARG_AF_FASTA) ) {
//          input_format = AF_FASTA;
//        } else if ( !strcmp(value, ARG_AF_DCSE) ) {
//          input_format = AF_DCSE;
//        } else if ( !strcmp(value, ARG_AF_GENBANK) ) {
//          input_format = AF_GENBANK;
//        } else {
//          //***ERROR!
//          Utilities::exit_partest(EX_USAGE);
//        }
			break;
		case ARG_USER_TREE:
			cout << "User Tree " << value << endl;
			strcpy(user_tree, value);
			break;
		case ARG_DATA_TYPE:
			if (!strcmp(value, ARG_DT_PROTEIC)) {
				data_type = DT_PROTEIC;
			} else if (!strcmp(value, ARG_DT_NUCLEIC)) {
				data_type = DT_NUCLEIC;
			} else {
				cerr << "[ERROR] \"-d " << value
						<< "\" is not a valid data type. Use one of the following:"
						<< endl;
				cerr << "  -d " << setw(8) << left << ARG_DT_PROTEIC
						<< "Protein sequence alignment" << endl;
				cerr << "  -d " << setw(8) << left << ARG_DT_NUCLEIC
						<< "DNA sequence alignment (DEFAULT)" << endl;
				Utilities::exit_partest(EX_USAGE);
			}
			break;
		case ARG_TOPOLOGY:
			if (!strcmp(value, ARG_TOPO_BIONJ)) {
				startingTopology = StartTopoBIONJ;
			} else if (!strcmp(value, ARG_TOPO_ML)) {
				startingTopology = StartTopoML;
			} else if (!strcmp(value, ARG_TOPO_FIXED)) {
				startingTopology = StartTopoFIXED;
			} else if (!strcmp(value, ARG_TOPO_USER)) {
				startingTopology = StartTopoUSER;
			} else {
				cerr << "[ERROR] \"-t " << value
						<< "\" is not a valid input topology. Use one of the following:"
						<< endl;
				cerr << "  -t " << setw(8) << left << ARG_TOPO_BIONJ
						<< "BIONJ topology" << endl;
				cerr << "  -t " << setw(8) << left << ARG_TOPO_FIXED
						<< "Fixed BIONJ topology for every model" << endl;
				cerr << "  -t " << setw(8) << left << ARG_TOPO_ML
						<< "Maximum Likelihood topology (DEFAULT, slowest)"
						<< endl;
				cerr << "  -t " << setw(8) << left << ARG_TOPO_USER
						<< "User-defined topology" << endl;
				Utilities::exit_partest(EX_USAGE);
			}
			break;
		case ARG_SEARCH_ALGORITHM:
			if (!strcmp(value, ARG_SEARCH_EXHAUSTIVE)) {
				searchAlgo = SearchExhaustive;
			} else if (!strcmp(value, ARG_SEARCH_RANDOM)) {
				searchAlgo = SearchRandom;
			} else if (!strcmp(value, ARG_SEARCH_GREEDY)) {
				searchAlgo = SearchGreedy;
			} else if (!strcmp(value, ARG_SEARCH_GREEDY_EXT)) {
							searchAlgo = SearchGreedyExtended;
			} else if (!strcmp(value, ARG_SEARCH_HIERARCHICAL)) {
				searchAlgo = SearchHCluster;
			} else {
				cerr << "[ERROR] \"-S " << value
						<< "\" is not a valid search algorithm. Use one of the following:"
						<< endl;
				cerr << "  -S " << setw(12) << left << ARG_SEARCH_EXHAUSTIVE
						<< "Exhaustive algorithm (horribly computationally expensive)"
						<< endl;
				cerr << "  -S " << setw(12) << left << ARG_SEARCH_RANDOM
						<< "Random walk algorithm (chinese restaurant process)"
						<< endl;
				cerr << "  -S " << setw(12) << left << ARG_SEARCH_GREEDY
						<< "Greedy hill climbing algorithm" << endl;
				cerr << "  -S " << setw(12) << left << ARG_SEARCH_GREEDY_EXT
						<< "Extended greedy hill climbing algorithm" << endl;
				cerr << "  -S " << setw(12) << left << ARG_SEARCH_HIERARCHICAL
						<< "Hierarchical clustering algorithm" << endl;
				Utilities::exit_partest(EX_USAGE);
			}
			break;
		case ARG_IC_TYPE:
			if (!strcmp(value, ARG_IC_AIC)) {
				ic_type = AIC;
			} else if (!strcmp(value, ARG_IC_BIC)) {
				ic_type = BIC;
			} else if (!strcmp(value, ARG_IC_AICC)) {
				ic_type = AICC;
			} else if (!strcmp(value, ARG_IC_DT)) {
				ic_type = DT;
			} else {
				cerr << "[ERROR] \"-s " << value
						<< "\" is not a valid criterion. Use one of the following:"
						<< endl;
				cerr << "  -s " << setw(8) << left << ARG_IC_AIC
						<< "Akaike Information Criterion" << endl;
				cerr << "  -s " << setw(8) << left << ARG_IC_BIC
						<< "Bayesian Information Criterion" << endl;
				cerr << "  -s " << setw(8) << left << ARG_IC_AICC
						<< "Corrected Akaike Information Criterion" << endl;
				cerr << "  -s " << setw(8) << left << ARG_IC_DT
						<< "Decision Theory" << endl;
				Utilities::exit_partest(EX_USAGE);
			}
			break;
		case ARG_SAMPLE_SIZE:
			sampleSizeValue = 0.0;
			if (!strcmp(value, ARG_SS_ALIGN)) {
				sampleSize = SS_ALIGNMENT;
			} else if (!strcmp(value, ARG_SS_SHANNON)) {
				sampleSize = SS_SHANNON;
			} else {
				sampleSizeValue = atof(value);
				if (!sampleSizeValue) {
					cerr << "[ERROR] \"-n " << value
							<< "\" is not a valid sample size. Use one of the following:"
							<< endl;
					cerr << "  -n " << setw(16) << left << ARG_SS_ALIGN
							<< "Alignment size" << endl;
					cerr << "  -n " << setw(16) << left << ARG_SS_SHANNON
							<< "Shannon entropy of the alignment" << endl;
					cerr << "  -n " << setw(16) << left << "[CUSTOM_VALUE]"
							<< "Custom sample size" << endl;
					Utilities::exit_partest(EX_USAGE);
				}
				sampleSize = SS_CUSTOM;
			}
			break;
		case ARG_FREQUENCIES:
			do_f = true;
			break;
		case ARG_INV:
			do_i = true;
			break;
		case ARG_GAMMA:
			do_g = true;
			break;
		case ARG_CONFIG_FILE:
			strcpy(config_file, value);
			break;
		case ARG_CONFIG_HELP:
			ConfigParser::printFormat();
			Utilities::exit_partest(0);
			break;
		case ARG_CONFIG_TEMPLATE:
			ConfigParser::createTemplate();
			Utilities::exit_partest(0);
			break;
		case ARG_NUM_PROCS:
#ifdef PTHREADS
			Globals::set_number_of_threads(atoi(value));
#else
			cerr
					<< "[ERROR] PThreads version is not available. You must recompile with PTHREADS flag."
					<< endl;
			Utilities::exit_partest(EX_USAGE);
#endif
			break;
		default:
			cerr << "[ERROR] \"" << argument
					<< "\" option not recognized.\n Run with --help for more information."
					<< endl;
			Utilities::exit_partest(EX_USAGE);
			break;
		}
	}

#ifdef DEBUG
	cout << "[TRACE] All arguments were parsed" << endl;
#endif

	// check required arguments
	if (!requiredInputFile) {
		cerr << "ERROR! Input File (-i) is required!" << endl;
		Utilities::exit_partest(EX_USAGE);
	}

	bitMask do_rate = RateVarM;
	if (do_f)
		do_rate |= (RateVarF);
	if (do_i)
		do_rate |= (RateVarI);
	if (do_g)
		do_rate |= (RateVarG);

#ifdef DEBUG
	cout << "[TRACE] Attempting to set ParTest Options  " << endl;
#endif

	options->set(input_file, data_type, do_rate, config_file, startingTopology,
			searchAlgo, ic_type, sampleSize, sampleSizeValue, user_tree);

#ifdef DEBUG
	cout << "[TRACE] ParTest options set" << endl;
#endif
}

} /* namespace partest */
