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

using namespace std;

namespace partest {

#define NUM_ARGUMENTS 18

void ArgumentParser::init() {
	index = 1;
	subindex = 1;
}

ArgumentParser::ArgumentParser(PartitionTest * ptest) :
		ptest(ptest), index(0), subindex(0) {

	option options_list[] = { { ARG_HELP, 'h', "help", false }, {
			ARG_INPUT_FILE, 'i', "input-file", true }, {
			ARG_INPUT_FORMAT, 'f', "input-format", true }, {
			ARG_USER_TREE, 'u', "user-tree", true }, {
			ARG_DATA_TYPE, 'd', "data-type", true }, {
			ARG_FREQUENCIES, 'F', "empirical-frequencies", false }, {
#ifndef _PLL
			ARG_INV, 'I', "invariant-sites", false }, {
			ARG_GAMMA, 'G', "gamma-rates", false }, {
#endif
			ARG_TOPOLOGY, 't', "topology", true }, {
			ARG_CONFIG_FILE, 'c', "config-file", true }, {
			ARG_SEARCH_ALGORITHM, 'S', "search", true }, {
			ARG_HCLUSTER_REPS, 'r', "replicates", true }, {
			ARG_OPTIMIZE, 'O', "optimize", true }, {
			ARG_OUTPUT, 'o', "output", true }, {
			ARG_NUM_PROCS, 'p', "num-procs", true }, {
			ARG_IC_TYPE, 's', "selection-criterion", true }, {
			ARG_SAMPLE_SIZE, 'n', "sample-size", true }, {
			ARG_CONFIG_HELP, 0, "config-help", false }, {
			ARG_CONFIG_TEMPLATE, 0, "config-template", false }, {
			ARG_NON_STOP, 0, "non-stop", false } };
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
									exit_partest(EX_CONFIG);
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
								exit_partest(EX_CONFIG);
							} else {
								index++;
								if (index >= argc || argv[index][0] == '-') {
									cerr << "ERROR! Argument " << arg[subindex]
											<< " requires a value." << endl;
									exit_partest(EX_CONFIG);
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
					exit_partest(EX_USAGE);
				}
			}
		} else {
			// NOT AN OPTION!
			cerr << "[ERROR] \"" << argv[index] << "\" is not a valid option."
					<< endl;
			cerr << "Arguments must start with '-' (short) or '--' (long)"
					<< endl;
			cerr << "Run with --help for more information." << endl;
			exit_partest(EX_USAGE);
		}

		return arg_index;
	}
}

void ArgumentParser::parse(int argc, char *argv[]) {

	bool requiredInputFile = false;
	init();
	int argument_index;
	char value[256];
	char argument[256];
	char input_file[256] = "";
	char user_tree[256] = "";
	char config_file[256] = "";
	char output_dir[256] = "";
	bool do_f = DEFAULT_DO_F;
	bool do_i = DEFAULT_DO_I;
	bool do_g = DEFAULT_DO_G;
	double sampleSizeValue = 0.0;

	DataType data_type = DT_DEFAULT;
	StartTopo startingTopology = DEFAULT_STARTING_TOPOLOGY;
	InformationCriterion ic_type = IC_DEFAULT;
	SampleSize sampleSize = SS_DEFAULT;
	SearchAlgo searchAlgo = SearchDefault;
	int maxSamples = 1;
	OptimizeMode optimize = OPT_DEFAULT;

	while ((argument_index = get_opt(argc, argv, argument, value)) != ARG_END) {
		switch (argument_index) {
		case ARG_HELP:
			exit_partest(EX_USAGE);
			break;
		case ARG_INPUT_FILE:
			strcpy(input_file, value);
			requiredInputFile = true;
			break;
		case ARG_USER_TREE:
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
				exit_partest(EX_CONFIG);
			}
			break;
		case ARG_TOPOLOGY:
			if (!strcmp(value, ARG_TOPO_MP)) {
				startingTopology = StartTopoMP;
			} else if (!strcmp(value, ARG_TOPO_FIXED)) {
				startingTopology = StartTopoFIXED;
			} else if (!strcmp(value, ARG_TOPO_USER)) {
				startingTopology = StartTopoUSER;
			} else {
				cerr << "[ERROR] \"-t " << value
						<< "\" is not a valid input topology. Use one of the following:"
						<< endl;
				cerr << "  -t " << setw(8) << left << ARG_TOPO_MP
						<< "MP topology" << endl;
				cerr << "  -t " << setw(8) << left << ARG_TOPO_FIXED
						<< "Fixed ML topology for every model" << endl;
				cerr << "  -t " << setw(8) << left << ARG_TOPO_USER
						<< "User-defined topology" << endl;
				exit_partest(EX_CONFIG);
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
						<< "Random walk algorithm (Chinese restaurant process)"
						<< endl;
				cerr << "  -S " << setw(12) << left << ARG_SEARCH_GREEDY
						<< "Greedy hill climbing algorithm" << endl;
				cerr << "  -S " << setw(12) << left << ARG_SEARCH_GREEDY_EXT
						<< "Extended greedy hill climbing algorithm" << endl;
				cerr << "  -S " << setw(12) << left << ARG_SEARCH_HIERARCHICAL
						<< "Hierarchical clustering algorithm" << endl;
				exit_partest(EX_CONFIG);
			}
			break;
		case ARG_HCLUSTER_REPS:
			for (int i = 0; value[i] != 0; i++)
				if (!isdigit(value[i])) {
					cerr << "[ERROR] \"-r " << value << "\" must be an integer."
							<< endl;
					exit_partest(EX_CONFIG);
				}
			maxSamples = atoi(value);
			break;
		case ARG_NON_STOP:
			non_stop = true;
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
				exit_partest(EX_CONFIG);
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
					exit_partest(EX_CONFIG);
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
		case ARG_OPTIMIZE:
			if (!strcmp(value, ARG_OPTIMIZE_BESTMODEL)) {
				optimize = OPT_SEARCH;
			} else if (!strcmp(value, ARG_OPTIMIZE_GTR)) {
				optimize = OPT_GTR;
			} else {
				cerr << "[ERROR] \"-n " << value
						<< "\" is not a valid optimize mode. Use one of the following:"
						<< endl;
				cerr << "  -O " << setw(16) << left << ARG_OPTIMIZE_BESTMODEL
						<< "\t Perform a model selection on the best partition"
						<< endl;
				cerr << "  -O " << setw(16) << left << ARG_OPTIMIZE_GTR
						<< "\t Optimize only GTR models on the best partition"
						<< endl;
				exit_partest(EX_CONFIG);
			}
			break;
		case ARG_OUTPUT:
			strcpy(output_dir, value);
			break;
		case ARG_CONFIG_FILE:
			strcpy(config_file, value);
			break;
		case ARG_CONFIG_HELP:
			ConfigParser::printFormat();
			exit_partest(0);
			break;
		case ARG_CONFIG_TEMPLATE:
			ConfigParser::createTemplate();
			exit_partest(0);
			break;
		case ARG_NUM_PROCS:
#ifdef PTHREADS
			number_of_threads = atoi(value);
#else
			cerr
					<< "[ERROR] PThreads version is not available. You must recompile with PTHREADS flag."
					<< endl;
			exit_partest(EX_CONFIG);
#endif
			break;
		default:
			cerr << "[ERROR] \"" << argument
					<< "\" option not recognized.\n Run with --help for more information."
					<< endl;
			exit_partest(EX_CONFIG);
			break;
		}
	}

	bitMask do_rate = RateVarM;
	if (do_f)
		do_rate |= (RateVarF);
	if (do_i)
		do_rate |= (RateVarI);
	if (do_g)
		do_rate |= (RateVarG);

	if (strcmp(input_file, ""))
		ptest->setInputFile(input_file);
	if (strcmp(config_file, ""))
		ptest->setConfigFile(config_file);
	if (strcmp(user_tree, ""))
		ptest->setUserTree(user_tree);
	if (strcmp(output_dir, ""))
		ptest->setOutputDir(output_dir);
	ptest->setDataType(data_type);
	ptest->setDoRate(do_rate);
	ptest->setStartingTopology(startingTopology);
	ptest->setMaxSamples(maxSamples);
	ptest->setOptimize(optimize);
	ptest->setSearchAlgo(searchAlgo);

}

} /* namespace partest */
