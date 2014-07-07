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
 * @file GlobalDefs.cpp
 * @author Diego Darriba
 */

#include "GlobalDefs.h"
#include "Utilities.h"
#include "PrintMeta.h"

#include <stdlib.h>
#include <time.h>
#include <sstream>
#include <iomanip>

using namespace std;

namespace partest {

#ifdef HAVE_MPI
	int myRank = 0;
	int numProcs = 0;
#endif

	/* checkpointing */
	bool ckpAvailable = true;
	string ckpPath;
	string ckpStartingTree = "starting_tree";
	string ckpFinalTree = "final_tree";

	string ** singleGeneNames;
	char * starting_tree = 0;

	/* configuration */
	int number_of_threads = 1;
	DataType data_type = DT_DEFAULT;
	bitMask do_rate = 0;
	StartTopo starting_topology;
	SearchAlgo search_algo;
	int max_samples;
	InformationCriterion ic_type;
	OptimizeMode optimize_mode;
	bool non_stop = false;
	bool compute_final_tree = false;
	bitMask protModels = Utilities::binaryPow(max(NUC_MATRIX_SIZE,PROT_MATRIX_SIZE)) - 1;

	/* input/output */
	string * input_file = 0;
	string * config_file = 0;
	string * user_tree = 0;
	string * output_dir = 0;
	string * models_logfile = 0;
	string * schemes_logfile = 0;
	string * results_logfile = 0;
	bool force_overriding = false;
	bool outputAvailable = true;

	/* data description */
	size_t num_taxa = 0;
	size_t seq_len = 0;
	size_t num_patterns = 0;
	size_t number_of_models = 0;
	size_t number_of_genes = 0;
	size_t number_of_schemes = 0;
	std::vector<t_partitioningScheme> * schemes = 0;

	/* data structures */
	pllQueue * pllPartsQueue = 0;
	partitionList * pllPartitions = 0;
	pllAlignmentData * phylip = 0;
	pllInstance * tree = 0;

	void exit_partest(int status) {
		/* free global variables */
		for (size_t i=0; i<number_of_genes; i++) {
			delete singleGeneNames[i];
		}
		free(singleGeneNames);

		if (user_tree)
			delete (user_tree);
		if (config_file)
			delete (config_file);
		if (input_file)
			delete (input_file);
		if (output_dir)
			delete (output_dir);
		if (models_logfile)
			delete (models_logfile);
		if (schemes_logfile)
			delete (schemes_logfile);
		if (results_logfile)
			delete (results_logfile);

		/* exit */
		switch(status) {
		case EX_USAGE:
			PrintMeta::print_usage(cout);
			break;
		case EX_SOFTWARE:
			cerr << " ... internal error that should NEVER raise ..." << endl;
			break;
		case EX_UNAVAILABLE:
			cerr << " ... attempting to use an unavailable feature ..." << endl;
			break;
		case EX_CONFIG:
			cerr << "Run with --help for more information." << endl;
			break;
		case EX_OK:
			/* print nothing */
			break;
		default:
			cerr << "... finishing with error status." << endl;
			break;
		}

		exit(status);
	}

	time_t start_time = time(NULL);

	string timestamp() {
		time_t seconds = time(NULL) - start_time;
		int hours = seconds / 3600;
		seconds %= 3600;
		int minutes = seconds / 60;
		seconds %= 60;
		stringstream ss;
		ss << setw(2) << setfill('0') << hours << ":" << setw(2) << minutes << ":" << setw(2) << seconds;
		return ss.str();
	}
}
