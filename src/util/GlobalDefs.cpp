/*
 * GlobalDefs.cpp
 *
 *  Created on: Apr 8, 2014
 *      Author: diego
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

#ifdef _MPI
	int myRank = 0;
	int numProcs = 0;
#endif

	int number_of_threads = 1;
	bool ckpAvailable = false;
	string ckpPath;
	string ckpStartingTree = "starting_tree";
	string ckpFinalTree = "final_tree";
	string ** singleGeneNames;
	char * starting_tree = 0;

	DataType data_type;
	bitMask do_rate;
	StartTopo starting_topology;
	SearchAlgo search_algo;
	int max_samples;
	InformationCriterion ic_type;
	OptimizeMode optimize_mode;
	bool non_stop = false;
	bool compute_final_tree = false;

	string * input_file = 0;
	string * config_file = 0;
	string * user_tree = 0;
	string * output_dir = 0;
	string * models_logfile = 0;
	string * schemes_logfile = 0;
	string * results_logfile = 0;

	bitMask protModels = Utilities::binaryPow(max(NUC_MATRIX_SIZE,PROT_MATRIX_SIZE)) - 1;

	pllQueue * pllPartsQueue = 0;
	partitionList * pllPartitions = 0;
	pllAlignmentData * phylip = 0;
	pllInstance * tree = 0;

	size_t num_taxa = 0;
	size_t seq_len = 0;
	size_t num_patterns = 0;
	size_t number_of_models = 0;
	size_t number_of_genes = 0;
	size_t number_of_schemes = 0;
	std::vector<t_partitioningScheme> * schemes = 0;

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
		if (status == EX_USAGE) {
			PrintMeta::print_usage(cout);
		}
		if (status == EX_SOFTWARE) {
			cerr << " ... internal error that should NEVER raise ..." << endl;
		}
		if (status == EX_UNAVAILABLE) {
			cerr << " ... attempting to use an unavailable feature ..." << endl;
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
