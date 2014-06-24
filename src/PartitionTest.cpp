/*
 * PartitionTest.cpp
 *
 *  Created on: Apr 8, 2014
 *      Author: diego
 */

#include "PartitionTest.h"
#include "search/SearchAlgorithm.h"
#include "search/HierarchicalClusteringSearchAlgorithm.h"
#include "search/GreedySearchAlgorithm.h"
#include "search/RandomSearchAlgorithm.h"
#include "indata/PartitioningScheme.h"
#include "indata/PartitionMap.h"
#include "util/PrintMeta.h"
#include "util/Utilities.h"
#include "parser/ArgumentParser.h"
#include "parser/ConfigParser.h"

#include <sys/stat.h>
#include <sys/types.h>
#include <fstream>
#include <string.h>

namespace partest {

bool PartitionTest::checkParameters(void) {
	/* check all non-default parameters have been set */
	return (data_type != DT_DEFAULT) & (search_algo != SearchDefault)
			& (optimize_mode != OPT_DEFAULT) & (ic_type != IC_DEFAULT)
			& (input_file != 0);
}

PartitionTest::PartitionTest() {
	phylip = 0;
	pllPartitions = 0;
	tree = 0;
	pllPartsQueue = 0;

	data_type = DT_DEFAULT;
	do_rate = 0;
	starting_topology = DEFAULT_STARTING_TOPOLOGY;
	search_algo = SearchDefault;
	max_samples = 1;
	optimize_mode = OPT_DEFAULT;
	ic_type = IC_DEFAULT;
}

PartitionTest::~PartitionTest() {
	cout << timestamp() << " Execution done." << endl;
	if (phylip) {
		pllAlignmentDataDestroy(phylip);
	}
	if (pllPartsQueue)
		pllQueuePartitionsDestroy(&pllPartsQueue);
	if (tree) {
		if (pllPartitions && pllPartitions->alphaList) {
//			pllPartitionsDestroy(tree, &pllPartitions);
		}
		if (tree->nameHash) {
			pllDestroyInstance(tree);
			tree = 0;
		}
	}

}

bool PartitionTest::configure(void) {

	if (!config_file)
		config_file = new string("");

	ConfigParser parser(config_file->c_str());

	if (optimize_mode == OPT_DEFAULT) {
		optimize_mode = parser.getOptimizeMode();
		if (optimize_mode == OPT_DEFAULT) {
			optimize_mode = DEFAULT_OPTIMIZE;
		}
	}

	if (search_algo == SearchDefault) {
		search_algo = DEFAULT_SEARCH_ALGO;
	}

	// check required arguments
	if (!input_file) {
		if (!strlen(parser.getInputFile())) {
			cerr << "ERROR! Input File (-i) is required!" << endl;
			exit_partest(EX_CONFIG);
		}
		input_file = new string(parser.getInputFile());
	}

	if (parser.getOutputBasePath().length() > 0) {
		output_dir = new string(parser.getOutputBasePath());
	} else {
		const char *baseName = basename(input_file->c_str());
		output_dir = new string("partest_");
		output_dir->append(baseName);
		output_dir->append(os_separator);
	}
	models_logfile = new string(*output_dir + "models");
	schemes_logfile = new string(*output_dir + "schemes");
	results_logfile = new string(*output_dir + "results");

	ckpPath = (*output_dir) + CKP_DIR;
	if (I_AM_ROOT) {
		int resMkdir = mkdir(output_dir->c_str(), 0777);
		if (resMkdir) {
			if (errno == EEXIST) {
				cerr << "[WARNING] Output directory " << (*output_dir)
						<< " already exists. Output files might be overwritten."
						<< endl;
			} else {
				cerr << "[WARNING] ***** WARNING *****" << endl;
				cerr << "[WARNING] Output directory " << (*output_dir)
						<< " cannot be created. No output files will be stored."
						<< endl;
				cerr << "[WARNING] ***** WARNING *****" << endl;
			}
		}

		if (!resMkdir || errno == EEXIST) {
			mkdir(ckpPath.c_str(), 0777);
			ckpAvailable = true;
		}
	}

#ifdef _MPI
	int tmpInt = ckpAvailable;
	MPI_Bcast(&tmpInt, 1, MPI_INT, 0, MPI_COMM_WORLD );
	ckpAvailable = tmpInt;
#endif

	pllInstanceAttr attr;
	attr.fastScaling = PLL_FALSE;
	attr.randomNumberSeed = 12345;
	attr.rateHetModel = PLL_GAMMA;
	attr.saveMemory = PLL_FALSE;
	attr.useRecom = PLL_FALSE;
	attr.numberOfThreads = number_of_threads;

	if (data_type == DT_DEFAULT) {
		data_type = parser.getDataType();
		if (data_type == DT_DEFAULT) {
			data_type = DEFAULT_DATA_TYPE;
		}
	}

	switch (optimize_mode) {
	case OPT_SEARCH:
		number_of_models =
				data_type == DT_NUCLEIC ? NUC_MATRIX_SIZE : PROT_MATRIX_SIZE;
		if (do_rate & RateVarF) {
			number_of_models *= 2;
		}
		break;
	case OPT_GTR:
		number_of_models = 1;
		break;
	case OPT_CUSTOM:
		number_of_models = Utilities::setbitsCount(protModels);
		if (do_rate & RateVarF) {
			number_of_models *= 2;
		}
		break;
	case OPT_DEFAULT:
		exit_partest(EX_SOFTWARE);
	}

	tree = pllCreateInstance(&attr);
	phylip = pllParseAlignmentFile(PLL_FORMAT_PHYLIP, input_file->c_str());
	pllPartitions = pllPartitionsCommit(pllPartsQueue, phylip);

	num_taxa = phylip->sequenceCount;
	seq_len = phylip->sequenceLength;

	num_patterns = phylip->sequenceLength;

	return EX_OK;
}

void PartitionTest::setDataType(DataType dataType) {
	data_type = dataType;
}

void PartitionTest::setDoRate(bitMask doRate) {
	do_rate = doRate;
}

void PartitionTest::setIcType(InformationCriterion icType) {
	ic_type = icType;
}

void PartitionTest::setMaxSamples(int maxSamples) {
	max_samples = maxSamples;
}

void PartitionTest::setOptimize(OptimizeMode optimize) {
	optimize_mode = optimize;
}

void PartitionTest::setSearchAlgo(SearchAlgo searchAlgo) {
	search_algo = searchAlgo;
}

void PartitionTest::setStartingTopology(StartTopo startingTopology) {
	starting_topology = startingTopology;
}

void PartitionTest::setConfigFile(const char * configFile) {
	if (config_file)
		exit_partest(EX_IOERR);
	config_file = new string(configFile);
}

void PartitionTest::setInputFile(const char * inputFile) {
	if (input_file)
		exit_partest(EX_IOERR);
	input_file = new string(inputFile);
}

void PartitionTest::setOutputDir(const char * outputDir) {
	if (output_dir)
		exit_partest(EX_IOERR);
	output_dir = new string(outputDir);
}

void PartitionTest::setUserTree(const char * userTree) {
	if (user_tree)
		exit_partest(EX_IOERR);
	user_tree = new string(userTree);
}

}

using namespace partest;
using namespace std;

int main(int argc, char * argv[]) {

#ifdef _MPI
	if (MPI_Init( &argc, &argv )) {
		cerr << "Error initializing MPI!!" << endl;
		exit_partest(EX_PROTOCOL);
	}

	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	MPI_Comm_rank( MPI_COMM_WORLD, &myRank );
#endif

	PartitionTest * ptest = new PartitionTest();

	if (I_AM_ROOT) {
		PrintMeta::print_header(cout);
	}

	ArgumentParser * parser = new ArgumentParser(ptest);
	parser->parse(argc, argv);

	ptest->configure();
	if (!config_file && !ptest->checkParameters()) {
		cerr << "ERROR: Configuration parameters are missing" << endl;
		exit_partest(EX_CONFIG);
	}

	if (I_AM_ROOT) {
		PrintMeta::print_options(cout);
	}

	SearchAlgorithm * searchAlgo = 0;
	switch (search_algo) {
	case SearchHCluster:
		searchAlgo = new HierarchicalClusteringSearchAlgorithm();
		break;
	case SearchGreedy:
		searchAlgo = new GreedySearchAlgorithm();
		break;
	case SearchRandom:
		searchAlgo = new RandomSearchAlgorithm();
		break;
	case SearchExhaustive:
	case SearchGreedyExtended:
	default:
		break;
	}

	PartitioningScheme * bestScheme = searchAlgo->start();

	if (I_AM_ROOT) {
		ModelOptimize mo;
		mo.buildFinalTree(bestScheme, true);

		if (ckpAvailable && results_logfile) {
			ofstream ofs(results_logfile->c_str(), ios::out);
			PrintMeta::print_results_xml(ofs, bestScheme);
			ofs.close();
		}
		PrintMeta::print_results(cout, bestScheme);
	}

	delete bestScheme;
	delete searchAlgo;

	PartitionMap::deleteInstance();

	delete parser;
	delete ptest;

#ifdef _MPI
	MPI_Finalize();
#endif

	exit_partest(EX_OK);
}
