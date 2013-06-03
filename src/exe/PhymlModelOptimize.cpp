/*
 * PhymlModelOptimize.cpp
 *
 *  Created on: Jan 10, 2013
 *      Author: diego
 */

#include "PhymlModelOptimize.h"
#include "model/ModelSet.h"
#include <time.h>
#include <string.h>
#include <alloca.h>

/* external C */
#ifndef PHYML_H
#define PHYML_H
extern "C" {
#include "partest_to_phyml.h"
//#include "utilities.h"
// void int compute_lk(int argc, char **argv);
typedef struct __Option { /*! mostly used in 'help.c' */
	struct __Model *mod; /*! pointer to a substitution model */
	struct __Tree *tree; /*! pointer to the current tree */
	struct __Align **data; /*! pointer to the uncompressed sequences */
	struct __Tree *cstr_tree; /*! pointer to a constraint tree (can be a multifurcating one) */
	struct __Calign *cdata; /*! pointer to the compressed sequences */
	struct __Super_Tree *st; /*! pointer to supertree */
	struct __Tnexcom **nex_com_list;
	struct __List_Tree *treelist; /*! list of trees. */

	int interleaved; /*! interleaved or sequential sequence file format ? */
	int in_tree; /*! =1 iff a user input tree is used as input */

	char *in_align_file; /*! alignment file name */
	FILE *fp_in_align; /*! pointer to the alignment file */

	char *in_tree_file; /*! input tree file name */
	FILE *fp_in_tree; /*! pointer to the input tree file */

	char *in_constraint_tree_file; /*! input constraint tree file name */
	FILE *fp_in_constraint_tree; /*! pointer to the input constraint tree file */

	char *out_tree_file; /*! name of the tree file */
	FILE *fp_out_tree;

	char *out_trees_file; /*! name of the tree file */
	FILE *fp_out_trees; /*! pointer to the tree file containing all the trees estimated using random starting trees */

	char *out_boot_tree_file; /*! name of the tree file */
	FILE *fp_out_boot_tree; /*! pointer to the bootstrap tree file */

	char *out_boot_stats_file; /*! name of the tree file */
	FILE *fp_out_boot_stats; /*! pointer to the statistics file */

	char *out_stats_file; /*! name of the statistics file */
	FILE *fp_out_stats;

	char *out_trace_file; /*! name of the file in which the likelihood of the model is written */
	FILE *fp_out_trace;

	char *out_lk_file; /*! name of the file in which the likelihood of the model is written */
	FILE *fp_out_lk;

	char *out_ps_file; /*! name of the file in which tree(s) is(are) written */
	FILE *fp_out_ps;

	char *aa_rate_mat_file;
	FILE *fp_aa_rate_mat;

	char *clade_list_file;

	int datatype; /*! 0->DNA, 1->AA */
	int print_boot_trees; /*! =1 if the bootstrapped trees are printed in output */
	int out_stats_file_open_mode; /*! opening file mode for statistics file */
	int out_tree_file_open_mode; /*! opening file mode for tree file */
	int n_data_sets; /*! number of data sets to be analysed */
	int n_trees; /*! number of trees */
	int init_len; /*! sequence length */
	int n_otu; /*! number of taxa */
	int n_data_set_asked; /*! number of bootstrap replicates */
	char *nt_or_cd; /*! nucleotide or codon data ? (not used) */
	int multigene; /*! if=1 -> analyse several partitions. */
	int config_multigene;
	int n_part; /*! number of data partitions */
	int curr_gt;
	int ratio_test; /*! from 1 to 4 for specific branch supports, 0 of not */
	int ready_to_go;
	int data_file_format; /*! Data format: Phylip or Nexus */
	int tree_file_format; /*! Tree format: Phylip or Nexus */

	int curr_interface;
	int r_seed; /*! random seed */
	int collapse_boot; /*! 0 -> branch length on bootstrap trees are not collapsed if too small */
	int random_boot_seq_order; /*! !0 -> sequence order in bootstrapped data set is random */
	int print_trace;
	int print_site_lnl;
	int m4_model;
	int rm_ambigu; /*! 0 is the default. 1: columns with ambiguous characters are discarded prior further analysis */
	int colalias;
	int appebr_run_ID;
	char *run_id_string;
	int quiet; /*! 0 is the default. 1: no interactive question (for batch mode) */
	int lk_approx; /* EXACT or NORMAL */
	char **alphabet;
	int codpos;
	int mutmap;

	char **long_tax_names;
	char **short_tax_names;
	int size_tax_names;

	phydbl *z_scores;
	phydbl *lat;
	phydbl *lon;

	int boot_prog_every;

	int mem_question;
	int do_alias_subpatt;
	struct __Tmcmc *mcmc;
	struct __T_Rate *rates;

};
}
#endif

namespace partest {

using namespace std;

PhymlModelOptimize::PhymlModelOptimize(ParTestOptions * options) :
		ModelOptimize(options) {

	// NOTHING
}

PhymlModelOptimize::~PhymlModelOptimize() {
	// NOTHING
}

#ifdef EXPERIMENTAL
void computeDistances(void) {
	for (int i = 1; i < 14; i++) {
		for (int j = 0; j < i; j++) {
			ProtMatrix m1 = static_cast<ProtMatrix>(i);
			ProtMatrix m2 = static_cast<ProtMatrix>(j);
			printf("%6.2f\t", Utilities::getEuclideanDistance(m1,m2));
		}
		printf("\n");
	}
	compute_distances();
}
#endif

int PhymlModelOptimize::optimizeModel(Model * model,
		PartitionElement * partitionElement, int index, int groupCount) {

	int dataType = 0, freqType = 0;
	int optimizeTopo = 0, optimizeBranchLengths = 0, optimizeRates = 0;

	double * frequencies = 0;
	double * rates = 0;

#ifdef BUILD_PHYML_ARGS
	string arguments = "-b 0 -i ";
	arguments += options->getInputFile();
#endif

	switch (options->getDataType()) {
	case DT_NUCLEIC:
		rates = (double *) alloca(6 * sizeof(double));
		frequencies = (double *) alloca(4 * sizeof(double));
		dataType = DATATYPE_NT;
#ifdef BUILD_PHYML_ARGS
		arguments += " -d nt";
		!(model->isPF()) ? (arguments += " -f m") : (arguments +=
				" -f 0.25,0.25,0.25,0.25");
#endif
		freqType = model->isPF() ? FREQTYPE_MODEL : FREQTYPE_CUSTOM;
		if (freqType == FREQTYPE_CUSTOM) {
			frequencies[0] = frequencies[1] = frequencies[2] = frequencies[3] =
					0.25;
		}
		break;
	case DT_PROTEIC:
		rates = 0;
		frequencies = (double *) alloca(20 * sizeof(double));
		dataType = DATATYPE_AA;
#ifdef BUILD_PHYML_ARGS
		arguments += " -d aa";
		(model->isPF()) ? (arguments += " -f e") : (arguments += " -f m");
#endif
		freqType = !(model->isPF()) ? FREQTYPE_MODEL : FREQTYPE_EMPIRICAL;
		break;
	}

	switch (options->getStartingTopology()) {
	case StartTopoBIONJ:
		optimizeTopo = false;
		optimizeBranchLengths = true;
		//optimizeRates = true;
#ifdef BUILD_PHYML_ARGS
		arguments += " -o lr";
#endif
		break;
	case StartTopoML:
		optimizeTopo = true;
		optimizeBranchLengths = true;
		//optimizeRates = true;
#ifdef BUILD_PHYML_ARGS
		arguments += " -o tlr";
#endif
		break;
	case StartTopoFIXED:
		//TODO: TEMPORARY THIS IS A REGULAR BIONJ!
		optimizeTopo = false;
		optimizeBranchLengths = false;
		//optimizeRates = true;
#ifdef BUILD_PHYML_ARGS
		arguments += " -o lr";
#endif
		break;
	case StartTopoUSER:
		optimizeTopo = false;
		optimizeBranchLengths = false;
		//optimizeRates = false;
#ifdef BUILD_PHYML_ARGS
		arguments += " -o lr -u ";
		arguments += options->getTreeFile();
#endif
		break;
	}
	optimizeRates = dataType == DATATYPE_NT;

#ifdef BUILD_PHYML_ARGS
	(model->isPInv()) ? (arguments += " -v e") : (arguments += " -v 0.0");
	(model->isGamma()) ? (arguments += " -a e") : (arguments += " -a 300 -c 1");
	arguments += " -m " + model->getMatrixName();

	char argsCStr[arguments.length()];
	strcpy(argsCStr, arguments.c_str());

	int numWords = Utilities::countWords(arguments);
	char *params[numWords];
	int curTok = 0;
	char* token = strtok(argsCStr, " ");

	while (token) {
		params[++curTok] = token;
		token = strtok(NULL, " ");
	}

#ifdef DEBUG
	cout << "[DEBUG] PhyML arguments: " << arguments << endl;
#endif
#endif

	notify_observers(MT_SINGLE_INIT, index, model, time(NULL), index,
			groupCount);

	phyml_indata indata;
	indata.ioFile = options->getInputFile().c_str();
	indata.treeFile =
			options->getStartingTopology() == StartTopoUSER ?
					options->getTreeFile() : 0;
	indata.dataType = dataType;
	indata.freqType = freqType;
	indata.frequencies = frequencies;
	indata.optimizeTopology = optimizeTopo;
	indata.optimizeBranchLengths = optimizeBranchLengths;
	indata.optimizeRates = optimizeRates;
	indata.optimizePInvar = model->isPInv();
	indata.optimizeGamma = model->isGamma();
	indata.model = model->getMatrixName().c_str();

	phyml_outdata * outdata = (phyml_outdata *) alloca(sizeof(phyml_outdata));
	outdata->frequencies = frequencies;
	outdata->rates = rates;

	struct __Option *io = build_options(indata);
	io->cdata =
			((PhymlAlignment *) partitionElement->getAlignment())->getCData();

	double lnL = phyml_lk(io, outdata);
	model->setLnL(lnL);

	model->setTree(outdata->tree);
	free(outdata->tree);

	if (model->isPInv()) {
		model->setpInv(outdata->pinv);
	}
	if (model->isGamma()) {
		model->setAlpha(outdata->alpha);
	}
	model->setFrequencies(outdata->frequencies);
	model->setRates(outdata->rates);

//	pt_free_io(io);

	notify_observers(MT_SINGLE_END, index, model, time(NULL), index,
			groupCount);

	return 0;
}

} /* namespace partest */
