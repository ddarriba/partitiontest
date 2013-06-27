/*
 * PLLModelOptimize.cpp
 *
 *  Created on: Jan 8, 2013
 *      Author: diego
 */

#include "PLLModelOptimize.h"
#include <iostream>
#if ! (defined(__ppc) || defined(__powerpc__) || defined(PPC))
#include <xmmintrin.h>
#endif

extern "C" {
#include "globalVariables.h"
#include "parser/phylip/phylip.h"
#include "utils.h"
#include "parser/partition/part.h"
}

namespace partest {

PLLModelOptimize::PLLModelOptimize(ParTestOptions * options) :
		ModelOptimize(options) {
	cout << "CREATING " << endl;
	alignment = static_cast<PLLAlignment *>(options->getAlignment());
	tr = alignment->getTree();

//  tr = options->getAlignment()->getTree();
//  empiricalFrequencies = (double **)malloc(sizeof(double *));
//  empiricalFrequencies[0] = (double *)malloc(4*sizeof(double));
//  empiricalFrequencies[0][0] = empiricalFrequencies[0][1] = empiricalFrequencies[0][2] = empiricalFrequencies[0][3] = 0.25;
//  cout << tr << " " << empiricalFrequencies << endl;
// // read_phylip_msa (tr[0], (char *) dataFileName.c_str(), PHYLIP|  _SEQUENTIAL, 0);
//
//#ifndef _MOCK_COMPUTATION
////  initModel(tr, empiricalFrequencies);
//    makeRandomTree(tr);
////  evaluateGeneric(tr, tr->start, TRUE);
//#endif
//  ModelSet modelset(options->getRateVariation(), options->getDataType(), options->getAlignment()->getNumSeqs());
//  cout << tr << " " << empiricalFrequencies << endl;
//  for (int i=0; i<modelset.getNumberOfModels(); i++) optimizeModel(modelset.getModel(i), i);
}

PLLModelOptimize::~PLLModelOptimize() {
}

int PLLModelOptimize::optimizePartitioningScheme(PartitioningScheme * scheme,
		bool forceRecomputation, int current_index, int max_index) {

	cout << "PLL OPTIMIZE STAGE 1" << endl;
	partitionList * partitions = alignment->getPartitions();
	struct pllPhylip * phylip = alignment->getPhylip();

	pllPhylipRemoveDuplicate(phylip, partitions);

	/* Set the topology of the PLL tree from a parsed newick tree */
	//pllTreeInitTopologyNewick (tr, newick, PLL_TRUE);
	/* Or instead of the previous function use the next commented line to create
	 a random tree topology
	 pllTreeInitTopologyRandom (tr, phylip->nTaxa, phylip->label); */

	pllTreeInitTopologyForAlignment(tr, phylip);
	cout << "PLL OPTIMIZE STAGE 2" << endl;
	/* Connect the alignment with the tree structure */
	if (!pllLoadAlignment(tr, phylip, partitions, PLL_DEEP_COPY)) {
		cerr << "ERROR: Incompatible tree/alignment combination" << endl;
		Utilities::exit_partest(EX_SOFTWARE);
	}

	/* Initialize the model TODO: Put the parameters in a logical order and change the TRUE to flags */
	pllInitModel(tr, PLL_TRUE, phylip, partitions);

	//pllPhylipDestroy(phylip);
	cout << "PLL OPTIMIZE STAGE 3" << endl;
	pllComputeRandomizedStepwiseAdditionParsimonyTree(tr, partitions);
	Tree2String(tr->tree_string, tr, partitions, tr->start->back, PLL_TRUE,
			PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH,
			PLL_FALSE, PLL_FALSE);
	printf("Tree: %s %d\n", tr->tree_string, tr->start->number);
	evaluateGeneric(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);

	/* now start the ML search algorithm */
	analdef *adef = (analdef*) rax_calloc(1, sizeof(analdef));
	adef->max_rearrange = 21;
	adef->stepwidth = 5;
	adef->initial = 10;
	adef->bestTrav = 10;
	adef->initialSet = PLL_FALSE;
	adef->mode = BIG_RAPID_MODE;
	adef->likelihoodEpsilon = 0.1;
	adef->permuteTreeoptimize = PLL_FALSE;
	adef->perGeneBranchLengths = PLL_FALSE;
	adef->useCheckpoint = PLL_FALSE;

	strcpy(resultFileName, "pll-result.out");
	strcpy(infoFileName, "pll-info.out");
	strcpy(logFileName, "pll-log.out");

	cout << "PARSED " << partitions->numberOfPartitions << " PARTITIONS:"
			<< endl;
	for (int i = 0; i < partitions->numberOfPartitions; i++) {
		cout << i + 1 << ": " << partitions->partitionData[i]->partitionName
				<< endl;
	}

	//  treeEvaluate(tr, 32);
	evaluate(tr, partitions, adef, PLL_TRUE);
	printf("tree evaluated: %f\n", tr->likelihood);
	evaluate(tr, partitions, adef, PLL_TRUE);
	printf("tree evaluated: %f\n", tr->likelihood);
	evaluate(tr, partitions, adef, PLL_TRUE);
	printf("tree evaluated: %f\n", tr->likelihood);

	/* fill models */
	for (int i = 0; i < scheme->getNumberOfElements(); i++) {
		PartitionElement * element = scheme->getElement(i);
		Model * model = element->getModelset()->getModel(0);
		model->setLnL(partitions->partitionData[i]->partitionLH);
		model->setFrequencies(partitions->partitionData[i]->frequencies);
		model->setAlpha(partitions->partitionData[i]->alpha);
		model->setRates(partitions->partitionData[i]->substRates);
		model->print(cout);
	}

	return 0;
}

int PLLModelOptimize::optimizeModel(Model * model,
		PartitionElement * partitionElement, int index, int groupCount) {

#ifdef _MOCK_COMPUTATION
	model->setLnL(-1 * (rand()%500 + 10000));
	if (model->isPInv()) {
		model->setpInv(0.5d);
	}
	if (model->isGamma()) {
		model->setAlpha(0.5d);
	}
	if (model->isPF()) {
		double *freqs = (double *)malloc(model->getNumberOfFrequencies() * sizeof(double));
		int i;
		for (i=0; i < model->getNumberOfFrequencies(); i++) {
			freqs[i] = 1.0/model->getNumberOfFrequencies();
		}
		model->setFrequencies(freqs);
		free(freqs);
	}
#else
	cerr << "ERROR: Only full schemes should be optimized with PLL" << endl;
	Utilities::exit_partest(EX_SOFTWARE);
#endif
}

} /* namespace partest */
