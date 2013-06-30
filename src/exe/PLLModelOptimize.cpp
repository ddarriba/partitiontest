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

}

PLLModelOptimize::~PLLModelOptimize() {
}

int PLLModelOptimize::optimizePartitioningScheme(PartitioningScheme * scheme,
		bool forceRecomputation, int current_index, int max_index) {

	for (int i = 0; i < scheme->getNumberOfElements(); i++) {
		PartitionElement * element = scheme->getElement(i);
		if (!element->getBestModel()) {
			optimizePartitionElement(element, i + 1,
					scheme->getNumberOfElements());
		}
	}
	return 0;
	partitionList * partitions = alignment->getPartitions();
	partitionList * nextPartitions;

	/* build new partitioning scheme */
	if (scheme->getNumberOfElements() == partitions->numberOfPartitions) {
		nextPartitions = partitions;
	} else {

	}

//	cout << partitions->partitionData[0]->partition << endl;
	if (tr->likelihood == 0) {

		//pllPhylipDestroy(phylip);

		pllComputeRandomizedStepwiseAdditionParsimonyTree(tr, partitions);
		Tree2String(tr->tree_string, tr, partitions, tr->start->back, PLL_TRUE,
				PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH,
				PLL_FALSE, PLL_FALSE);
		evaluateGeneric(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
	}
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

	cout << "PLL OPTIMIZE STAGE 1" << endl;
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
		if (model->isGamma())
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

	PLLAlignment * alignment =
			static_cast<PLLAlignment *>(partitionElement->getAlignment());
	pllInstance * tree = alignment->getTree();
	partitionList * partitions = alignment->getPartitions();
	//pllPhylipDestroy(phylip);
	pllComputeRandomizedStepwiseAdditionParsimonyTree(tree, partitions);
//	Tree2String(tree->tree_string, tree, partitions, tree->start->back, PLL_TRUE,
//			PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH,
//			PLL_FALSE, PLL_FALSE);
//	printf("Tree: %s\n", tree->tree_string);

	strcpy(resultFileName, "pll-result.out");
	strcpy(infoFileName, "pll-info.out");
	strcpy(logFileName, "pll-log.out");

	evaluateGeneric(tree, partitions, tree->start, PLL_TRUE, PLL_FALSE);
	//  treeEvaluate(tr, 32);
	evaluate(tree, partitions, adef, PLL_TRUE);
	model->setLnL(tree->likelihood);
	model->setFrequencies(partitions->partitionData[0]->frequencies);
	if (model->isGamma())
		model->setAlpha(partitions->partitionData[0]->alpha);
	model->setRates(partitions->partitionData[0]->substRates);

	return 0;
	cerr << "ERROR: Only full schemes should be optimized with PLL" << endl;
	Utilities::exit_partest(EX_SOFTWARE);
#endif
}

} /* namespace partest */
