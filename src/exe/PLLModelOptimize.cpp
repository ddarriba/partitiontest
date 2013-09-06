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
#include "parser/newick/newick.h"
#include "parser/partition/part.h"

extern double treeOptimizeRapid(pllInstance *tr, partitionList *pr, int mintrav,
		int maxtrav, analdef *adef, bestlist *bt, infoList *iList);
}

namespace partest {

void PLLModelOptimize::initializeStructs(pllInstance * tree,
		partitionList * partitions, pllAlignmentData * phylip) {

	pllTreeInitTopologyForAlignment(tree, phylip);

	/* Connect the alignment with the tree structure */
	if (!pllLoadAlignment(tree, phylip, partitions, PLL_SHALLOW_COPY)) {
		cerr << "ERROR: Incompatible tree/alignment combination" << endl;
		Utilities::exit_partest(EX_SOFTWARE);
	}

	/* Initialize the model TODO: Put the parameters in a logical order and change the TRUE to flags */
	pllInitModel(tree, partitions, phylip);

}

double PLLModelOptimize::evaluateNNI(pllInstance * tr, partitionList *pr,
		bool estimateModel) {

#ifdef DEBUG
	cout << "[TRACE] START EVALUATING NNI WITH PLL" << endl;
#endif
	if (pllNniSearch(tr, pr, estimateModel)) {
#ifdef DEBUG
		cout << "[TRACE] DONE EVALUATING NNI WITH PLL " << endl;
		cout << "[TRACE] LIKELIHOOD IS " << tr->likelihood << endl;
#endif
		return tr->likelihood;
	} else {
#ifdef DEBUG
		cout << "[TRACE] ERROR EVALUATING NNI WITH PLL" << endl;
		exit(-1);
#endif
		return 0;
	}

}

double PLLModelOptimize::evaluateSPR(pllInstance * tr,
		partitionList *partitions, bool estimateModel) {

	pllListSPR * bestList;
	int i;

	/* TODO: evaluate likelihood, create interface calls */
	evaluateGeneric(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
	//printf("Likelihood: %f\n\n", tr->likelihood);

	/* another eval*/
	double computed_lh = tr->likelihood;
	evaluateGeneric(tr, partitions, tr->start, PLL_FALSE, PLL_FALSE);
	assert(computed_lh == tr->likelihood);
	int numBranches =
			partitions->perGeneBranchLengths ?
					partitions->numberOfPartitions : 1;

	tr->thoroughInsertion = 1;
	//printf(	"Computing the best 20 SPRs in a radius (1,20)     [thoroughInsertion = enabled]\n");
	bestList = pllComputeSPR(tr, partitions, tr->nodep[tr->mxtips + 1], 1, 20,
			20);

	//printf("Number of SPRs computed : %d\n", bestList->entries);

//	for (i = 0; i < bestList->entries; ++i) {
//		printf("\t bestList->sprInfo[%2d].likelihood     = %f\n", i,
//				bestList->sprInfo[i].likelihood);
//	}

//	printf(
//			"Committing bestList->sprInfo[0]                   [thoroughInsertion = disabled]\n");
	tr->thoroughInsertion = 0;
	pllCommitSPR(tr, partitions, &(bestList->sprInfo[0]), PLL_TRUE);

	evaluateGeneric(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
//	printf("New likelihood: %f\n\n", tr->likelihood);

	tr->thoroughInsertion = PLL_FALSE;
	pllDestroyListSPR(&bestList);
//	printf(
//			"Computing the best 20 SPRs in a radius (1,30)     [thoroughInsertion = enabled]\n");
	bestList = pllComputeSPR(tr, partitions, tr->nodep[tr->mxtips + 1], 1, 30,
			20);

//	printf("Number of SPRs computed : %d\n", bestList->entries);
//	for (i = 0; i < bestList->entries; ++i) {
//		printf("\t bestList->sprInfo[2%d].likelihood     = %f\n", i,
//				bestList->sprInfo[i].likelihood);
//	}

//	printf(
//			"Committing bestList->sprInfo[0]                   [thoroughInsertion = false]\n");
	tr->thoroughInsertion = 0;
	pllCommitSPR(tr, partitions, &(bestList->sprInfo[0]), PLL_TRUE);

	evaluateGeneric(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
//	printf("New likelihood: %f\n\n", tr->likelihood);
//
//	printf("Rolling back...\n");
	pllRollbackSPR(tr, partitions);
	evaluateGeneric(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
//	printf("New likelihood: %f\n\n", tr->likelihood);
//
//	printf("Rolling back...\n");
	pllRollbackSPR(tr, partitions);
	evaluateGeneric(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
//	printf("New likelihood: %f\n\n", tr->likelihood);

	pllDestroyListSPR(&bestList);

	return tr->likelihood;
}

PLLModelOptimize::PLLModelOptimize(ParTestOptions * options) :
		ModelOptimize(options) {
	alignment = static_cast<PLLAlignment *>(options->getAlignment());
	tr = alignment->getTree();
}

PLLModelOptimize::~PLLModelOptimize() {
}

int PLLModelOptimize::optimizePartitioningSchemeAtOnce(
		PartitioningScheme * scheme) {

	/* build the whole alignment */
	char * pllPartitionsFile = strdup("/tmp/tmpfileXXXXXX");
	mkstemp(pllPartitionsFile);

	ofstream pllOutputStream(pllPartitionsFile);

	for (int i = 0; i < scheme->getNumberOfElements(); i++) {
		PartitionElement * element = scheme->getElement(i);
		pllOutputStream << "DNA, PART" << i << "=" << element->getStart(0)
				<< "-" << element->getEnd(0);
		for (int j = 1; j < element->getNumberOfSections(); j++) {
			pllOutputStream << "," << element->getStart(j) << "-"
					<< element->getEnd(j);
		}
		pllOutputStream << endl;
	}
	pllOutputStream.close();

#ifdef DEBUG
	cout << "[TRACE] PLLModelOptimize - Parsing partitions" << endl;
#endif
	struct pllQueue * parts = pllPartitionParse(pllPartitionsFile);

#ifdef DEBUG
	cout << "[TRACE] PLLModelOptimize - Removing temporary file" << endl;
#endif
	if (remove(pllPartitionsFile) != 0)
		cerr << "Error deleting temporary file" << endl;

#ifdef DEBUG
	cout << "[TRACE] PLLModelOptimize - Committing partitions" << endl;
#endif
	partitionList * partitions = pllPartitionsCommit(parts,
			alignment->getPhylip());
	pllQueuePartitionsDestroy(&parts);

#ifdef DEBUG
	cout << "[TRACE] PLLModelOptimize - Creating tree instance" << endl;
#endif
	pllInstanceAttr * attr = (pllInstanceAttr *) rax_malloc(
			sizeof(pllInstanceAttr));
	attr->rateHetModel = GAMMA;
	attr->fastScaling = PLL_FALSE;
	attr->saveMemory = PLL_FALSE;
	attr->useRecom = PLL_FALSE;
	attr->randomNumberSeed = 12345;
	attr->numberOfThreads = 1;
	pllInstance * tr = pllCreateInstance(attr);
	rax_free(attr);

#ifdef DEBUG
	cout << "[TRACE] PLLModelOptimize - Initializing topology" << endl;
#endif
	pllTreeInitTopologyForAlignment(tr, alignment->getPhylip());
	/* Connect the alignment with the tree structure */
#ifdef DEBUG
	cout << "[TRACE] PLLModelOptimize - Connecting tree with alignment" << endl;
#endif
	if (!pllLoadAlignment(tr, alignment->getPhylip(), partitions,
	PLL_SHALLOW_COPY)) {
		cerr << "ERROR: Incompatible tree/alignment combination" << endl;
		Utilities::exit_partest(EX_SOFTWARE);
	}

#ifdef DEBUG
	cout << "[TRACE] PLLModelOptimize - Initializing model" << endl;
#endif
	/* Initialize the model TODO: Put the parameters in a logical order and change the TRUE to flags */
	pllInitModel(tr, partitions, alignment->getPhylip());
#ifdef DEBUG
	cout << "[TRACE] PLLModelOptimize - Initialized model" << endl;
#endif

	/* set best-fit models */
	for (int i = 0; i < scheme->getNumberOfElements(); i++) {
		Model * model = scheme->getElement(i)->getBestModel()->getModel();

		cout << "PART" << i << ": " << model->getName() << " ";

		for (int f = 0; f < 4; f++) {
			cout << model->getFrequencies()[f] << " ";
		}
		for (int r = 0; r < 6; r++) {
			cout << model->getRates()[r] << " ";
		}
		cout << model->getAlpha() << endl;
		pInfo * current_part = partitions->partitionData[i];
		current_part->optimizeBaseFrequencies = PLL_FALSE;
		current_part->alpha = model->getAlpha();
		memcpy(current_part->frequencies, model->getFrequencies(),
				4 * sizeof(double));
		memcpy(current_part->substRates, model->getRates(), 6 * sizeof(double));
		initReversibleGTR(tr, partitions, i);
		makeGammaCats(current_part->alpha, current_part->gammaRates, 4,
				tr->useMedian);
	}

	if (options->getTreeString() == 0) {
		rax_free(tr->ti);
		pllComputeRandomizedStepwiseAdditionParsimonyTree(tr, partitions);
	} else {
		rax_free(tr->nodep);
		rax_free(tr->td[0].ti);
		rax_free(tr->td[0].parameterValues);
		rax_free(tr->td[0].executeModel);
		rax_free(tr->nodeBaseAddress);
		for (int j = 0; j < tr->nameHash->tableSize; ++j) {
			rax_free(tr->nameHash->table[j]);
		}
		rax_free(tr->nameHash);
		tr->nameHash = NULL;
		pllNewickTree * nt = pllNewickParseString(options->getTreeString());
		pllTreeInitTopologyNewick(tr, nt, PLL_TRUE);
	}

	evaluateGeneric(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
	evaluateNNI(tr, partitions, false);

	Tree2String(tr->tree_string, tr, partitions, tr->start->back, PLL_TRUE,
			PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH,
			PLL_FALSE, PLL_FALSE);

	cout << "FINAL TREE: " << tr->tree_string << endl;

	pllPartitionsDestroy(tr, &partitions);

	return 0;
}

int PLLModelOptimize::optimizePartitioningScheme(PartitioningScheme * scheme,
		bool forceRecomputation, int current_index, int max_index) {

	for (int i = 0; i < scheme->getNumberOfElements(); i++) {
		PartitionElement * element = scheme->getElement(i);
		if (!element->getBestModel()) {

			PLLAlignment * alignment =
					static_cast<PLLAlignment *>(element->getAlignment());
			pllInstance * tree = alignment->getTree();
			partitionList * partitions = alignment->getPartitions();
			pllAlignmentData * phylip = alignment->getPhylip();

			initializeStructs(tree, partitions, phylip);

			if (options->getTreeString() == 0) {
				pllComputeRandomizedStepwiseAdditionParsimonyTree(tree,
						partitions);
			} else {
				rax_free(tree->nodep);
				rax_free(tree->td[0].ti);
				rax_free(tree->td[0].parameterValues);
				rax_free(tree->td[0].executeModel);
				rax_free(tree->nodeBaseAddress);
				for (int j = 0; j < tree->nameHash->tableSize; ++j) {
					rax_free(tree->nameHash->table[j]);
				}
				rax_free(tree->nameHash);
				tree->nameHash = NULL;

				pllNewickTree * nt = pllNewickParseString(
						options->getTreeString());
				pllTreeInitTopologyNewick(tree, nt, PLL_FALSE);
				pllNewickParseDestroy(&nt);
			}

			optimizePartitionElement(element, i + 1,
					scheme->getNumberOfElements());
		}
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

	PLLAlignment * alignment =
			static_cast<PLLAlignment *>(partitionElement->getAlignment());
	pllInstance * tree = alignment->getTree();
	partitionList * partitions = alignment->getPartitions();

	pllNewickTree * nt = pllNewickParseString(options->getTreeString());
	pllTreeInitTopologyNewick(tree, nt, PLL_FALSE);
	pllNewickParseDestroy(&nt);
	Tree2String(tree->tree_string, tree, partitions, tree->start->back,
			PLL_TRUE, PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE,
			PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);

	const char * m = model->getMatrixName().c_str();
	char * symmetryPar = (char *) malloc(12 * sizeof(char));
	symmetryPar[0] = m[0];
	symmetryPar[11] = '\0';
	for (int i = 1; i < 6; i++) {
		symmetryPar[(i - 1) * 2 + 1] = ',';
		symmetryPar[i * 2] = m[i];
	}

	pllSetSubstitutionRateMatrixSymmetries(symmetryPar, partitions, 0);

	partitions->partitionData[0]->optimizeBaseFrequencies = model->isPF();
	if (!model->isPF()) {
		partitions->partitionData[0]->optimizeBaseFrequencies = PLL_FALSE;
		for (int i = 0; i < 4; i++) {
			partitions->partitionData[0]->frequencies[i] = 0.25;
		}
		for (int i = 0; i < 6; i++) {
			partitions->partitionData[0]->substRates[i] = 1;
		}
	}

	// TODO: I think this means CAT categories, and not GAMMA
	partitions->partitionData[0]->alpha = 100;
	if (model->isGamma()) {
		partitions->partitionData[0]->numberOfCategories = 4;
		tree->categories = 4;
		tree->maxCategories = 4;
	} else {
		partitions->partitionData[0]->numberOfCategories = 1;
		tree->categories = 1;
		tree->maxCategories = 1;
	}

	free(symmetryPar);

	initReversibleGTR(tree, partitions, 0);
	evaluateGeneric(tree, partitions, tree->start, PLL_TRUE, PLL_FALSE);
	evaluateSPR(tree, partitions);

	Tree2String(tree->tree_string, tree, partitions, tree->start->back,
			PLL_TRUE, PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE,
			PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);

	model->setLnL(tree->likelihood);
	model->setTree(tree->tree_string);
	model->setFrequencies(partitions->partitionData[0]->frequencies);
	if (model->isGamma())
		model->setAlpha(partitions->partitionData[0]->alpha);
	model->setRates(partitions->partitionData[0]->substRates);

	return 0;
#endif
}

} /* namespace partest */
