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
#include "parser/newick/newick.h"
#include "utils.h"
#include "parser/partition/part.h"

extern double treeOptimizeRapid(pllInstance *tr, partitionList *pr, int mintrav,
		int maxtrav, analdef *adef, bestlist *bt, infoList *iList);
}

namespace partest {

double PLLModelOptimize::evaluate(pllInstance * tr, partitionList *pr,
		analdef * adef, bool estimateModel) {
	int i, impr, bestTrav = 0, rearrangementsMax = 0, rearrangementsMin = 0,
			thoroughIterations = 0, fastIterations = 0;

	double lh = PLL_UNLIKELY, previousLh = PLL_UNLIKELY, difference, epsilon;

	bestlist *bestT, *bt;
	infoList *iList = (infoList*) rax_malloc(sizeof(infoList));

	/* initialize two lists of size 1 and size 20 that will keep track of the best
	 and 20 best tree topologies respectively */

	bestT = (bestlist *) rax_malloc(sizeof(bestlist));
	bestT->ninit = 0;
	initBestTree(bestT, 1, tr->mxtips);

	bt = (bestlist *) rax_malloc(sizeof(bestlist));
	bt->ninit = 0;
	initBestTree(bt, 20, tr->mxtips);

	/* initialize an additional data structure used by the search algo, all of this is pretty
	 RAxML-specific and should probably not be in the library */

	iList->n = 50;
	iList->valid = 0;
	iList->list = (bestInfo *) rax_malloc(sizeof(bestInfo) * (size_t) 50);

	for (int i = 0; i < 50; i++) {
		iList->list[i].node = (nodeptr) NULL;
		iList->list[i].likelihood = PLL_UNLIKELY;
	}

	/* some pretty atbitrary thresholds */

	difference = 10.0;
	epsilon = 0.01;

	/* Thorough = 0 means that we will do fast SPR inbsertions without optimizing the
	 three branches adjacent to the subtree insertion position via Newton-Raphson
	 */

	tr->thoroughInsertion = PLL_FALSE;
	if (estimateModel)
		modOpt(tr, pr, 10.0);
	else
		treeEvaluate(tr, pr, 64);

	/* save the current tree (which is the input tree parsed via -t in the bestT list */
	saveBestTree(bestT, tr,
			pr->perGeneBranchLengths ? pr->numberOfPartitions : 1);

	bestTrav = adef->bestTrav = adef->initial;

	/* optimize model params more thoroughly or just optimize branch lengths */
	if (estimateModel)
		modOpt(tr, pr, 5.0);
	else
		treeEvaluate(tr, pr, 32);   // 32 * 1

	if (tr->doCutoff)
		tr->itCount = 0;
	impr = 1;
	while (impr) {
		recallBestTree(bestT, 1, tr, pr);
		if (tr->searchConvergenceCriterion) {
			int bCounter = 0;
			cleanupHashTable(tr->h, (fastIterations % 2));
			bitVectorInitravSpecial(tr->bitVectors, tr->nodep[1]->back,
					tr->mxtips, tr->vLength, tr->h, fastIterations % 2,
					BIPARTITIONS_RF, (branchInfo *) NULL, &bCounter, 1,
					PLL_FALSE, PLL_FALSE, tr->threadID);
			{
				char *buffer = (char*) rax_calloc((size_t) tr->treeStringLength,
						sizeof(char));
				Tree2String(buffer, tr, pr, tr->start->back, PLL_FALSE,
						PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE,
						PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);
				if (fastIterations % 2 == 0)
					memcpy(tr->tree0, buffer,
							tr->treeStringLength * sizeof(char));
				else
					memcpy(tr->tree1, buffer,
							tr->treeStringLength * sizeof(char));
				rax_free(buffer);
			}
			assert(bCounter == tr->mxtips - 3);
			if (fastIterations > 0) {
				double rrf = convergenceCriterion(tr->h, tr->mxtips);
				if (rrf <= 0.01) /* 1% cutoff */
				{
					cleanupHashTable(tr->h, 0);
					cleanupHashTable(tr->h, 1);
					goto cleanup_fast;
				}
			}
		}

		/* count how many fast iterations with so-called fast SPR moves we have executed */
		fastIterations++;

		/* optimize branch lengths */
		treeEvaluate(tr, pr, 32);  // 32 * 1 = 32

		/* save the tree with those branch lengths again */
		saveBestTree(bestT, tr,
				pr->perGeneBranchLengths ? pr->numberOfPartitions : 1);

		/* update the current best likelihood */
		lh = previousLh = tr->likelihood;

		/* in here we actually do a cycle of SPR moves */
		treeOptimizeRapid(tr, pr, 1, bestTrav, adef, bt, iList);

		/* set impr to 0 since in the immediately following for loop we check if the SPR moves above have generated
		 a better tree */

		impr = 0;

		/* loop over the 20 best trees generated by the fast SPR moves, and check if they improve the likelihood after all of their branch lengths
		 have been optimized */

		for (i = 1; i <= bt->nvalid; i++) {
			recallBestTree(bt, i, tr, pr);
			treeEvaluate(tr, pr, 8); // 0.25 * 32
			difference = (
					(tr->likelihood > previousLh) ?
							tr->likelihood - previousLh :
							previousLh - tr->likelihood);
			if (tr->likelihood > lh && difference > epsilon) {
				impr = 1;
				lh = tr->likelihood;
				saveBestTree(bestT, tr,
						pr->perGeneBranchLengths ? pr->numberOfPartitions : 1);
			}
		}
	}

	if (tr->searchConvergenceCriterion) {
		cleanupHashTable(tr->h, 0);
		cleanupHashTable(tr->h, 1);
	}

	cleanup_fast: tr->thoroughInsertion = PLL_TRUE;
	impr = 1;
	recallBestTree(bestT, 1, tr, pr);
	evaluateGeneric(tr, pr, tr->start, PLL_TRUE, PLL_FALSE);

	/* optimize model params (including branch lengths) or just
	 optimize branch lengths and leave the other model parameters (GTR rates, alhpa)
	 alone */
//	if (estimateModel)
//		modOpt(tr, pr, 1.0);
//	else
//		treeEvaluate(tr, pr, 32); //32 * 1
	while (1) {
		recallBestTree(bestT, 1, tr, pr);

		if (impr) {
			rearrangementsMin = 1;
			rearrangementsMax = adef->stepwidth;
			if (tr->searchConvergenceCriterion) {
				int bCounter = 0;
				if (thoroughIterations > 1)
					cleanupHashTable(tr->h, (thoroughIterations % 2));
				bitVectorInitravSpecial(tr->bitVectors, tr->nodep[1]->back,
						tr->mxtips, tr->vLength, tr->h, thoroughIterations % 2,
						BIPARTITIONS_RF, (branchInfo *) NULL, &bCounter, 1,
						PLL_FALSE, PLL_FALSE, tr->threadID);
				char *buffer = (char*) rax_calloc((size_t) tr->treeStringLength,
						sizeof(char));

				Tree2String(buffer, tr, pr, tr->start->back, PLL_FALSE,
						PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE,
						PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);

				if (thoroughIterations % 2 == 0)
					memcpy(tr->tree0, buffer,
							tr->treeStringLength * sizeof(char));
				else
					memcpy(tr->tree1, buffer,
							tr->treeStringLength * sizeof(char));
				rax_free(buffer);
				assert(bCounter == tr->mxtips - 3);
				if (thoroughIterations > 0) {
					double rrf = convergenceCriterion(tr->h, tr->mxtips);

					if (rrf <= 0.01) {/* 1% cutoff */
						goto cleanup;
					}
				}
			}
			thoroughIterations++;
		} else {
			rearrangementsMax += adef->stepwidth;
			rearrangementsMin += adef->stepwidth;
			if (rearrangementsMax > adef->max_rearrange) {
				goto cleanup;
			}
		}

		/* optimize branch lengths of best tree */
		treeEvaluate(tr, pr, 32); // 32 * 1

		/* do some bokkeeping and printouts again */
		previousLh = lh = tr->likelihood;
		saveBestTree(bestT, tr,
				pr->perGeneBranchLengths ? pr->numberOfPartitions : 1);
		/* do a cycle of thorough SPR moves with the minimum and maximum rearrangement radii */
		treeOptimizeRapid(tr, pr, rearrangementsMin, rearrangementsMax, adef,
				bt, iList);
		impr = 0;
		for (i = 1; i <= bt->nvalid; i++) {
			recallBestTree(bt, i, tr, pr);

			treeEvaluate(tr, pr, 8); // 0.25	* 32
			difference = (
					(tr->likelihood > previousLh) ?
							tr->likelihood - previousLh :
							previousLh - tr->likelihood);
			if (tr->likelihood > lh && difference > 0.01) { // epsilon) {
				impr = 1;
				lh = tr->likelihood;
				saveBestTree(bestT, tr,
						pr->perGeneBranchLengths ? pr->numberOfPartitions : 1);
			}
		}
	}

	cleanup:
	/* do a final full tree traversal, not sure if this is required here */
	{
		evaluateGeneric(tr, pr, tr->start, PLL_TRUE, PLL_FALSE);
	}

	/* free data structures */
	if (tr->searchConvergenceCriterion) {
		freeBitVectors(tr->bitVectors, 2 * tr->mxtips);
		rax_free(tr->bitVectors);
		freeHashTable(tr->h);
		rax_free(tr->h);
	}

	freeBestTree(bestT);
	rax_free(bestT);
	freeBestTree(bt);
	rax_free(bt);
	rax_free(iList->list);
	rax_free(iList);

	return 0;
}

PLLModelOptimize::PLLModelOptimize(ParTestOptions * options) :
		ModelOptimize(options) {
	cout << "CREATING " << endl;
	alignment = static_cast<PLLAlignment *>(options->getAlignment());
	tr = alignment->getTree();

}

PLLModelOptimize::~PLLModelOptimize() {
}

int PLLModelOptimize::optimizePartitioningSchemeAtOnce(
		PartitioningScheme * scheme) {

	/* build the whole alignment */
	string pllPartitionsFile("pllConfig2.tmp");
	ofstream pllOutputStream(pllPartitionsFile.c_str());

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

	struct pllQueue * parts = pllPartitionParse("pllConfig2.tmp");
	partitionList * partitions = pllPartitionsCommit(parts,
			alignment->getPhylip());
	pllInstance * tr = pllCreateInstance(GAMMA, PLL_FALSE, PLL_FALSE, PLL_FALSE,
			12345);
	pllTreeInitTopologyForAlignment(tr, alignment->getPhylip());
	/* Connect the alignment with the tree structure */
	if (!pllLoadAlignment(tr, alignment->getPhylip(), partitions,
	PLL_SHALLOW_COPY)) {
		cerr << "ERROR: Incompatible tree/alignment combination" << endl;
		Utilities::exit_partest(EX_SOFTWARE);
	}

	/* Initialize the model TODO: Put the parameters in a logical order and change the TRUE to flags */
	pllInitModel(tr, PLL_TRUE, alignment->getPhylip(), partitions);

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
		makeGammaCats(current_part->alpha,
				current_part->gammaRates, 4, tr->useMedian);
	}

	if (options->getTreeString() == 0) {
		pllComputeRandomizedStepwiseAdditionParsimonyTree(tr, partitions);
	} else {
		struct pllNewickTree * nt = pllNewickParseString (options->getTreeString());
		pllTreeInitTopologyNewick (tr, nt, PLL_TRUE);
	}

	analdef *adef = (analdef*) rax_calloc(1, sizeof(analdef));
	adef->max_rearrange = 100;
	adef->stepwidth = 5;
	adef->initial = 10;
	adef->bestTrav = 10;
	adef->initialSet = PLL_FALSE;
	adef->mode = BIG_RAPID_MODE;
	adef->likelihoodEpsilon = 0.1;
	adef->permuteTreeoptimize = PLL_FALSE;
	adef->perGeneBranchLengths = PLL_FALSE;
	adef->useCheckpoint = PLL_FALSE;

	evaluateGeneric(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
	evaluate(tr, partitions, adef, false);

	rax_free(adef);

	Tree2String(tr->tree_string, tr, partitions, tr->start->back, PLL_TRUE,
			PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH,
			PLL_FALSE, PLL_FALSE);

	cout << "FINAL TREE: " << tr->tree_string << endl;
	return 0;
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
	adef->max_rearrange = 100;
	adef->stepwidth = 50;
	adef->initial = 1;
	adef->bestTrav = 1;
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

	if (options->getTreeString() == 0) {
		pllComputeRandomizedStepwiseAdditionParsimonyTree(tree, partitions);
	} else {
		struct pllNewickTree * nt = pllNewickParseString (options->getTreeString());
		pllTreeInitTopologyNewick (tree, nt, PLL_FALSE);
		pllNewickParseDestroy(&nt);
	}

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
	evaluate(tree, partitions, adef);

	rax_free(adef);

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
