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

typedef struct {
	pllInstance * tr;
	nodeptr p;
	int nniType;
	double z[NUM_BRANCHES]; // optimize branch lengths
	double z0[NUM_BRANCHES]; // unoptimized branch lengths
	double likelihood;
	double deltaLH;
} nniMove;

void PLLModelOptimize::initializeStructs(pllInstance * tree,
		partitionList * partitions, pllPhylip * phylip) {

	pllTreeInitTopologyForAlignment(tree, phylip);

	/* Connect the alignment with the tree structure */
	if (!pllLoadAlignment(tree, phylip, partitions, PLL_SHALLOW_COPY)) {
		cerr << "ERROR: Incompatible tree/alignment combination" << endl;
		Utilities::exit_partest(EX_SOFTWARE);
	}

	/* Initialize the model TODO: Put the parameters in a logical order and change the TRUE to flags */
	pllInitModel(tree, PLL_TRUE, phylip, partitions);

}

double doOneNNI(pllInstance* tr, partitionList *pr, nodeptr p, int swap,
		int optBran) {
	nodeptr q;
	nodeptr tmp;

	q = p->back;
	//printTopology(tr, TRUE);
	assert(!isTip(q->number, tr->mxtips));
	assert(!isTip(p->number, tr->mxtips));

	int pNum = p->number;
	int qNum = q->number;

	if (swap == 1) {
		tmp = p->next->back;
		hookup(p->next, q->next->back, q->next->z, pr->numberOfPartitions);
		hookup(q->next, tmp, tmp->z, pr->numberOfPartitions);
		//hookupDefault(p->next, q->next->back, tr->numBranches);
		//hookupDefault(q->next, tmp, tr->numBranches);
	} else {
		tmp = p->next->next->back;
		hookup(p->next->next, q->next->back, q->next->z,
				pr->numberOfPartitions);
		hookup(q->next, tmp, tmp->z, pr->numberOfPartitions);
		//hookup(p->next->next, q->next->back, q->next->z, tr->numBranches);
		//hookup(q->next, tmp, tmp->z, tr->numBranches);
	}

	assert(pNum == p->number);
	assert(qNum == q->number);

	if (optBran) {
		newviewGeneric(tr, pr, p, PLL_FALSE);
		newviewGeneric(tr, pr, q, PLL_FALSE);
		update(tr, pr, p);
		//printf("New branch length %f \n", getBranchLength(tr, 0, p) );
		evaluateGeneric(tr, pr, p, PLL_FALSE, PLL_FALSE);
		return tr->likelihood;
	} else {
		//newviewGeneric(tr, p, FALSE);
		//newviewGeneric(tr, q, FALSE);
		return -1.0;
	}
}

nniMove getBestNNIForBran(pllInstance* tr, partitionList *pr, nodeptr p,
		double curLH) {
	nodeptr q = p->back;
	assert( ! isTip(p->number, tr->mxtips));
	assert( ! isTip(q->number, tr->mxtips));
#ifdef DEBUG_MAX
	Tree2String(tr->tree_string, tr, tr->start->back, TRUE, FALSE, 0, 0, 0, SUMMARIZE_LH, 0,0);
	fprintf(stderr, "%s\n", tr->tree_string);
#endif

	/* Backup the current branch length */
	double z0[NUM_BRANCHES];
	int i;
	for (i = 0; i < pr->numberOfPartitions; i++) {
		z0[i] = p->z[i];
	}
#ifdef DEBUG_MAX
	double lhOld = tr->likelihood;
	printf("lhOld: %f \n", lhOld);
#endif
	//update(tr, p);
	//evaluateGeneric(tr, p, FALSE);
	//printf("Current tree LH = %f \n", tr->likelihood);
	//double lh0 = tr->likelihood;
	double lh0 = curLH;
	//printf("zNew: %f \n", getBranchLength(tr, 0, p));

#ifdef DEBUG_MAX
	printf("lh0: %f \n", lh0);
#endif
	nniMove nni0; // nni0 means no NNI move is done
	nni0.p = p;
	nni0.nniType = 0;
	nni0.deltaLH = 0;
	for (i = 0; i < pr->numberOfPartitions; i++) {
		nni0.z[i] = p->z[i];
	}

	/* TODO Save the likelihood vector at node p and q */
	//saveLHVector(p, q, p_lhsave, q_lhsave);
	/* Save the scaling factor */
	// Now try to do an NNI move of type 1
	double lh1 = doOneNNI(tr, pr, p, 1, PLL_TRUE);
	nniMove nni1;
	nni1.p = p;
	nni1.nniType = 1;
	// Store the optimized und unoptimized central branch length
	for (i = 0; i < pr->numberOfPartitions; i++) {
		nni1.z[i] = p->z[i];
		nni1.z0[i] = z0[i];
	}
	nni1.likelihood = lh1;
	nni1.deltaLH = lh1 - lh0;
#ifdef DEBUG_MAX
	printf("Delta likelihood of the 1.NNI move: %f\n", nni1.deltaLH);
	//printTopology(tr, TRUE);
#endif

	/* Restore previous NNI move */
	doOneNNI(tr, pr, p, 1, PLL_FALSE);
	/* Restore the old branch length */
	for (i = 0; i < pr->numberOfPartitions; i++) {
		p->z[i] = z0[i];
		p->back->z[i] = z0[i];
	}

#ifdef DEBUG_MAX
	printf("Restore topology\n");
	Tree2String(tr->tree_string, tr, tr->start->back, TRUE, FALSE, 0, 0, 0, SUMMARIZE_LH, 0,0);
	fprintf(stderr, "%s\n", tr->tree_string);
	evaluateGeneric(tr, tr->start, TRUE);
	printf("Likelihood after restoring from NNI 1: %f\n", tr->likelihood);
#endif
	/* Try to do an NNI move of type 2 */
	double lh2 = doOneNNI(tr, pr, p, 2, PLL_TRUE);
	// Create the nniMove struct to store this move
	nniMove nni2;
	nni2.p = p;
	nni2.nniType = 2;
	// Store the optimized and unoptimized central branch length
	for (i = 0; i < pr->numberOfPartitions; i++) {
		nni2.z[i] = p->z[i];
		nni2.z0[i] = z0[i];
	}
	nni2.likelihood = lh2;
	nni2.deltaLH = lh2 - lh0;
#ifdef DEBUG_MAX
	printf("Delta likelihood of the 2.NNI move: %f\n", nni2.deltaLH);
	//printTopology(tr, TRUE);
#endif

	/* Restore previous NNI move */
	doOneNNI(tr, pr, p, 2, PLL_FALSE);
	newviewGeneric(tr, pr, p, PLL_FALSE);
	newviewGeneric(tr, pr, p->back, PLL_FALSE);
	/* Restore the old branch length */
	for (i = 0; i < pr->numberOfPartitions; i++) {
		p->z[i] = z0[i];
		p->back->z[i] = z0[i];
	}
	if (nni1.deltaLH > 0 && nni1.deltaLH >= nni2.deltaLH) {
		return nni1;
	} else if (nni1.deltaLH > 0 && nni1.deltaLH < nni2.deltaLH) {
		return nni2;
	} else if (nni1.deltaLH < 0 && nni2.deltaLH > 0) {
		return nni2;
	} else {
		return nni0;
	}
	/******************** NNI part *******************************/

	/* Restore the likelihood vector */
	//restoreLHVector(p,q, p_lhsave, q_lhsave);
}

void evalNNIForSubtree(pllInstance* tr, partitionList *pr, nodeptr p,
		nniMove* nniList, int* cnt, int* cnt_nni, double curLH) {
	if (!isTip(p->number, tr->mxtips)) {
		//newviewGeneric(tr, p, FALSE);
		//newviewGeneric(tr, p->back, FALSE);
		nniList[*cnt] = getBestNNIForBran(tr, pr, p, curLH);
		if (nniList[*cnt].deltaLH != 0.0) {
			*cnt_nni = *cnt_nni + 1;
		}
		*cnt = *cnt + 1;
		nodeptr q = p->next;
		while (q != p) {
			evalNNIForSubtree(tr, pr, q->back, nniList, cnt, cnt_nni, curLH);
			q = q->next;
		}
	}
}

int cmp_nni(const void* nni1, const void* nni2) {
	nniMove* myNNI1 = (nniMove*) nni1;
	nniMove* myNNI2 = (nniMove*) nni2;
	return (int) (1000000.f * myNNI1->deltaLH - 1000000.f * myNNI2->deltaLH);
}

double PLLModelOptimize::evaluateSPR(pllInstance * tr, partitionList *pr, analdef * adef,
			bool estimateModel) {

	int i, impr, bestTrav = 0, rearrangementsMax = 0, rearrangementsMin = 0,
			thoroughIterations = 0, fastIterations = 0;
	time_t t0, t1;
	t0 = time(NULL);
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

	t1 = time(NULL);
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

	t1 = time(NULL);
	if (estimateModel)
		modOpt(tr, pr, 5.0);
	else
		treeEvaluate(tr, pr, 32);   // 32 * 1

	if (tr->doCutoff)
		tr->itCount = 0;
	impr = 1;

	int iters = 0;
	int imprIter = 0;
	t1 = time(NULL);

	while (impr) {
		iters++;
		recallBestTree(bestT, 1, tr, pr);
		if (tr->searchConvergenceCriterion) {
			imprIter++;
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
	iters = 0;
	imprIter = 0;
	t1 = time(NULL);
	while (1) {
		iters++;
		recallBestTree(bestT, 1, tr, pr);

		if (impr) {
			imprIter++;
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
				cout << "EVALUATE: DIFF " << tr->likelihood - previousLh << endl;
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

	return tr->likelihood;
}

double PLLModelOptimize::evaluateNNI(pllInstance * tr, partitionList *pr, bool estimateModel) {

	double curScore = tr->likelihood;

	/* Initialize the NNI list */
	nniMove* nniList = (nniMove*) malloc((tr->mxtips - 3) * sizeof(nniMove));
	int i;
	/* fill up the NNI list */
	nodeptr p = tr->start->back;
	nodeptr q = p->next;
	int cnt = 0; // number of visited internal branches during NNI evaluation
	int cnt_nni = 0; // number of positive NNI found
	while (q != p) {
		evalNNIForSubtree(tr, pr, q->back, nniList, &cnt, &cnt_nni, curScore);
		q = q->next;
	}
	if (cnt_nni == 0)
		return 0.0;

	nniMove* impNNIList = (nniMove*) malloc(cnt_nni * sizeof(nniMove));
	int j = 0;
	for (i = 0; i < tr->mxtips - 3; i++) {
		if (nniList[i].deltaLH > 0.0) {
			impNNIList[j] = nniList[i];
			j++;
		}
	}
	// sort impNNIList
	qsort(impNNIList, cnt_nni, sizeof(nniMove), cmp_nni);

	// creating a list of non-conflicting positive NNI
	nniMove* nonConfNNIList = (nniMove*) calloc(cnt_nni, sizeof(nniMove));

	// the best NNI will always be taken
	nonConfNNIList[0] = impNNIList[cnt_nni - 1];

	// Filter out conflicting NNI
	int numNonConflictNNI = 1; // size of the non-conflicting NNI list;
	int k;
	for (k = cnt_nni - 2; k >= 0; k--) {
		int conflict = PLL_FALSE;
		int j;
		for (j = 0; j < numNonConflictNNI; j++) {
			if (impNNIList[k].p->number == nonConfNNIList[j].p->number
					|| impNNIList[k].p->number
							== nonConfNNIList[j].p->back->number) {
				conflict = PLL_TRUE;
				break;
			}
		}
		if (conflict) {
			continue;
		} else {
			nonConfNNIList[numNonConflictNNI] = impNNIList[k];
			numNonConflictNNI++;
		}
	}

	// Applying non-conflicting NNI moves
	double delta = 1.0; // portion of NNI moves to apply
	int notImproved;
	do {
		notImproved = PLL_FALSE;
		//printf("numNonConflictNNI = %d \n", numNonConflictNNI);
		int numNNI2Apply = ceil(numNonConflictNNI * delta);
		//printf("numNNI2Apply = %d \n", numNNI2Apply);
		for (i = 0; i < numNNI2Apply; i++) {
			// Just do the topological change
			doOneNNI(tr, pr, nonConfNNIList[i].p, nonConfNNIList[i].nniType,
					PLL_FALSE);
			newviewGeneric(tr, pr, nonConfNNIList[i].p, PLL_FALSE);
			newviewGeneric(tr, pr, nonConfNNIList[i].p->back, PLL_FALSE);
			// Apply the store branch length
			int j;
			for (j = 0; j < pr->numberOfPartitions; j++) {
				nonConfNNIList[i].p->z[j] = nonConfNNIList[i].z[j];
				nonConfNNIList[i].p->back->z[j] = nonConfNNIList[i].z[j];
			}
		}
		// Re-optimize all branches
		smoothTree(tr, pr, 2);
		evaluateGeneric(tr, pr, tr->start, PLL_FALSE, PLL_FALSE);
		if (estimateModel) {
			modOpt(tr, pr, 0.1);
		}
		evaluateGeneric(tr, pr, tr->start, PLL_FALSE, PLL_FALSE);
		if (tr->likelihood < curScore) {
			printf("Tree likelihood gets worse after applying NNI\n");
			printf("curScore = %30.20f\n", curScore);
			printf("newScore = %30.20f\n", tr->likelihood);
			printf("Rolling back the tree\n");
			for (i = 0; i < numNNI2Apply; i++) {
				doOneNNI(tr, pr, nonConfNNIList[i].p, nonConfNNIList[i].nniType,
						PLL_FALSE);
				// Restore the branch length
				int j;
				for (j = 0; j < pr->numberOfPartitions; j++) {
					nonConfNNIList[i].p->z[j] = nonConfNNIList[i].z0[j];
					nonConfNNIList[i].p->back->z[j] = nonConfNNIList[i].z0[j];
				}
			}
			evaluateGeneric(tr, pr, tr->start, PLL_FALSE, PLL_FALSE);
			printf("Tree likelihood after rolling back = %f \n",
					tr->likelihood);
			notImproved = PLL_TRUE & (numNNI2Apply > 1);
			delta = delta * 0.5;
		}
	} while (notImproved);
	free(nniList);
	free(impNNIList);
	free(nonConfNNIList);

	return 0;
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

	struct pllQueue * parts = pllPartitionParse(pllPartitionsFile);

	if (remove(pllPartitionsFile) != 0)
		cerr << "Error deleting temporary file" << endl;

	partitionList * partitions = pllPartitionsCommit(parts,
			alignment->getPhylip());
	pllQueuePartitionsDestroy(&parts);

	pllInstance * tr = pllCreateInstance(GAMMA, PLL_FALSE, PLL_FALSE, PLL_FALSE,
			12345);
	pllTreeInitTopologyForAlignment(tr, alignment->getPhylip());
	/* Connect the alignment with the tree structure */
	if (!pllLoadAlignment(tr, alignment->getPhylip(), partitions,
	PLL_DEEP_COPY)) {
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
		struct pllNewickTree * nt = pllNewickParseString(
				options->getTreeString());
		pllTreeInitTopologyNewick(tr, nt, PLL_TRUE);
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
	evaluateNNI(tr, partitions, false);

	rax_free(adef);

	Tree2String(tr->tree_string, tr, partitions, tr->start->back, PLL_TRUE,
			PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH,
			PLL_FALSE, PLL_FALSE);

	cout << "FINAL TREE: " << tr->tree_string << endl;

	pllPartitionsDestroy(&partitions, partitions->numberOfPartitions,
			tr->mxtips);
	pllTreeDestroy(tr);

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
			pllPhylip * phylip = alignment->getPhylip();

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

				struct pllNewickTree * nt = pllNewickParseString(
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
	evaluateNNI(tree, partitions);

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
