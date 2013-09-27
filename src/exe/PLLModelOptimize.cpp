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

char convert3(unsigned char c) {
	switch (c) {
	case 1:
		return 'A';
	case 2:
		return 'C';
	case 4:
		return 'G';
	case 8:
		return 'T';
	default:
		return 'X';
	}
}

void PLLModelOptimize::initializeStructs(pllInstance * tree,
		partitionList * partitions, pllAlignmentData * phylip) {

#ifdef DEBUG
	cout << "[TRACE] PLLModelOptimize - Initializing topology" << endl;
#endif

	pllTreeInitTopologyForAlignment(tree, phylip);

#ifdef DEBUG
	cout << "[TRACE] PLLModelOptimize - Loading alignment" << endl;
#endif
	/* Connect the alignment with the tree structure */
	if (!pllLoadAlignment(tree, phylip, partitions, PLL_SHALLOW_COPY)) {
		cerr << "ERROR: Incompatible tree/alignment combination" << endl;
		Utilities::exit_partest(EX_SOFTWARE);
	}

#ifdef DEBUG
	cout << "[TRACE] PLLModelOptimize - Initializing model" << endl;
#endif
	pllInitModel(tree, partitions, phylip);
#ifdef DEBUG
	cout << "[TRACE] PLLModelOptimize - Initialized model" << endl;
#endif

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

double PLLModelOptimize::optimizeParameters(pllInstance * tr,
		partitionList *partitions, bool estimateModel,
		bool estimateBranchLengths, bool estimateTopology) {

	double lk;
	double epsilon = 0.1;

	evaluateGeneric(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);

	do {
		lk = tr->likelihood;

		// TODO: Optimize topology?

		if (estimateModel)
			pllOptimizeModelParameters(tr, partitions, 0.1);

		if (estimateBranchLengths)
			pllTreeEvaluate(tr, partitions, 200);

		evaluateGeneric(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
	} while (fabs(lk - tr->likelihood) > epsilon);

	return tr->likelihood;
}
double PLLModelOptimize::evaluateSPR(pllInstance * tr,
		partitionList *partitions, bool estimateModel, bool estimateTopology) {

	pllRearrangeList * bestList = pllCreateRearrangeList(1);
	int i;
	double lk = PLL_UNLIKELY, epsilon = 0.01;

	time_t t0 = time(NULL);

	cout << endl << "[SPR] " << Utilities::timeToString(time(NULL) - t0)
			<< " START" << endl;

	evaluateGeneric(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);

	cout << "[SPR] " << Utilities::timeToString(time(NULL) - t0)
			<< " Initial Likelihood: " << tr->likelihood << endl;

	Tree2String(tr->tree_string, tr, partitions, tr->start->back, PLL_TRUE,
			PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH,
			PLL_FALSE, PLL_FALSE);

	cout << "[SPR] FIRST TREE -> " << tr->tree_string;

	tr->thoroughInsertion = PLL_TRUE;
	double modelLkThreshold = 0.1;
	int SPRdistance = 2 * tr->mxtips - 2;
	int startNode = 1;

	do {
		lk = tr->likelihood;

		cout << "[SPR] " << Utilities::timeToString(time(NULL) - t0)
				<< " Computing SPR best moves..." << endl;

		if (estimateTopology) {
			pllRearrangeSearch(tr, partitions, PLL_REARRANGE_SPR,
					tr->nodep[tr->mxtips + 1], 1, SPRdistance, bestList);
//
			pllRearrangeCommit(tr, partitions, bestList->rearr, PLL_TRUE);
		}

		pllTreeEvaluate(tr, partitions, 200);

		//pllCommitSPR(tr, partitions, &(bestList->sprInfo[0]), PLL_TRUE);

		tr->thoroughInsertion = PLL_FALSE;
		//pllDestroyListSPR(&bestList);

		evaluateGeneric(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);

		cout << "[SPR] " << Utilities::timeToString(time(NULL) - t0)
				<< " New Likelihood (TOPO): " << tr->likelihood << endl;

		if (lk > tr->likelihood) {
			if (estimateTopology)
				pllRearrangeRollback(tr, partitions);
			//pllRollbackSPR(tr, partitions);
			evaluateGeneric(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
			cout << "[SPR] " << Utilities::timeToString(time(NULL) - t0)
					<< " Rollback (TOPO): " << tr->likelihood << endl;
			break;
		} else {
			Tree2String(tr->tree_string, tr, partitions, tr->start->back,
					PLL_TRUE, PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE,
					PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);
		}

		if (estimateModel) {
			cout << "[SPR] " << Utilities::timeToString(time(NULL) - t0)
					<< " Optimizing model parameters..." << endl;
			pllOptimizeModelParameters(tr, partitions, modelLkThreshold);
			evaluateGeneric(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
			cout << "[SPR] " << Utilities::timeToString(time(NULL) - t0)
					<< " New Likelihood (MODEL): " << tr->likelihood << endl;
		}

		cout << "[SPR] " << Utilities::timeToString(time(NULL) - t0)
				<< " Loop diff: " << fabs(lk - tr->likelihood) << endl;

	} while (fabs(lk - tr->likelihood) > epsilon);

	pllDestroyRearrangeList(&bestList);

	cout << "[SPR] " << Utilities::timeToString(time(NULL) - t0) << " Done loop"
			<< endl;

	Tree2String(tr->tree_string, tr, partitions, tr->start->back, PLL_TRUE,
			PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH,
			PLL_FALSE, PLL_FALSE);

	cout << "[SPR] LAST TREE -> " << tr->tree_string;
	cout << "[SPR] " << Utilities::timeToString(time(NULL) - t0) << " END"
			<< endl << endl;

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

#ifdef DEBUG
	cout << "[TRACE] PLLModelOptimize - Constructing partitions structure" << endl;
#endif

	struct pllQueue * parts;
	struct pllPartitionRegion * pregion;
	struct pllPartitionInfo * pinfo;

	pllQueueInit(&parts);

	int numSeqs = alignment->getNumSeqs();
	int numSites = alignment->getNumSites();
	pllAlignmentData * phylip = pllInitAlignmentData(numSeqs, numSites);
	phylip->sequenceCount = numSeqs;
	phylip->sequenceLength = numSites;
	for (int i = 0; i < numSeqs; i++) {
		phylip->sequenceLabels[i + 1] = strdup(
				alignment->getPhylip()->sequenceLabels[i + 1]);
	}
	phylip->siteWeights = (int *) malloc(numSites * sizeof(int));
	int firstSite = 1;
	int nextSite = 0;
	for (int i = 0; i < scheme->getNumberOfElements(); i++) {
		PartitionElement * element = scheme->getElement(i);

		pinfo = (pllPartitionInfo *) malloc(sizeof(struct pllPartitionInfo));
		pllQueueInit(&(pinfo->regionList));
		pllQueueAppend(parts, (void *) pinfo);

		pinfo->partitionName = (char *) malloc(
				(element->getName().size() + 1) * sizeof(char));
		strcpy(pinfo->partitionName, element->getName().c_str());
		pinfo->partitionModel = (char *) malloc(5 * sizeof(char));
		if (element->getBestModel()->getModel()->isPF()) {
			strcpy(pinfo->partitionModel, "DNAX");
			pinfo->optimizeBaseFrequencies = PLL_TRUE;
		} else {
			strcpy(pinfo->partitionModel, "DNA");
			pinfo->optimizeBaseFrequencies = PLL_FALSE;
		}
		pinfo->protModels = -1;
		pinfo->protFreqs = -1;
		pinfo->dataType = DNA_DATA;

		for (int j = 0; j < element->getNumberOfSections(); j++) {
			for (int site = element->getStart(j) - 1; site < element->getEnd(j);
					site++) {
				phylip->siteWeights[nextSite] = 1;
				for (int nk = 0; nk < numSeqs; nk++) {
					phylip->sequenceData[nk + 1][nextSite] =
							alignment->getPhylip()->sequenceData[nk + 1][site];
				}
				nextSite++;
			}
		}

		pregion = (struct pllPartitionRegion *) malloc(
				sizeof(struct pllPartitionRegion));
		pregion->start = firstSite;
		pregion->end = nextSite;
		pregion->stride = 1;
		pllQueueAppend(pinfo->regionList, (void *) pregion);

		firstSite = nextSite + 1;
	}

#ifdef DEBUG
	cout << "[TRACE] PLLModelOptimize - Committing partitions" << endl;
#endif
	partitionList * partitions = pllPartitionsCommit(parts, phylip);

	pllQueuePartitionsDestroy(&parts);

	for (int i = 0; i < partitions->numberOfPartitions; i++) {
		if (partitions->partitionData[i]->lower
				== partitions->partitionData[i]->upper) {
			exit(-1);
		}
	}

	pllPhylipRemoveDuplicate(phylip, partitions);

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
	pllTreeInitTopologyForAlignment(tr, phylip);
	/* Connect the alignment with the tree structure */

#ifdef DEBUG
	cout << "[TRACE] PLLModelOptimize - Connecting tree with alignment" << endl;
#endif
	if (!pllLoadAlignment(tr, phylip, partitions,
	PLL_SHALLOW_COPY)) {
		cerr << "ERROR: Incompatible tree/alignment combination" << endl;
		Utilities::exit_partest(EX_SOFTWARE);
	}

#ifdef DEBUG
	cout << "[TRACE] PLLModelOptimize - Initializing model" << endl;
#endif
	pllInitModel(tr, partitions, phylip);
#ifdef DEBUG
	cout << "[TRACE] PLLModelOptimize - Initialized model" << endl;
#endif

	/* set best-fit models */
	for (int i = 0; i < scheme->getNumberOfElements(); i++) {
		Model * model = scheme->getElement(i)->getBestModel()->getModel();
		setModelParameters(model, tr, partitions, i);
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
		pllTreeInitTopologyNewick(tr, nt, PLL_FALSE);
		Tree2String(tr->tree_string, tr, partitions, tr->start->back, PLL_TRUE,
				PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH,
				PLL_FALSE, PLL_FALSE);
	}

	evaluateGeneric(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
	evaluateSPR(tr, partitions, options->getOptimizeMode() == OPT_GTR, false);

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

			cout << "** ELEMENT " << element->getName() << endl;

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

	if (options->getStartingTopology() == StartTopoFIXED) {
		pllNewickTree * nt = pllNewickParseString(options->getTreeString());
		pllTreeInitTopologyNewick(tree, nt, PLL_FALSE);
		pllNewickParseDestroy(&nt);
	} else if (options->getStartingTopology() == StartTopoMP) {
		pllComputeRandomizedStepwiseAdditionParsimonyTree(tree, partitions);
	}

	Tree2String(tree->tree_string, tree, partitions, tree->start->back,
			PLL_TRUE, PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE,
			PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);

	setModelParameters(model, tree, partitions, 0, false);

	evaluateGeneric(tree, partitions, tree->start, PLL_TRUE, PLL_FALSE);

	optimizeParameters(tree, partitions, true, true, false);

	Tree2String(tr->tree_string, tr, partitions, tr->start->back, PLL_TRUE,
			PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH,
			PLL_FALSE, PLL_FALSE);

	model->setLnL(tree->likelihood);
	model->setTree(tree->tree_string);
	model->setFrequencies(partitions->partitionData[0]->frequencies);
	if (model->isGamma())
		model->setAlpha(partitions->partitionData[0]->alpha);
	model->setRates(partitions->partitionData[0]->substRates);

	return 0;
#endif
}

void PLLModelOptimize::setModelParameters(Model * model, pllInstance * tr,
		partitionList * partitions, int index, bool setAlphaFreqs) {

	const char * m = model->getMatrixName().c_str();
	char * symmetryPar = (char *) malloc(12 * sizeof(char));
	symmetryPar[0] = m[0];
	symmetryPar[11] = '\0';
	for (int j = 1; j < 6; j++) {
		symmetryPar[(j - 1) * 2 + 1] = ',';
		symmetryPar[j * 2] = m[j];
	}

	pllSetSubstitutionRateMatrixSymmetries(symmetryPar, partitions, index);

	pInfo * current_part = partitions->partitionData[index];
	current_part->optimizeBaseFrequencies = model->isPF();
	current_part->alpha = model->getAlpha();

	if (setAlphaFreqs) {
		memcpy(current_part->frequencies, model->getFrequencies(),
				4 * sizeof(double));
		memcpy(current_part->substRates, model->getRates(), 6 * sizeof(double));
	} else {
		if (!model->isPF()) {
			partitions->partitionData[index]->optimizeBaseFrequencies =
					PLL_FALSE;
			for (int i = 0; i < 4; i++) {
				partitions->partitionData[index]->frequencies[i] = 0.25;
			}
		}
		for (int i = 0; i < 6; i++) {
			partitions->partitionData[index]->substRates[i] = 1;
		}

		partitions->partitionData[index]->alpha = 100;
	}

	free(symmetryPar);
	initReversibleGTR(tr, partitions, index);
	if (setAlphaFreqs) {
		makeGammaCats(current_part->alpha, current_part->gammaRates, 4,
				tr->useMedian);
	}
}

} /* namespace partest */
