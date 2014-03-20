/*
 * PLLModelOptimize.cpp
 *
 *  Created on: Jan 8, 2013
 *      Author: diego
 */

#include "PLLModelOptimize.h"
#include <iostream>
#include <sstream>
#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

namespace partest {

void PLLModelOptimize::initializeStructs(pllInstance * tree,
		partitionList * partitions, pllAlignmentData * phylip) {

#ifdef DEBUG
	cout << "[TRACE] PLLModelOptimize - Loading alignment" << endl;
#endif
	/* Connect the alignment with the tree structure */
	if (!pllLoadAlignment(tree, phylip, partitions, PLL_SHALLOW_COPY)) {
		cerr << "ERROR: Incompatible tree/alignment combination" << endl;
		Utilities::exit_partest(EX_SOFTWARE);
	}
	tree->start = tree->nodep[1];

#ifdef DEBUG
	cout << "[TRACE] PLLModelOptimize - Building initial MP topology" << endl;
#endif
	pllComputeRandomizedStepwiseAdditionParsimonyTree(tree, partitions);

#ifdef DEBUG
	cout << "[TRACE] PLLModelOptimize - Initializing model (str)" << endl;
#endif
	if (partitions->partitionData[0]->dataType == PLL_AA_DATA) {
		for (int i = 0; i < partitions->numberOfPartitions; i++) {
			partitions->partitionData[i]->dataType = PLL_AA_DATA;
			partitions->partitionData[i]->states = 20;
			partitions->partitionData[i]->protModels = PLL_DAYHOFF;
		}
	}
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

	tr->thoroughInsertion = PLL_FALSE;

	pllEvaluateLikelihood(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);

	do {
		lk = tr->likelihood;
		if (estimateModel)
			pllOptimizeModelParameters(tr, partitions, 1);
		if (estimateBranchLengths)
			pllOptimizeBranchLengths(tr, partitions, 10);
		pllEvaluateLikelihood(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
	} while (fabs(lk - tr->likelihood) > epsilon);

	return tr->likelihood;
}

#ifdef _WIN32
string ExePath() {
	char buffer[255];
	GetModuleFileName( NULL, buffer, 255 );
	string::size_type pos = string( buffer ).find_last_of( "\\/" );
	return string( buffer ).substr( 0, pos);
}
#endif

double PLLModelOptimize::evaluateSPR(pllInstance * tr,
		partitionList *partitions, bool estimateModel, bool estimateTopology) {

	pllRaxmlSearchAlgorithm(tr, partitions, estimateModel);

	return tr->likelihood;
}

double prevSPR(pllInstance * tr, partitionList *partitions, bool estimateModel,
		bool estimateTopology) {

	pllRearrangeList * bestList = pllCreateRearrangeList(1);
	double lk = PLL_UNLIKELY, epsilon = 0.01;

#ifdef DEBUG
	time_t t0 = time(NULL);
	cout << endl << "[SPR] " << Utilities::timeToString(time(NULL) - t0)
			<< " START" << endl;
#endif

	pllEvaluateLikelihood(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);

#ifdef DEBUG
	cout << "[SPR] " << Utilities::timeToString(time(NULL) - t0)
			<< " Initial Likelihood: " << tr->likelihood << endl;

	pllTreeToNewick(tr->tree_string, tr, partitions, tr->start->back, PLL_TRUE,
	PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH,
	PLL_FALSE, PLL_FALSE);

	cout << "[SPR] FIRST TREE -> " << tr->tree_string;
#endif

	tr->thoroughInsertion = PLL_TRUE;
	double modelLkThreshold = 0.1;
	int SPRdistance = 2 * tr->mxtips - 2;

	do {
#ifndef DEBUG
		cout << ".";
		cout.flush();
#endif
		lk = tr->likelihood;

#ifdef DEBUG
		cout << "[SPR] " << Utilities::timeToString(time(NULL) - t0)
				<< " Computing SPR best moves..." << endl;
#endif

		if (estimateTopology) {
			pllRearrangeSearch(tr, partitions, PLL_REARRANGE_SPR,
					tr->nodep[tr->mxtips + 1], 1, SPRdistance, bestList);
			pllRearrangeCommit(tr, partitions, bestList->rearr, PLL_TRUE);
		}

		pllOptimizeBranchLengths(tr, partitions, 200);

		tr->thoroughInsertion = PLL_FALSE;

		pllEvaluateLikelihood(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);

#ifdef DEBUG
		cout << "[SPR] " << Utilities::timeToString(time(NULL) - t0)
				<< " New Likelihood (TOPO): " << tr->likelihood << endl;
#endif

		if (lk > tr->likelihood) {
			if (estimateTopology) {
				pllRearrangeRollback(tr, partitions);
			}
			pllEvaluateLikelihood(tr, partitions, tr->start, PLL_TRUE,
					PLL_FALSE);

#ifdef DEBUG
			cout << "[SPR] " << Utilities::timeToString(time(NULL) - t0)
					<< " Rollback (TOPO): " << tr->likelihood << endl;
#endif
			break;
		} else {
			pllTreeToNewick(tr->tree_string, tr, partitions, tr->start->back,
			PLL_TRUE, PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE,
			PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);
		}

		if (estimateModel) {

#ifdef DEBUG
			cout << "[SPR] " << Utilities::timeToString(time(NULL) - t0)
					<< " Optimizing model parameters..." << endl;
#endif

			pllOptimizeModelParameters(tr, partitions, modelLkThreshold);
			pllEvaluateLikelihood(tr, partitions, tr->start, PLL_TRUE,
					PLL_FALSE);

#ifdef DEBUG
			cout << "[SPR] " << Utilities::timeToString(time(NULL) - t0)
					<< " New Likelihood (MODEL): " << tr->likelihood << endl;
#endif
		}

#ifdef DEBUG
		cout << "[SPR] " << Utilities::timeToString(time(NULL) - t0)
				<< " Loop diff: " << fabs(lk - tr->likelihood) << endl;
#endif

	} while (fabs(lk - tr->likelihood) > epsilon);

	pllDestroyRearrangeList(&bestList);

	pllTreeToNewick(tr->tree_string, tr, partitions, tr->start->back, PLL_TRUE,
	PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH,
	PLL_FALSE, PLL_FALSE);

#ifdef DEBUG
	cout << "[SPR] " << Utilities::timeToString(time(NULL) - t0) << " Done loop"
			<< endl;
	cout << "[SPR] LAST TREE -> " << tr->tree_string;
	cout << "[SPR] " << Utilities::timeToString(time(NULL) - t0) << " END"
			<< endl << endl;
#endif

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

	char * tree = getMlTree(scheme, options->getInputFile());
	scheme->setTree(tree);

	return 0;
}

char * PLLModelOptimize::getMlTree(PartitioningScheme * scheme,
		string inputFile) {
#ifdef DEBUG
	cout << "[TRACE] PLLModelOptimize - Constructing partitions structure"
			<< endl;
#endif

	pllQueue * parts;
	pllPartitionRegion * pregion;
	pllPartitionInfo * pinfo;

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

#ifdef DEBUG
		cout << "[TRACE] PLLModelOptimize - " << element->getName() << endl;
#endif

		pinfo = (pllPartitionInfo *) malloc(sizeof(pllPartitionInfo));
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
		if (options->getDataType()==DT_NUCLEIC) {
			pinfo->protModels = -1;
			pinfo->protFreqs = -1;
		} else {
			ProteicModel * pModel = ((ProteicModel *) dynamic_cast<ProteicModel *>(element->getBestModel()->getModel()));
			pinfo->protModels = pModel->getMatrix();
			pinfo->protFreqs = pModel->isPF();
		}
		pinfo->dataType = options->getDataType()==DT_NUCLEIC?PLL_DNA_DATA:PLL_AA_DATA;

#ifdef DEBUG
		cout << "[TRACE] PLLModelOptimize - Loop in sections! "
				<< element->getNumberOfSections() << endl;
#endif

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
		pregion = (pllPartitionRegion *) malloc(sizeof(pllPartitionRegion));
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

	//pllPhylipRemoveDuplicate(phylip, partitions);

#ifdef DEBUG
	cout << "[TRACE] PLLModelOptimize - Creating tree instance" << endl;
#endif
	pllInstanceAttr * attr = (pllInstanceAttr *) malloc(
			sizeof(pllInstanceAttr));
	attr->rateHetModel = PLL_GAMMA;
	attr->fastScaling = PLL_FALSE;
	attr->saveMemory = PLL_FALSE;
	attr->useRecom = PLL_FALSE;
	attr->randomNumberSeed = 12345;
	attr->numberOfThreads = 1;
	pllInstance * tr = pllCreateInstance(attr);
	free(attr);

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

	if (options->getTreeString() == 0) {
		pllComputeRandomizedStepwiseAdditionParsimonyTree(tr, partitions);
	} else {
		pllNewickTree * nt = pllNewickParseString(options->getTreeString());
		pllTreeInitTopologyNewick(tr, nt, PLL_FALSE);
		pllTreeToNewick(tr->tree_string, tr, partitions, tr->start->back,
				PLL_TRUE,
				PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH,
				PLL_FALSE, PLL_FALSE);
	}


	pllInitModel(tr, partitions, phylip);

	/* set best-fit models */
	for (int i = 0; i < scheme->getNumberOfElements(); i++) {
		Model * model = scheme->getElement(i)->getBestModel()->getModel();
		setModelParameters(model, tr, partitions, i);
		pllInitReversibleGTR(tr, partitions, i);
		pllMakeGammaCats(partitions->partitionData[i]->alpha, partitions->partitionData[i]->gammaRates, 4, tr->useMedian);
	}

	pllEvaluateLikelihood(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
	//evaluateSPR(tr, partitions, options->getOptimizeMode() == OPT_GTR, true);
	evaluateSPR(tr, partitions, true, true);

	pllTreeToNewick(tr->tree_string, tr, partitions, tr->start->back, PLL_TRUE,
	PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH,
	PLL_FALSE, PLL_FALSE);

	cout << "final likelihood = " << tr->likelihood << endl;

	char * newick = (char *)malloc(strlen(tr->tree_string));
	strcpy(newick, tr->tree_string);
	newick[strlen(newick)-1] = '\0';

	pllAlignmentDataDestroy(phylip);

	pllPartitionsDestroy(tr, &partitions);
	partitions = 0;

	pllDestroyInstance(tr);
	tr = 0;

	return newick;
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

	pllTreeToNewick(tree->tree_string, tree, partitions, tree->start->back,
	PLL_TRUE, PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE,
	PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);

	//evaluateSPR(tree, partitions, PLL_TRUE, PLL_TRUE);
	setModelParameters(model, tree, partitions, 0, false);

	pllInitReversibleGTR(tr, partitions, 0);

	optimizeParameters(tree, partitions, true, true, false);

	pllTreeToNewick(tree->tree_string, tree, partitions, tree->start->back,
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

void PLLModelOptimize::setModelParameters(Model * model, pllInstance * tr,
		partitionList * partitions, int index, bool setAlphaFreqs) {

	DataType dataType = options->getDataType();
	pInfo * current_part = partitions->partitionData[index];

	if (dataType == DT_NUCLEIC) {
		current_part->optimizeBaseFrequencies = model->isPF();
		current_part->alpha = model->getAlpha();
		current_part->nonGTR = PLL_FALSE;
		current_part->dataType = PLL_DNA_DATA;
		const char * m = model->getMatrixName().c_str();
		char * symmetryPar = (char *) malloc(12 * sizeof(char));
		symmetryPar[0] = m[0];
		symmetryPar[11] = '\0';
		for (int j = 1; j < 6; j++) {
			symmetryPar[(j - 1) * 2 + 1] = ',';
			symmetryPar[j * 2] = m[j];
		}

		pllSetSubstitutionRateMatrixSymmetries(symmetryPar, partitions, index);

		if (setAlphaFreqs) {
			memcpy(current_part->frequencies, model->getFrequencies(),
					4 * sizeof(double));
			memcpy(current_part->substRates, model->getRates(),
					6 * sizeof(double));
			current_part->alpha = model->getAlpha();
		} else {
			current_part->optimizeBaseFrequencies =
					model->isPF();
			if (!model->isPF()) {
				for (int i = 0; i < 4; i++) {
					current_part->frequencies[i] = 0.25;
				}
			}
			for (int i = 0; i < 6; i++) {
				current_part->substRates[i] = 1;
			}

			current_part->alpha = 100;
		}

		free(symmetryPar);
//		initReversibleGTR(tr, partitions, index);
//		if (setAlphaFreqs) {
//			makeGammaCats(current_part->alpha, current_part->gammaRates, 4,
//					tr->useMedian);
//		}
	} else {
		ProteicModel * pModel = static_cast<ProteicModel *>(model);
		current_part->dataType = PLL_AA_DATA;
		current_part->protFreqs = pModel->isPF();
		current_part->optimizeBaseFrequencies = PLL_FALSE;
		current_part->protModels = pModel->getMatrix();
		//initReversibleGTR(tr, partitions, index);
		if (pModel->isPF()) {
			memcpy(current_part->frequencies,
					current_part->empiricalFrequencies,
					NUM_AA_STATES * sizeof(double));
		}
		current_part->alpha = model->getAlpha();
	}
	//pllInitModel(tree, partitions, alignment->getPhylip());
}

} /* namespace partest */
