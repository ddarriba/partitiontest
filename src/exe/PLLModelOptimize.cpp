/*
 * PLLModelOptimize.cpp
 *
 *  Created on: Jan 8, 2013
 *      Author: diego
 */

#include "PLLModelOptimize.h"
#include <iostream>
#include <sstream>
//#if ! (defined(__ppc) || defined(__powerpc__) || defined(PPC))
//#include <xmmintrin.h>
//#endif
#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

extern "C" {
#include "globalVariables.h"
#include "newick.h"

#define MAX_FILE_LEN 255

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
	//pllTreeInitTopologyForAlignment(tree, phylip);

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
	cout << "[TRACE] PLLModelOptimize - Initializing model (str)" << endl;
#endif

	if (partitions->partitionData[0]->dataType == PLL_AA_DATA)
		for (int i=0; i<partitions->numberOfPartitions; i++)
			partitions->partitionData[i]->protModels = PLL_DAYHOFF;

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

	pllEvaluateGeneric(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);

	do {
		lk = tr->likelihood;

		// TODO: Optimize topology?

		if (estimateModel)
			pllOptimizeModelParameters(tr, partitions, 0.1);

		if (estimateBranchLengths)
			pllTreeEvaluate(tr, partitions, 200);

		pllEvaluateGeneric(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
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

char * PLLModelOptimize::getMlTree(PartitioningScheme * scheme,
		string inputFile) {

	char controlFilename[MAX_FILE_LEN];
	strcpy(controlFilename, tmpnam(NULL));
	char treeFilename[MAX_FILE_LEN];
	strcpy(treeFilename, tmpnam(NULL));

	string suffix(controlFilename);
	string::size_type pos = string( suffix ).find_last_of( "\\/" );
	suffix = suffix.substr( pos+1 );
	const char *cSuffix = suffix.c_str();

#ifdef _WIN32
	/* NOT SURE IF THIS WORKS */
	char currentDir[255];
	GetModuleFileName( NULL, currentDir, 255 );
	string::size_type pos = string( currentDir ).find_last_of( "\\/" );
#else
	char * currentDir = NULL;
	size_t tmpS;
	currentDir = getcwd(currentDir, tmpS);
#endif

	char outputTreeFilename[MAX_FILE_LEN];
	outputTreeFilename[0] = '\0';
	strcat(outputTreeFilename, options->getOutputTmpPath().c_str());
	strcat(outputTreeFilename, "RAxML_bestTree.");
	strcat(outputTreeFilename, cSuffix);

	char outputInfoFilename[MAX_FILE_LEN];
	outputInfoFilename[0] = '\0';
	strcat(outputInfoFilename, options->getOutputTmpPath().c_str());
	strcat(outputInfoFilename, "RAxML_info.");
	strcat(outputInfoFilename, cSuffix);

	char outputLogFilename[MAX_FILE_LEN];
	outputLogFilename[0] = '\0';
	strcat(outputLogFilename, options->getOutputTmpPath().c_str());
	strcat(outputLogFilename, "RAxML_log.");
	strcat(outputLogFilename, cSuffix);

	char outputResultsFilename[MAX_FILE_LEN];
	outputResultsFilename[0] = '\0';
	strcat(outputResultsFilename, options->getOutputTmpPath().c_str());
	strcat(outputResultsFilename, "RAxML_result.");
	strcat(outputResultsFilename, cSuffix);

	char outputParsimonyFilename[MAX_FILE_LEN];
	outputParsimonyFilename[0] = '\0';
	strcat(outputParsimonyFilename, options->getOutputTmpPath().c_str());
	strcat(outputParsimonyFilename, "RAxML_parsimonyTree.");
	strcat(outputParsimonyFilename, cSuffix);

	remove(outputTreeFilename);
	remove(outputInfoFilename);
	remove(outputLogFilename);
	remove(outputResultsFilename);
	remove(outputParsimonyFilename);

	ofstream controlFile;

	controlFile.open(controlFilename);
	for (int i = 0; i < scheme->getNumberOfElements(); i++) {
		PartitionElement * element = scheme->getElement(i);
		const char * name = element->getName().c_str();
		controlFile << "DNA," << name << "=";
		for (int j = 0; j < element->getNumberOfSections(); j++) {
			int start = element->getStart(j);
			int end = element->getEnd(j);
			controlFile << start << "-" << end;
			if (j < element->getNumberOfSections() - 1) {
				controlFile << ",";
			}
		}
		controlFile << endl;
	}
	controlFile.close();

	stringstream command;
	command << "bin/raxmlHPC -s " << options->getInputFile() << " -q "
			<< controlFilename << " -n " << cSuffix << " -m GTRGAMMA -p 1";

	if (options->getOutputTmpPath().length()) {
		if (options->getOutputTmpPath().c_str()[0] == char_separator) {
			command << " -w " << options->getOutputTmpPath();
		} else {
			char outputTmpFilename[MAX_FILE_LEN];
			outputTmpFilename[0] = '\0';
			strcat(outputTmpFilename, currentDir);
			strcat(outputTmpFilename, os_separator.c_str());
			strcat(outputTmpFilename, options->getOutputTmpPath().c_str());
			command << " -w " << outputTmpFilename;
		}
	}

	if (strlen(options->getTreeFile())) {
		command << " -t " << options->getTreeFile();
	} else if (options->getTreeString()) {
		ofstream stTreeFile;
		stTreeFile.open(treeFilename);
		stTreeFile << options->getTreeString() << endl;
		stTreeFile.close();
		command << " -t " << treeFilename;
	}
	cout << "COMMAND = " << command.str() << endl;

	if (system(command.str().c_str()) != 0) {
		cerr << "There was an error finding the ML topology. Check the logfiles ";
		if (options->getOutputTmpPath().length() > 0) {
			cerr << "(" << options->getOutputTmpPath() << ")" << endl;
		}
		cerr << endl;
		exit(0);
	}

	cout << "OUTPUT FILE = " << outputTreeFilename << endl;

	ifstream treeFile(outputTreeFilename);
	string line;
	if (treeFile.is_open()) {
		getline(treeFile, line);
		treeFile.close();
	}

	char * treeString = (char *) malloc(line.size() + 1);
	strcpy(treeString, line.c_str());

	int currentPartition = 0;
	ifstream resultsFile(outputInfoFilename);
	if (resultsFile.is_open()) {
		while (getline(resultsFile, line)) {
			if (!line.compare(0, 17, "Base frequencies:")) {
				double frequencies[4];
				frequencies[0] = atof(line.substr(18, 5).c_str());
				frequencies[1] = atof(line.substr(24, 5).c_str());
				frequencies[2] = atof(line.substr(30, 5).c_str());
				frequencies[3] = atof(line.substr(36, 5).c_str());
				scheme->getElement(currentPartition)->getModelset()->getModel(0)->setFrequencies(
						frequencies);
				currentPartition++;
			}
			if (!line.compare(0, 8, "alpha[0]")) {
				currentPartition = 0;
				istringstream iss(line);
				do {
					string token;
					iss >> token;
					if (token.compare(0, 5, "alpha")) {
						iss >> token;
						scheme->getElement(currentPartition)->getModelset()->getModel(
								0)->setAlpha(atof(token.c_str()));
					}
					double rates[6];
					if (token.compare(0, 5, "rates")) {
						for (int i = 0; i < 6; i++) {
							iss >> token;
						}
						for (int i = 0; i < 6; i++) {
							iss >> token;
							rates[i] = atof(token.c_str());
						}
					}
					scheme->getElement(currentPartition)->getModelset()->getModel(
							0)->setRates(rates);
					currentPartition++;
				} while (currentPartition < scheme->getNumberOfElements());
			}
		}
		resultsFile.close();
	}

	remove(controlFilename);
	remove(treeFilename);
	remove(outputTreeFilename);
	remove(outputInfoFilename);
	remove(outputLogFilename);
	remove(outputResultsFilename);
	remove(outputParsimonyFilename);

	return treeString;

}

double PLLModelOptimize::evaluateSPR(pllInstance * tr,
		partitionList *partitions, bool estimateModel, bool estimateTopology) {

	time_t t0 = time(NULL);

	pllRearrangeList * bestList = pllCreateRearrangeList(5);
	int i;
	double lk = PLL_UNLIKELY, prevLK = PLL_UNLIKELY, epsilon = 0.01, difference;

	tr->thoroughInsertion = PLL_FALSE;

	cout << "[SPR] " << Utilities::timeToString(time(NULL) - t0)
			<< " Initial Likelihood: " << tr->likelihood << endl;

	Tree2String(tr->tree_string, tr, partitions, tr->start->back, PLL_TRUE,
			PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH,
			PLL_FALSE, PLL_FALSE);

	cout << "[SPR] treeMP<- read.tree(text=\"" << tr->tree_string << "\")"
			<< endl;

	pllOptimizeModelParameters(tr, partitions, 1000);
	int bestMove = 0;
	int sprDistance = 2;

	do {
		difference = 0.0;

		pllRearrangeSearch(tr, partitions, PLL_REARRANGE_SPR,
				tr->nodep[tr->mxtips + 1], 1, sprDistance, bestList);

		for (i = 0; i < bestList->entries; i++) {
			pllRearrangeCommit(tr, partitions, &(bestList->rearr[i]), PLL_TRUE);
			pllTreeEvaluate(tr, partitions, 64);
			pllOptimizeModelParameters(tr, partitions, 10);
			pllEvaluateGeneric(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);

			cout << "[SPR] " << Utilities::timeToString(time(NULL) - t0)
					<< "LK " << tr->likelihood << endl;

			if (tr->likelihood > lk) {

				Tree2String(tr->tree_string, tr, partitions, tr->start->back,
						PLL_TRUE, PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE,
						PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);
				cout << "[SPR] treeSPR" << i << "<- read.tree(text=\""
						<< tr->tree_string << "\")" << endl;

				difference = tr->likelihood - lk;
				bestMove = i;
				lk = tr->likelihood;
			}

			pllRearrangeRollback(tr, partitions);
			pllEvaluateGeneric(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
		}
		cout << "[SPR] " << Utilities::timeToString(time(NULL) - t0) << " DIFF "
				<< difference << endl;
		sprDistance *= 2;

		if (difference > epsilon && sprDistance < (2 * tr->mxtips - 1)) {
			/* We don't need the rearrange list anymore */
			pllDestroyRearrangeList(&bestList);

			/* Now let's create another list and compute 30 rearrangement moves */
			bestList = pllCreateRearrangeList(5);
		}

	} while (difference > epsilon && sprDistance < (2 * tr->mxtips - 1));

	pllRearrangeCommit(tr, partitions, &(bestList->rearr[bestMove]), PLL_FALSE);
	pllTreeEvaluate(tr, partitions, 64);
	pllOptimizeModelParameters(tr, partitions, 10);
	pllEvaluateGeneric(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);

	tr->thoroughInsertion = PLL_TRUE;

	/* We don't need the rearrange list anymore */
	pllDestroyRearrangeList(&bestList);

	/* Now let's create another list and compute 30 rearrangement moves */
	bestList = pllCreateRearrangeList(1);

	pllRearrangeSearch(tr, partitions, PLL_REARRANGE_SPR,
			tr->nodep[tr->mxtips + 1], 1, 30, bestList);

	pllRearrangeCommit(tr, partitions, &(bestList->rearr[0]), PLL_TRUE);
	pllTreeEvaluate(tr, partitions, 64);
	pllOptimizeModelParameters(tr, partitions, 10);
	pllEvaluateGeneric(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);

	Tree2String(tr->tree_string, tr, partitions, tr->start->back, PLL_TRUE,
			PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH,
			PLL_FALSE, PLL_FALSE);

	cout << "[SPR] " << Utilities::timeToString(time(NULL) - t0) << "LK "
			<< tr->likelihood << endl;
	cout << "[SPR] treeLAST <- read.tree(text=\"" << tr->tree_string << "\")"
			<< endl;

	return tr->likelihood;
}

double prevSPR(pllInstance * tr, partitionList *partitions, bool estimateModel,
		bool estimateTopology) {

	pllRearrangeList * bestList = pllCreateRearrangeList(1);
	int i;
	double lk = PLL_UNLIKELY, epsilon = 0.01;

	time_t t0 = time(NULL);

#ifdef DEBUG
	cout << endl << "[SPR] " << Utilities::timeToString(time(NULL) - t0)
	<< " START" << endl;
#endif

	pllEvaluateGeneric(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);

#ifdef DEBUG
	cout << "[SPR] " << Utilities::timeToString(time(NULL) - t0)
	<< " Initial Likelihood: " << tr->likelihood << endl;

	Tree2String(tr->tree_string, tr, partitions, tr->start->back, PLL_TRUE,
			PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH,
			PLL_FALSE, PLL_FALSE);

	cout << "[SPR] FIRST TREE -> " << tr->tree_string;
#endif

	tr->thoroughInsertion = PLL_TRUE;
	double modelLkThreshold = 0.1;
	int SPRdistance = 2 * tr->mxtips - 2;
	int startNode = 1;

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

		pllTreeEvaluate(tr, partitions, 200);

		tr->thoroughInsertion = PLL_FALSE;

		pllEvaluateGeneric(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);

#ifdef DEBUG
		cout << "[SPR] " << Utilities::timeToString(time(NULL) - t0)
		<< " New Likelihood (TOPO): " << tr->likelihood << endl;
#endif

		if (lk > tr->likelihood) {
			if (estimateTopology) {
				pllRearrangeRollback(tr, partitions);
			}
			pllEvaluateGeneric(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);

#ifdef DEBUG
			cout << "[SPR] " << Utilities::timeToString(time(NULL) - t0)
			<< " Rollback (TOPO): " << tr->likelihood << endl;
#endif
			break;
		} else {
			Tree2String(tr->tree_string, tr, partitions, tr->start->back,
					PLL_TRUE, PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE,
					PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);
		}

		if (estimateModel) {

#ifdef DEBUG
			cout << "[SPR] " << Utilities::timeToString(time(NULL) - t0)
			<< " Optimizing model parameters..." << endl;
#endif

			pllOptimizeModelParameters(tr, partitions, modelLkThreshold);
			pllEvaluateGeneric(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);

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

	Tree2String(tr->tree_string, tr, partitions, tr->start->back, PLL_TRUE,
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
	partitionList *partitions = (partitionList *) malloc(sizeof(partitionList));
	char * tree = getMlTree(scheme, options->getInputFile());

	cout << "FINAL TREE: " << tree << endl << endl;

	return 0;
}

int oldPartitioningSchemeAtOnce(PartitioningScheme * scheme) {
#ifdef DEBUG
	cout << "[TRACE] PLLModelOptimize - Constructing partitions structure" << endl;
#endif

//	pllQueue * parts;
//	pllPartitionRegion * pregion;
//	pllPartitionInfo * pinfo;
//
//	pllQueueInit(&parts);
//
//	int numSeqs = alignment->getNumSeqs();
//	int numSites = alignment->getNumSites();
//	pllAlignmentData * phylip = pllInitAlignmentData(numSeqs, numSites);
//	phylip->sequenceCount = numSeqs;
//	phylip->sequenceLength = numSites;
//	for (int i = 0; i < numSeqs; i++) {
//		phylip->sequenceLabels[i + 1] = strdup(
//				alignment->getPhylip()->sequenceLabels[i + 1]);
//	}
//	phylip->siteWeights = (int *) malloc(numSites * sizeof(int));
//	int firstSite = 1;
//	int nextSite = 0;
//	for (int i = 0; i < scheme->getNumberOfElements(); i++) {
//		PartitionElement * element = scheme->getElement(i);
//
//		pinfo = (pllPartitionInfo *) malloc(sizeof(pllPartitionInfo));
//		pllQueueInit(&(pinfo->regionList));
//		pllQueueAppend(parts, (void *) pinfo);
//
//		pinfo->partitionName = (char *) malloc(
//				(element->getName().size() + 1) * sizeof(char));
//		strcpy(pinfo->partitionName, element->getName().c_str());
//		pinfo->partitionModel = (char *) malloc(5 * sizeof(char));
//		if (element->getBestModel()->getModel()->isPF()) {
//			strcpy(pinfo->partitionModel, "DNAX");
//			pinfo->optimizeBaseFrequencies = PLL_TRUE;
//		} else {
//			strcpy(pinfo->partitionModel, "DNA");
//			pinfo->optimizeBaseFrequencies = PLL_FALSE;
//		}
//		pinfo->protModels = -1;
//		pinfo->protFreqs = -1;
//		pinfo->dataType = PLL_DNA_DATA;
//
//		for (int j = 0; j < element->getNumberOfSections(); j++) {
//			for (int site = element->getStart(j) - 1; site < element->getEnd(j);
//					site++) {
//				phylip->siteWeights[nextSite] = 1;
//				for (int nk = 0; nk < numSeqs; nk++) {
//					phylip->sequenceData[nk + 1][nextSite] =
//							alignment->getPhylip()->sequenceData[nk + 1][site];
//				}
//				nextSite++;
//			}
//		}
//
//		pregion = (pllPartitionRegion *) malloc(sizeof(pllPartitionRegion));
//		pregion->start = firstSite;
//		pregion->end = nextSite;
//		pregion->stride = 1;
//		pllQueueAppend(pinfo->regionList, (void *) pregion);
//
//		firstSite = nextSite + 1;
//	}
//
//#ifdef DEBUG
//	cout << "[TRACE] PLLModelOptimize - Committing partitions" << endl;
//#endif
//	partitionList * partitions = pllPartitionsCommit(parts, phylip);
//
//	pllQueuePartitionsDestroy(&parts);
//
//	for (int i = 0; i < partitions->numberOfPartitions; i++) {
//		if (partitions->partitionData[i]->lower
//				== partitions->partitionData[i]->upper) {
//			exit(-1);
//		}
//	}
//
//	//pllPhylipRemoveDuplicate(phylip, partitions);
//
//#ifdef DEBUG
//	cout << "[TRACE] PLLModelOptimize - Creating tree instance" << endl;
//#endif
//	pllInstanceAttr * attr = (pllInstanceAttr *) rax_malloc(
//			sizeof(pllInstanceAttr));
//	attr->rateHetModel = PLL_GAMMA;
//	attr->fastScaling = PLL_FALSE;
//	attr->saveMemory = PLL_FALSE;
//	attr->useRecom = PLL_FALSE;
//	attr->randomNumberSeed = 12345;
//	attr->numberOfThreads = 1;
//	pllInstance * tr = pllCreateInstance(attr);
//	rax_free(attr);
//
//#ifdef DEBUG
//	cout << "[TRACE] PLLModelOptimize - Initializing topology" << endl;
//#endif
//	pllTreeInitTopologyForAlignment(tr, phylip);
//	/* Connect the alignment with the tree structure */
//
//#ifdef DEBUG
//	cout << "[TRACE] PLLModelOptimize - Connecting tree with alignment" << endl;
//#endif
//	if (!pllLoadAlignment(tr, phylip, partitions,
//	PLL_SHALLOW_COPY)) {
//		cerr << "ERROR: Incompatible tree/alignment combination" << endl;
//		Utilities::exit_partest(EX_SOFTWARE);
//	}
//
//#ifdef DEBUG
//	cout << "[TRACE] PLLModelOptimize - Initializing model (sch)" << endl;
//#endif
//	pllInitModel(tr, partitions, phylip);
//#ifdef DEBUG
//	cout << "[TRACE] PLLModelOptimize - Initialized model" << endl;
//#endif
//
//	/* set best-fit models */
//	for (int i = 0; i < scheme->getNumberOfElements(); i++) {
//		Model * model = scheme->getElement(i)->getBestModel()->getModel();
//		setModelParameters(model, tr, partitions, i);
//	}
//
//	if (options->getTreeString() == 0) {
//		pllComputeRandomizedStepwiseAdditionParsimonyTree(tr, partitions);
//	} else {
//		pllNewickTree * nt = pllNewickParseString(options->getTreeString());
//		pllTreeInitTopologyNewick(tr, nt, PLL_FALSE);
//		Tree2String(tr->tree_string, tr, partitions, tr->start->back, PLL_TRUE,
//				PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH,
//				PLL_FALSE, PLL_FALSE);
//	}
//
//	pllEvaluateGeneric(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
//	evaluateSPR(tr, partitions, options->getOptimizeMode() == OPT_GTR, false);
//
//	Tree2String(tr->tree_string, tr, partitions, tr->start->back, PLL_TRUE,
//			PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH,
//			PLL_FALSE, PLL_FALSE);
//	cout << "FINAL TREE: " << tr->tree_string << endl;
//
//	pllPartitionsDestroy(tr, &partitions);
//	partitions = 0;
//
//	pllDestroyInstance(tr);
//	tr = 0;

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
	pllEvaluateGeneric(tree, partitions, tree->start, PLL_TRUE, PLL_FALSE);
	optimizeParameters(tree, partitions, true, true, false);

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

void PLLModelOptimize::setModelParameters(Model * model, pllInstance * tr,
		partitionList * partitions, int index, bool setAlphaFreqs) {

	DataType dataType = options->getDataType();
	pInfo * current_part = partitions->partitionData[index];

	if (dataType == DT_NUCLEIC) {
		current_part->optimizeBaseFrequencies = model->isPF();
		current_part->alpha = model->getAlpha();
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
	} else {
		ProteicModel * pModel = static_cast<ProteicModel *>(model);
		current_part->dataType = PLL_AA_DATA;
		current_part->protFreqs = PLL_FALSE;
		current_part->optimizeBaseFrequencies = PLL_FALSE;
		current_part->protModels = pModel->getMatrix();
		initReversibleGTR(tr, partitions, index);
		if (pModel->isPF()) {
			memcpy(current_part->frequencies, current_part->empiricalFrequencies, NUM_AA_STATES*sizeof(double));
		}
	}
}

} /* namespace partest */
