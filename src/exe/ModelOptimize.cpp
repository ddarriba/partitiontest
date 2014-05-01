/*
 * ModelOptimize.cpp
 *
 *  Created on: Apr 9, 2014
 *      Author: diego
 */

#include "ModelOptimize.h"

#include "exe/ModelSelector.h"
#include "util/Utilities.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <pll.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <parsePartition.h>

using namespace std;

namespace partest {

ModelOptimize::ModelOptimize() {

}

ModelOptimize::~ModelOptimize() {

}

string ModelOptimize::buildStartingTree() {
	bool loadedTree = false;

	if (ckpAvailable) {
		fstream ofs((ckpPath + os_separator + ckpStartingTree).c_str(),
				ios::in);

		if (ofs) {
			cout << timestamp() << " Loading topology from checkpoint..."
					<< endl;
			ofs.seekg(0);
			int treeLen;
			ofs.read((char *) &(treeLen), sizeof(int));
			starting_tree = (char *) malloc(treeLen + 1);
			ofs.read((char *) starting_tree, treeLen);
			ofs.close();

			pllNewickTree * nt;
			nt = pllNewickParseString(starting_tree);
			pllTreeInitTopologyNewick(tree, nt, PLL_FALSE);
			pllNewickParseDestroy(&nt);
			strcpy(tree->tree_string, starting_tree);
			free(starting_tree);
			loadedTree = true;
		}
	}

	if (!loadedTree) {
		cout << timestamp() << " Computing fixed topology..." << endl;
		pllAlignmentData * alignData = 0;
		Utilities::duplicateAlignmentData(&alignData, phylip);

		partitionList * compParts = pllPartitionsCommit(pllPartsQueue,
				alignData);

		pllAlignmentRemoveDups(alignData, compParts);

		pllTreeInitTopologyForAlignment(tree, alignData);
		pllLoadAlignment(tree, alignData, compParts, PLL_DEEP_COPY);

		pllComputeRandomizedStepwiseAdditionParsimonyTree(tree, compParts);
		tree->start = tree->nodep[1];

		switch (data_type) {
		case DT_PROTEIC:
			for (int cur_part = 0; cur_part < compParts->numberOfPartitions;
					cur_part++) {
				pInfo * current_part = compParts->partitionData[cur_part];
				current_part->dataType = PLL_AA_DATA;
				current_part->states = 20;
				current_part->protFreqs = PLL_FALSE;
				current_part->optimizeBaseFrequencies = PLL_FALSE;
				current_part->optimizeAlphaParameter = PLL_TRUE;
				current_part->optimizeSubstitutionRates = PLL_FALSE;
				current_part->protModels = PLL_AUTO;
				current_part->alpha = 0.0;
				cout << timestamp() << " Loading AUTO models" << endl;
			}
			break;
		case DT_NUCLEIC:
			for (int cur_part = 0; cur_part < compParts->numberOfPartitions;
					cur_part++) {
				pInfo * current_part = compParts->partitionData[cur_part];
				current_part->dataType = PLL_DNA_DATA;
				current_part->states = 4;
				current_part->optimizeBaseFrequencies = PLL_TRUE;
				current_part->optimizeAlphaParameter = PLL_TRUE;
				current_part->optimizeSubstitutionRates = PLL_TRUE;
			}
			cout << timestamp() << " Loading GTR models" << endl;
			break;
		default:
			cerr << "Unknown datatype " << data_type << endl;
			exit_partest(EX_SOFTWARE);
		}

		pllInitModel(tree, compParts, alignData);

		cout << timestamp() << " Building ML topology" << endl;
		pllRaxmlSearchAlgorithm(tree, compParts, PLL_FALSE);

		pllTreeToNewick(tree->tree_string, tree, compParts, tree->start->back,
		PLL_TRUE,
		PLL_TRUE,
		PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH, PLL_FALSE,
		PLL_FALSE);
		tree->tree_string[tree->treeStringLength - 1] = '\0';

		pllPartitionsDestroy(tree, &compParts);
		pllAlignmentDataDestroy(alignData);

		if (ckpAvailable) {
			/* store tree */
			fstream ofs((ckpPath + os_separator + ckpStartingTree).c_str(),
					ios::out);
			ofs.seekg(0);
			int treeLen = strlen(tree->tree_string) + 1;
			ofs.write((char *) &treeLen, sizeof(int));
			ofs.write((char *) tree->tree_string, treeLen);
			ofs.close();
		}
	}

	starting_tree = tree->tree_string;
	cout << timestamp() << " Starting tree loaded " << starting_tree << endl;
	return string(tree->tree_string);
}

string ModelOptimize::buildFinalTreeLinking(PartitioningScheme * finalScheme,
		bool reoptimizeParameters) {
	cout << timestamp() << " Computing fixed topology..." << endl;
	pllAlignmentData * alignData = 0;
	Utilities::duplicateAlignmentData(&alignData, phylip);

	partitionList * compParts = pllPartitionsCommit(pllPartsQueue, alignData);

	pllAlignmentRemoveDups(alignData, compParts);

	pllTreeInitTopologyForAlignment(tree, alignData);
	pllLoadAlignment(tree, alignData, compParts, PLL_DEEP_COPY);

	pllComputeRandomizedStepwiseAdditionParsimonyTree(tree, compParts);
	tree->start = tree->nodep[1];

	for (unsigned int i = 0; i < finalScheme->getNumberOfElements(); i++) {
		PartitionElement * pe = finalScheme->getElement(i);
		switch (data_type) {
		case DT_PROTEIC:
			for (int cur_part = 0; cur_part < compParts->numberOfPartitions;
					cur_part++) {
				pe = finalScheme->getElement(cur_part);
				ProteicModel * pModel =
						static_cast<ProteicModel *>(pe->getBestModel()->getModel());
				int matrix = pModel->getMatrix();
				pInfo * current_part = compParts->partitionData[cur_part];
				current_part->dataType = PLL_AA_DATA;
				current_part->states = 20;
				current_part->protFreqs = pModel->isPF();
				current_part->optimizeBaseFrequencies = PLL_FALSE;
				current_part->optimizeAlphaParameter =
						reoptimizeParameters ? PLL_TRUE : PLL_FALSE;
				current_part->optimizeSubstitutionRates = PLL_FALSE;
				current_part->protModels = matrix;
				current_part->alpha = pModel->getAlpha();
			}
			break;
		case DT_NUCLEIC: {
			for (int cur_part = 0; cur_part < compParts->numberOfPartitions;
					cur_part++) {
				pe = finalScheme->getElement(cur_part);
				NucleicModel * nModel =
						static_cast<NucleicModel *>(pe->getBestModel()->getModel());
				//int matrix = nModel->getMatrix();
				pInfo * current_part = compParts->partitionData[cur_part];
				current_part->dataType = PLL_DNA_DATA;
				current_part->states = 4;
				current_part->optimizeBaseFrequencies = nModel->isPF();
				current_part->optimizeAlphaParameter =
						reoptimizeParameters ? PLL_TRUE : PLL_FALSE;
				current_part->optimizeSubstitutionRates =
						reoptimizeParameters ? PLL_TRUE : PLL_FALSE;
			}
			break;
		}
		default:
			exit_partest(EX_UNAVAILABLE);
		}
	}

	pllInitModel(tree, compParts, alignData);

	pllRaxmlSearchAlgorithm(tree, compParts, PLL_FALSE);

	pllTreeToNewick(tree->tree_string, tree, compParts, tree->start->back,
	PLL_TRUE,
	PLL_TRUE,
	PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH, PLL_FALSE,
	PLL_FALSE);
	tree->tree_string[tree->treeStringLength - 1] = '\0';

	pllPartitionsDestroy(tree, &compParts);
	pllAlignmentDataDestroy(alignData);

	starting_tree = tree->tree_string;
	cout << timestamp() << " Starting tree loaded" << endl << endl;
	return string(tree->tree_string);
}

string ModelOptimize::buildFinalTree(PartitioningScheme * finalScheme,
		bool reoptimizeParameters) {
	bool loadedTree = false;

	char * final_tree = 0;

	if (ckpAvailable) {
		fstream ofs((ckpPath + os_separator + ckpFinalTree).c_str(), ios::in);

		if (ofs) {
			cout << endl << timestamp()
					<< " Loading topology from checkpoint..." << endl;
			ofs.seekg(0);
			int treeLen;
			ofs.read((char *) &(treeLen), sizeof(int));
			final_tree = (char *) malloc(treeLen + 1);
			ofs.read((char *) final_tree, treeLen);
			ofs.close();

			loadedTree = true;
		}
	}

	if (!loadedTree) {
		PartitionElement * pe;

		pllInstanceAttr attr;
		attr.fastScaling = PLL_FALSE;
		attr.randomNumberSeed = 12345;
		attr.rateHetModel = PLL_GAMMA;
		attr.saveMemory = PLL_FALSE;
		attr.useRecom = PLL_FALSE;
		attr.numberOfThreads = number_of_threads;

		pllInstance * fTree = pllCreateInstance(&attr);

		pllQueue * parts;
		pllPartitionRegion * pregion;
		pllPartitionInfo * pinfo;

		pllQueueInit(&parts);
		for (unsigned int i = 0; i < finalScheme->getNumberOfElements(); i++) {
			pe = finalScheme->getElement(i);
			pinfo = (pllPartitionInfo *) malloc(sizeof(pllPartitionInfo));
			pllQueueInit(&(pinfo->regionList));
			pllQueueAppend(parts, (void *) pinfo);

			pinfo->partitionName = (char *) malloc(5);
			strcpy(pinfo->partitionName, "PART");
			pinfo->partitionModel = (char *) malloc(1);
			pinfo->protModels = -1;
			pinfo->protFreqs = -1;
			pinfo->dataType =
					(data_type == DT_NUCLEIC) ? PLL_DNA_DATA : PLL_AA_DATA;
			pinfo->optimizeBaseFrequencies = PLL_TRUE;
			for (int j = 0; j < pe->getNumberOfSections(); j++) {
				pregion = (pllPartitionRegion *) malloc(
						sizeof(pllPartitionRegion));
				pregion->start = pe->getSection(j).start;
				pregion->end = pe->getSection(j).end;
				pregion->stride = 1;
				pllQueueAppend(pinfo->regionList, (void *) pregion);
			}
		}
		partitionList * compParts = pllPartitionsCommit(parts, phylip);
		pllQueuePartitionsDestroy(&parts);

		cout << endl << timestamp() << " Computing final topology... " << endl;

		pllAlignmentRemoveDups(phylip, compParts);

		pllTreeInitTopologyForAlignment(fTree, phylip);
		pllLoadAlignment(fTree, phylip, compParts, PLL_SHALLOW_COPY);

		switch (starting_topology) {
		case StartTopoMP:
			pllComputeRandomizedStepwiseAdditionParsimonyTree(fTree, compParts);
			fTree->start = fTree->nodep[1];
			break;
		case StartTopoFIXED: {
			pllNewickTree * nt;
			nt = pllNewickParseString(starting_tree);
			pllTreeInitTopologyNewick(fTree, nt, PLL_FALSE);
			pllNewickParseDestroy(&nt);
			break;
		}
		case StartTopoUSER:
			cerr << "User Topo Not Available" << endl;
			exit_partest(EX_UNAVAILABLE);
			break;
		}

		switch (data_type) {
		case DT_PROTEIC:
			for (int cur_part = 0; cur_part < compParts->numberOfPartitions;
					cur_part++) {
				pe = finalScheme->getElement(cur_part);
				ProteicModel * pModel =
						static_cast<ProteicModel *>(pe->getBestModel()->getModel());
				int matrix = pModel->getMatrix();
				pInfo * current_part = compParts->partitionData[cur_part];
				current_part->dataType = PLL_AA_DATA;
				current_part->states = 20;
				current_part->protFreqs = pModel->isPF();
				current_part->optimizeBaseFrequencies = PLL_FALSE;
				current_part->optimizeAlphaParameter =
						reoptimizeParameters ? PLL_TRUE : PLL_FALSE;
				current_part->optimizeSubstitutionRates = PLL_FALSE;
				current_part->protModels = matrix;
				current_part->alpha = pModel->getAlpha();
				if (!reoptimizeParameters) {
					current_part->alpha = pModel->getAlpha();
					memcpy(current_part->frequencies, pModel->getFrequencies(),
					NUM_PROT_FREQS);
					memcpy(current_part->substRates, pModel->getRates(),
					NUM_AA_RATES);
				}
			}
			break;
		case DT_NUCLEIC: {
			for (int cur_part = 0; cur_part < compParts->numberOfPartitions;
					cur_part++) {
				pe = finalScheme->getElement(cur_part);
				NucleicModel * nModel =
						static_cast<NucleicModel *>(pe->getBestModel()->getModel());
				//int matrix = nModel->getMatrix();
				pInfo * current_part = compParts->partitionData[cur_part];
				current_part->dataType = PLL_DNA_DATA;
				current_part->states = 4;
				current_part->optimizeBaseFrequencies = nModel->isPF();
				current_part->optimizeAlphaParameter =
						reoptimizeParameters ? PLL_TRUE : PLL_FALSE;
				current_part->optimizeSubstitutionRates =
						reoptimizeParameters ? PLL_TRUE : PLL_FALSE;
				if (!reoptimizeParameters) {
					current_part->alpha = nModel->getAlpha();
					memcpy(current_part->frequencies, nModel->getFrequencies(),
					NUM_NUC_FREQS);
					memcpy(current_part->substRates, nModel->getRates(),
					NUM_DNA_RATES);
					if (!nModel->isPF()) {
						for (int i = 0; i < 4; i++) {
							current_part->frequencies[i] = 0.25;
						}
					}
				}
			}
			break;
		}
		default:
			exit_partest(EX_UNAVAILABLE);
		}

		pllInitModel(fTree, compParts, phylip);

		if (data_type == DT_NUCLEIC) {
			for (int cur_part = 0; cur_part < compParts->numberOfPartitions;
					cur_part++) {
				pe = finalScheme->getElement(cur_part);
				NucleicModel * nModel =
						static_cast<NucleicModel *>(pe->getBestModel()->getModel());
				const char * m = nModel->getMatrixName().c_str();
				char * symmetryPar = (char *) malloc(12 * sizeof(char));
				symmetryPar[0] = m[0];
				symmetryPar[11] = '\0';
				for (int j = 1; j < 6; j++) {
					symmetryPar[(j - 1) * 2 + 1] = ',';
					symmetryPar[j * 2] = m[j];
				}
				pllSetSubstitutionRateMatrixSymmetries(symmetryPar, compParts,
						cur_part);

			}
		}
		pllRaxmlSearchAlgorithm(fTree, compParts, PLL_FALSE);

		pllTreeToNewick(fTree->tree_string, fTree, compParts,
				fTree->start->back,
				PLL_TRUE,
				PLL_TRUE,
				PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH, PLL_FALSE,
				PLL_FALSE);
		fTree->tree_string[strlen(fTree->tree_string) - 1] = '\0';

		pllPartitionsDestroy(fTree, &compParts);

		int treeLen = strlen(fTree->tree_string) + 1;
		final_tree = (char *) malloc(treeLen);
		strcpy(final_tree, fTree->tree_string);

		if (ckpAvailable) {
			/* store tree */
			fstream ofs((ckpPath + os_separator + ckpFinalTree).c_str(),
					ios::out);
			ofs.seekg(0);

			ofs.write((char *) &treeLen, sizeof(int));
			ofs.write((char *) final_tree, treeLen);
			ofs.close();
		}

		cout << timestamp() << " Final tree lnL: " << fTree->likelihood << endl;
		pllDestroyInstance(fTree);
	}

	cout << timestamp() << " Final tree: " << final_tree << endl;

	return string(final_tree);
}

int ModelOptimize::optimizePartitioningScheme(PartitioningScheme * scheme,
		int index, int limit) {

	cout << timestamp() << " - scheme " << setw(Utilities::iDecLog(limit) + 1)
			<< setfill('0') << right << index + 1 << "/" << limit
			<< setfill(' ') << endl;

	/* check number of elements to optimize */
	int elementsToOptimize = 0;
	for (unsigned int cur_element = 0;
			cur_element < scheme->getNumberOfElements(); cur_element++) {
		if (!scheme->getElement(cur_element)->isOptimized()) {
			elementsToOptimize++;
		}
	}

	int currentElementToOptimize = 0;
	for (unsigned int cur_element = 0;
			cur_element < scheme->getNumberOfElements(); cur_element++) {
		PartitionElement * element = scheme->getElement(cur_element);
		if (!element->isOptimized()) {
			optimizePartitionElement(element, currentElementToOptimize,
					elementsToOptimize);
			currentElementToOptimize++;
		}
	}

	return EX_OK;
}

int ModelOptimize::optimizePartitionElement(PartitionElement * element,
		int index, int limit) {

	if (element->isOptimized()) {
		return EX_OK;
	}

	cout << timestamp() << " - - - element "
			<< setw(Utilities::iDecLog(limit) + 1) << setfill('0') << right
			<< index + 1 << "/" << limit << setfill(' ') << endl;
	element->setupStructures();

	for (int modelIndex = 0; modelIndex < element->getNumberOfModels();
			modelIndex++) {
		optimizeModel(element, modelIndex, element->getNumberOfModels());
	}

	ModelSelector ms(element, BIC, element->getSampleSize());
	//ms.print(cout);

	element->destroyStructures();

	return EX_OK;
}

void ModelOptimize::setModelParameters(Model * _model, pllInstance * _tree,
		partitionList * _partitions, pllAlignmentData * _alignData, int index,
		bool setAlphaFreqs) {

	pInfo * current_part = _partitions->partitionData[index];

	if (data_type == DT_NUCLEIC) {
		const char * m = _model->getMatrixName().c_str();
		char * symmetryPar = (char *) malloc(12 * sizeof(char));
		symmetryPar[0] = m[0];
		symmetryPar[11] = '\0';
		for (int j = 1; j < 6; j++) {
			symmetryPar[(j - 1) * 2 + 1] = ',';
			symmetryPar[j * 2] = m[j];
		}
		pllSetSubstitutionRateMatrixSymmetries(symmetryPar, _partitions, index);

		current_part->optimizeBaseFrequencies = _model->isPF();
		if (!_model->isPF()) {
			for (int i = 0; i < 4; i++) {
				current_part->frequencies[i] = 0.25;
			}
		}
		current_part->alpha = _model->getAlpha();
		//current_part->nonGTR = PLL_FALSE;
		current_part->dataType = PLL_DNA_DATA;

		pllInitReversibleGTR(_tree, _partitions, index);
		if (setAlphaFreqs) {
			memcpy(current_part->frequencies, _model->getFrequencies(),
					4 * sizeof(double));
			memcpy(current_part->substRates, _model->getRates(),
					6 * sizeof(double));
			current_part->alpha = _model->getAlpha();
		} else {
			current_part->optimizeBaseFrequencies = _model->isPF();
			if (!_model->isPF()) {
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
	} else {
		ProteicModel * pModel = static_cast<ProteicModel *>(_model);
		current_part->dataType = PLL_AA_DATA;
		current_part->protFreqs = pModel->isPF();
		current_part->optimizeBaseFrequencies = PLL_FALSE;
		current_part->protModels = pModel->getMatrix();
		current_part->alpha = pModel->getAlpha();
	}

	//TODO: This works if partitions has one single partitions
	assert(_partitions->numberOfPartitions == 1);
	pllInitReversibleGTR(_tree, _partitions, 0);

	double **ef = pllBaseFrequenciesGTR(_partitions, _alignData);
	//initModel(_tree, ef, _partitions);
	free(*ef);
	free(ef);
}

void ModelOptimize::optimizeModel(PartitionElement * element,
		unsigned int modelIndex, int limit) {

	pllInstance * _tree = element->getTree();
	partitionList * _partitions = element->getPartitions();
	pllAlignmentData * _alignData = element->getAlignData();
	Model * model = element->getModel(modelIndex);

	/* set parameters for single partition element */
	setModelParameters(model, _tree, _partitions, _alignData, 0, false);\
	double lk;
	double epsilon = 0.1;

	_tree->thoroughInsertion = PLL_FALSE;

	pllEvaluateLikelihood(_tree, _partitions, _tree->start, PLL_TRUE,
	PLL_FALSE);

	/* main optimization loop */
	do {
		lk = _tree->likelihood;
		pllOptimizeModelParameters(_tree, _partitions, 1);
		pllOptimizeBranchLengths(_tree, _partitions, 10);
		pllEvaluateLikelihood(_tree, _partitions, _tree->start, PLL_TRUE,
		PLL_FALSE);
	} while (fabs(lk - _tree->likelihood) > epsilon);

	/* construct newick tree for optimized model */
	pllTreeToNewick(_tree->tree_string, _tree, _partitions, _tree->start->back,
	PLL_TRUE, PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE,
	PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);
	_tree->tree_string[strlen(_tree->tree_string) - 1] = '\0';
	model->setLnL(_tree->likelihood);
	model->setTree(_tree->tree_string);

	model->setFrequencies(_partitions->partitionData[0]->frequencies);
	if (model->isGamma())
		model->setAlpha(_partitions->partitionData[0]->alpha);
	model->setRates(_partitions->partitionData[0]->substRates);

	if (data_type == DT_PROTEIC && optimize_mode == OPT_GTR) {
		/* set chosen model */
		model->setName(
				Utilities::getProtMatrixName(
						static_cast<ProtMatrix>(_partitions->partitionData[0]->autoProtModels)));
	}

	cout << timestamp() << " - - - - - " << setw(Utilities::iDecLog(limit) + 1)
			<< setfill('0') << right << modelIndex + 1 << "/" << limit << " "
			<< model->getName() << " (" << _tree->likelihood << ")"
			<< setfill(' ') << endl;
}

} /* namespace partest */
