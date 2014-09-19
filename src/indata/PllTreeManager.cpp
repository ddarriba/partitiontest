/*
 * PLLTreeManager.cpp
 *
 *  Created on: Sep 18, 2014
 *      Author: diego
 */

#include "PllTreeManager.h"
#include "util/Utilities.h"
#include "indata/PartitionMap.h"
#include "model/NucleicModel.h"
#include "model/ProteicModel.h"

#include <cassert>
#include <cstring>
#include <iostream>

#include <parsePartition.h>

#define ESTIMATE_PARAMETERS 0

namespace partest {

pllInstance * buildTree() {
	pllInstanceAttr attr;
	attr.fastScaling = PLL_FALSE;
	attr.randomNumberSeed = 12345;
	attr.rateHetModel = PLL_GAMMA;
	attr.saveMemory = PLL_FALSE;
	attr.useRecom = PLL_FALSE;
	attr.numberOfThreads = number_of_threads;
	return (pllCreateInstance(&attr));
}

PllTreeManager::PllTreeManager(const t_partitionElementId id,
		const pllAlignmentData * phylip, const vector<PEsection> & sections,
		size_t numberOfSites) :
		TreeManager(id, numberOfSites, numberOfSites) {

	_tree = buildTree();
	_alignData = pllInitAlignmentData(phylip->sequenceCount, numberOfSites);
	_alignData->siteWeights = (int *) malloc(numberOfSites * sizeof(int));
	for (int seq = 1; seq <= _alignData->sequenceCount; seq++) {
		_alignData->sequenceLabels[seq] = (char *) malloc(
				strlen(phylip->sequenceLabels[seq]) + 1);
		strcpy(_alignData->sequenceLabels[seq], phylip->sequenceLabels[seq]);
		int nextSite = 0;
		for (size_t i = 0; i < sections.size(); i++) {
			size_t part = sections[i].id;
			int lower = pllPartitions->partitionData[part]->lower;
			int width = pllPartitions->partitionData[part]->width;
			memcpy(&(_alignData->sequenceData[seq][nextSite]),
					&(phylip->sequenceData[seq][lower]),
					width * sizeof(unsigned char));
			nextSite += width;
		}
	}

	for (size_t site = 0; site < numberOfSites; site++) {
		_alignData->siteWeights[site] = 1;
	}
	pllQueue * partsQueue;
	pllQueueInit(&partsQueue);

	pllPartitionInfo * pinfo = (pllPartitionInfo *) malloc(
			sizeof(pllPartitionInfo));
	pllQueueInit(&(pinfo->regionList));
	pinfo->partitionModel = (char *) malloc(1);
	pinfo->protModels = -1;
	pinfo->protUseEmpiricalFreqs = -1;
	pinfo->dataType = data_type == DT_NUCLEIC ? PLL_DNA_DATA : PLL_AA_DATA;
	pinfo->optimizeBaseFrequencies = PLL_TRUE;
	pinfo->partitionName = (char *) malloc(8 * sizeof(char));
	pinfo->ascBias = PLL_FALSE;
	strcpy(pinfo->partitionName, "NewGene");

	int nextStart = 1;
	for (int part = 0; part <= pllPartitions->numberOfPartitions; part++) {
		for (size_t i = 0; i < sections.size(); i++) {
			if (sections[i].id == part) {
				pllPartitionRegion * pregion = (pllPartitionRegion *) malloc(
						sizeof(pllPartitionRegion));
				pregion->start = nextStart;
				pregion->end = nextStart
						+ pllPartitions->partitionData[part]->width - 1;
				pregion->stride = 1;
				nextStart = pregion->end + 1;
				pllQueueAppend(pinfo->regionList, (void *) pregion);
				break;
			}
		}
	}
	pllQueueAppend(partsQueue, (void *) pinfo);

	_partitions = pllPartitionsCommit(partsQueue, _alignData);

	pllQueuePartitionsDestroy(&partsQueue);

	assert(_alignData->sequenceLength == (int ) numberOfSites);
	pllAlignmentRemoveDups(_alignData, _partitions);
	numberOfPatterns = _alignData->sequenceLength;

	pllTreeInitTopologyForAlignment(_tree, _alignData);
	pllLoadAlignment(_tree, _alignData, _partitions);

	switch (starting_topology) {
	case StartTopoML:
		/* set ML optimization parameters */
		_tree->doCutoff = ML_PARAM_CUTOFF;
		if (epsilon == AUTO_EPSILON) {
			_tree->likelihoodEpsilon = 0.1;
		} else {
			_tree->likelihoodEpsilon = epsilon;
		}
		_tree->stepwidth = ML_PARAM_STEPWIDTH;
		_tree->max_rearrange = ML_PARAM_MAXREARRANGE;
		_tree->initial = tree->bestTrav = ML_PARAM_BESTTRAV;
		_tree->initialSet = ML_PARAM_INITIALSET;
		pllComputeRandomizedStepwiseAdditionParsimonyTree(_tree, _partitions);
		_tree->start = _tree->nodep[1];
		break;
	case StartTopoMP:
		pllComputeRandomizedStepwiseAdditionParsimonyTree(_tree, _partitions);
		_tree->start = _tree->nodep[1];
		break;
	case StartTopoFIXED: {
		pllNewickTree * nt;
		if (pergene_branch_lengths) {
			if (pergene_branch_lengths && sections.size() == 1) {
				nt = pllNewickParseString(
						pergene_starting_tree[sections[0].id]);
			} else {
				t_partitionElementId id(sections.size());
				for (size_t i = 0; i < sections.size(); i++) {
					id[i] = sections[i].id;
				}
				nt = Utilities::averageBranchLengths(id);
			}
		} else {
			nt = pllNewickParseString(starting_tree);
		}
		_tree->fracchange = 1.0;
		pllTreeInitTopologyNewick(_tree, nt, PLL_FALSE);
		pllNewickParseDestroy(&nt);
		break;
	}
	case StartTopoUSER:
		cerr << "User Topo Not Available" << endl;
		exit_partest(EX_UNAVAILABLE);
		break;
	default:
		cerr << "ERROR: Undefined starting topology" << endl;
		exit_partest(EX_SOFTWARE);
		break;
	}

	pllInitModel(_tree, _partitions);
}

PllTreeManager::~PllTreeManager() {
	if (_tree) {
		if (_partitions) {
			pllPartitionsDestroy(_tree, &_partitions);
			_partitions = 0;
		}
		pllDestroyInstance(_tree);
		_tree = 0;
	}
	if (_alignData) {
		pllAlignmentDataDestroy(_alignData);
		_alignData = 0;
	}
}

double * PllTreeManager::getBranchLengths(void) {
	double * bls = (double *) malloc(
			(size_t) Utilities::numberOfBranches(num_taxa) * sizeof(double));
	for (int i = 0; i < Utilities::numberOfBranches(num_taxa); i++) {
		if (isnan(_tree->nodep[i + 1]->z[0])) {
			cout << "ERROR: NAN branch in " << i + 1 << endl;
			exit_partest(EX_SOFTWARE);
		}
		bls[i] = _tree->nodep[i + 1]->z[0];
	}
	return bls;
}

void PllTreeManager::setModelParameters(Model * _model, int index,
		bool setAlphaFreqs) {

	pInfo * current_part = _partitions->partitionData[index];
	current_part->ascBias = PLL_FALSE;

	if (data_type == DT_NUCLEIC) {
		const char * m = _model->getMatrixName().c_str();
		char * symmetryPar = (char *) malloc(12 * sizeof(char));
		symmetryPar[0] = m[0];
		symmetryPar[11] = '\0';
		for (int j = 1; j < NUM_DNA_RATES; j++) {
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
		current_part->dataType = PLL_DNA_DATA;
		pllInitReversibleGTR(_tree, _partitions, index);
		if (setAlphaFreqs || _id.size() == 1) {
			current_part->optimizeBaseFrequencies = _model->isPF();
			memcpy(current_part->frequencies, _model->getFrequencies(),
					4 * sizeof(double));
			for (int i = 0; i < NUM_NUC_FREQS; i++) {
				current_part->freqExponents[i] = 1.0;
			}
			memcpy(current_part->substRates, _model->getRates(),
			NUM_DNA_RATES * sizeof(double));

			current_part->alpha = _model->getAlpha();
			pllMakeGammaCats(current_part->alpha, current_part->gammaRates, 4,
					_tree->useMedian);
		} else {
			PartitionMap * partitionMap = PartitionMap::getInstance();
			if (ESTIMATE_PARAMETERS) {
				double alpha = 0.0;
				double frequencies[NUM_NUC_FREQS];
				for (int j = 0; j < NUM_NUC_FREQS; j++) {
					frequencies[j] = 0.25;
				}
				double rates[NUM_DNA_RATES];
				for (int j = 0; j < NUM_DNA_RATES; j++) {
					rates[j] = 0.0;
				}
				for (size_t i = 0; i < _id.size(); i++) {
					t_partitionElementId singleId(1);
					singleId.at(0) = _id.at(i);
					PartitionElement * element =
							partitionMap->getPartitionElement(singleId);
//					double * elementFreqs =
//							element->getBestModel()->getModel()->getFrequencies();
//					for (int j = 0; j < NUM_NUC_FREQS; j++) {
//						frequencies[j] += elementFreqs[j] / id.size();
//					}
					double * elementRates =
							element->getBestModel()->getModel()->getRates();
					for (int j = 0; j < NUM_DNA_RATES; j++) {
						rates[j] += elementRates[j] / _id.size();
					}
					alpha += element->getBestModel()->getModel()->getAlpha()
							/ _id.size();
				}

//				if (!_model->isPF()) {
//					for (int i = 0; i < NUM_NUC_FREQS; i++) {
//						frequencies[i] = 0.25;
//					}
//				}

				pllSetFixedBaseFrequencies(frequencies, NUM_NUC_FREQS, index,
						_partitions, _tree);
				current_part->optimizeBaseFrequencies = _model->isPF();

				if (strcmp(_model->getMatrixName().c_str(), "012345")) {
					for (int i = 0; i < NUM_DNA_RATES; i++) {
						rates[i] = 1;
					}
				}
				pllSetFixedSubstitutionMatrix(rates, NUM_DNA_RATES, index,
						_partitions, _tree);
				current_part->optimizeSubstitutionRates = PLL_TRUE;

				current_part->alpha = alpha;
				current_part->optimizeAlphaParameter = PLL_TRUE;
				pllMakeGammaCats(current_part->alpha, current_part->gammaRates,
						4, _tree->useMedian);

			} else {
				if (!_model->isPF()) {
					for (int i = 0; i < NUM_NUC_FREQS; i++) {
						current_part->freqExponents[i] = 1.0;
						current_part->frequencies[i] = 0.25;
					}
				}
				for (int i = 0; i < NUM_DNA_RATES; i++) {
					current_part->substRates[i] = 1;
				}
				current_part->alpha = 1.0;
				current_part->optimizeBaseFrequencies = _model->isPF();
				pllMakeGammaCats(current_part->alpha, current_part->gammaRates,
						4, _tree->useMedian);
			}
		}

		free(symmetryPar);
	} else {
		ProteicModel * pModel = static_cast<ProteicModel *>(_model);
		current_part->dataType = PLL_AA_DATA;
		current_part->protUseEmpiricalFreqs = pModel->isPF();
		current_part->optimizeBaseFrequencies = PLL_FALSE;
		current_part->optimizeSubstitutionRates = PLL_FALSE;
		current_part->optimizeAlphaParameter = PLL_TRUE;
		current_part->protModels = pModel->getMatrix();
		current_part->alpha = 1.0;
		for (int i = 0; i < NUM_PROT_FREQS; i++) {
			current_part->freqExponents[i] = 1.0;
		}
		pllMakeGammaCats(current_part->alpha, current_part->gammaRates, 4,
				_tree->useMedian);

	}

	pllInitReversibleGTR(_tree, _partitions, index);

	_tree->thoroughInsertion = PLL_FALSE;
	pllEvaluateLikelihood(_tree, _partitions, _tree->start, PLL_TRUE,
	PLL_FALSE);
}

double PllTreeManager::searchMlTopology(bool estimateModel) {
	pllRaxmlSearchAlgorithm(_tree, _partitions, estimateModel);
	return _tree->likelihood;
}

double PllTreeManager::evaluateLikelihood(bool fullTraversal) {
	pllEvaluateLikelihood(_tree, _partitions, _tree->start, fullTraversal,
	PLL_FALSE);
	return _tree->likelihood;
}

double PllTreeManager::getLikelihood() {
	return _tree->likelihood;
}

double * PllTreeManager::getFrequencies(size_t partition) {
	return _partitions->partitionData[partition]->frequencies;
}

double * PllTreeManager::getRates(size_t partition) {
	return _partitions->partitionData[partition]->substRates;
}

double PllTreeManager::getAlpha(size_t partition) {
	return _partitions->partitionData[partition]->alpha;
}

int PllTreeManager::getAutoProtModel(size_t partition) {
	return _partitions->partitionData[partition]->autoProtModels;
}

void PllTreeManager::optimizeBranchLengths(int smoothIterations) {
	pllOptimizeBranchLengths(_tree, _partitions, smoothIterations);
}

void PllTreeManager::optimizeBaseFreqs(double epsilon) {
	pllOptBaseFreqs(_tree, _partitions, epsilon, _partitions->freqList);
}

void PllTreeManager::optimizeRates(double epsilon) {
	pllOptRatesGeneric(_tree, _partitions, epsilon, _partitions->rateList);
}

void PllTreeManager::optimizeAlphas(double epsilon) {
	pllOptAlphasGeneric(_tree, _partitions, epsilon, _partitions->alphaList);
}

const char * PllTreeManager::getNewickTree() {
	pllTreeToNewick(_tree->tree_string, _tree, _partitions, _tree->start->back,
	PLL_TRUE, PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE,
	PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);
	return _tree->tree_string;
}

void PllTreeManager::optimizeModelParameters(double epsilon) {
	pllOptimizeModelParameters(_tree, _partitions, epsilon);
}

}

