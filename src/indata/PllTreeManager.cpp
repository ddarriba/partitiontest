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

#include <pll/parsePartition.h>
#include <cassert>
#include <cstring>
#include <iostream>
#include <iomanip>


#define ESTIMATE_PARAMETERS 0

using namespace std;

namespace partest {

static pllInstance * buildTree() {
	pllInstanceAttr attr;
	attr.fastScaling = PLL_FALSE;
	attr.randomNumberSeed = 0x54321;
	attr.rateHetModel = PLL_GAMMA;
	attr.saveMemory = PLL_FALSE;
	attr.useRecom = PLL_FALSE;
	attr.numberOfThreads = number_of_threads;
	return (pllCreateInstance(&attr));
}

PllTreeManager::PllTreeManager(const t_partitionElementId id,
		const pllAlignmentData * _phylip, const vector<PEsection> & sections,
		size_t numberOfSites) :
		TreeManager(id, numberOfSites, numberOfSites) {

	_tree = buildTree();
	_alignData = pllInitAlignmentData(_phylip->sequenceCount, (int) numberOfSites);
	_alignData->siteWeights = (int *) malloc(numberOfSites * sizeof(int));
	for (int seq = 1; seq <= _alignData->sequenceCount; seq++) {
		_alignData->sequenceLabels[seq] = (char *) malloc(
				strlen(_phylip->sequenceLabels[seq]) + 1);
		strcpy(_alignData->sequenceLabels[seq], _phylip->sequenceLabels[seq]);
		int nextSite = 0;
		for (size_t i = 0; i < sections.size(); i++) {
			size_t part = sections[i].id;
			size_t lower = (size_t) pllPartitions->partitionData[part]->lower;
			size_t width = (size_t) pllPartitions->partitionData[part]->width;
			memcpy(&(_alignData->sequenceData[seq][nextSite]),
					&(_phylip->sequenceData[seq][lower]),
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
			if (sections[i].id == (size_t) part) {
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
	numberOfPatterns = (size_t) _alignData->sequenceLength;

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
			if (sections.size() == 1) {
				nt = pllNewickParseString(
						pergene_starting_tree[sections[0].id]);
			} else {
				t_partitionElementId nextId(sections.size());
				for (size_t i = 0; i < sections.size(); i++) {
					nextId[i] = sections[i].id;
				}
				nt = Utilities::averageBranchLengths(nextId);
			}
		} else {
			nt = pllNewickParseString(starting_tree);
		}
		pllTreeInitTopologyNewick(_tree, nt, PLL_FALSE);
		pllNewickParseDestroy(&nt);
		break;
	}
	case StartTopoUSER:
		pllNewickTree * nt;
		if (pergene_branch_lengths) {
			if (sections.size() == 1) {
				nt = pllNewickParseFile(user_tree->c_str());
			} else {
				t_partitionElementId nextI(sections.size());
				for (size_t i = 0; i < sections.size(); i++) {
					nextI[i] = sections[i].id;
				}
				nt = Utilities::averageBranchLengths(nextI);
			}
		} else {
			nt = pllNewickParseFile(user_tree->c_str());
		}
		pllTreeInitTopologyNewick(_tree, nt, PLL_FALSE);
		pllNewickParseDestroy(&nt);
		break;
	default:
		assert(0);
		break;
	}

	if (data_type == DT_PROTEIC)
			_partitions->partitionData[0]->protModels = PLL_AUTO;
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
			(size_t) Utilities::numberOfBranches((int)num_taxa) * sizeof(double));
	for (int i = 0; i < Utilities::numberOfBranches((int)num_taxa); i++) {
		assert(!isnan(_tree->nodep[i + 1]->z[0]));
		bls[i] = _tree->nodep[i + 1]->z[0];
	}
	return bls;
}

void PllTreeManager::setModelParameters(const Model * _model, int index,
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
			#if(ESTIMATE_PARAMETERS)
				PartitionMap * partitionMap = PartitionMap::getInstance();
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

			#else
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
			#endif
		}

		free(symmetryPar);
	} else {
		const ProteicModel * pModel = static_cast<const ProteicModel *>(_model);
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

		if (pModel->isPF()) {
			double ** freqs = pllBaseFrequenciesInstance(_tree,_partitions);
			Utilities::smoothFrequencies(freqs[index], NUM_PROT_FREQS);
			memcpy(current_part->empiricalFrequencies, freqs[index], 20*sizeof(double));
			free(freqs);
		}
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

#define MIN_MULTIPLIER 0.001

void PllTreeManager::scaleBranchLengthsSymmetric( int smoothIterations ) {
	double blScaler = branchLengthMultiplier;
	double epsMultiplier = 0.5;
	double lkUpper,lkLower;
	double lkIni = getLikelihood();
	double lkEnd = -DOUBLE_INF;
	for (int i=0; i<smoothIterations; i++) {
		lkUpper = scaleBranchLengths(blScaler + epsMultiplier);
		lkLower = scaleBranchLengths(max(blScaler - epsMultiplier, MIN_MULTIPLIER));
		if (lkUpper > lkLower) {
			lkEnd = lkUpper;
			blScaler = blScaler + epsMultiplier;
		} else {
			lkEnd = lkLower;
			blScaler = blScaler - epsMultiplier;
		}
		epsMultiplier /= 2;
	}
	if (lkEnd > lkIni) {
		branchLengthMultiplier = blScaler;
	}
	lkEnd = scaleBranchLengths(branchLengthMultiplier);
}

void PllTreeManager::optimizeBranchLengths(int smoothIterations) {
	if (reoptimize_branch_lengths) {
		pllOptimizeBranchLengths(_tree, _partitions, smoothIterations);
	} else {
		scaleBranchLengthsSymmetric(smoothIterations);
	}
}

void PllTreeManager::optimizeBaseFreqs(double _epsilon) {
	pllOptBaseFreqs(_tree, _partitions, _epsilon, _partitions->freqList);
}

void PllTreeManager::optimizeRates(double _epsilon) {
	pllOptRatesGeneric(_tree, _partitions, _epsilon, _partitions->rateList);
}

void PllTreeManager::optimizeAlphas(double _epsilon) {
	pllOptAlphasGeneric(_tree, _partitions, _epsilon, _partitions->alphaList);
}

const char * PllTreeManager::getNewickTree() {
	pllTreeToNewick(_tree->tree_string, _tree, _partitions, _tree->start->back,
	PLL_TRUE, PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE,
	PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);
	return _tree->tree_string;
}

void PllTreeManager::optimizeModelParameters(double _epsilon) {
	pllOptimizeModelParameters(_tree, _partitions, _epsilon);
}

static double fixZ(double z) {
	if (z > PLL_ZMAX)
		return PLL_ZMAX;

	if (z < PLL_ZMIN)
		return PLL_ZMIN;

	return z;
}

double PllTreeManager::scaleBranchLengths(double multiplier) {
	int nodes = _tree->mxtips + _tree->mxtips - 2;
	assert(_partitions->numberOfPartitions == 1);
	size_t count;

	evaluateLikelihood(PLL_TRUE);

	if (storedBranchLengths.size() == 0) {
		/* store original branch lengths */
		storedBranchLengths.resize(((2 * (size_t)_tree->mxtips - 3) * 2));
		count = 0;
		for (int i = 1; i <= nodes; i++) {
			storedBranchLengths[count] = -log(_tree->nodep[i]->z[0]);
			count++;
			if (i > _tree->mxtips) {
				storedBranchLengths[count] = -log(_tree->nodep[i]->next->z[0]);
				count++;
				storedBranchLengths[count] = -log(_tree->nodep[i]->next->next->z[0]);
				count++;
			}
		}
		assert(count == (2 * (size_t) _tree->mxtips - 3) * 2);
	}

		count = 0;
		for (int i = 1; i <= nodes; i++) {
			double z = multiplier * storedBranchLengths[count];
			z = exp(-z);
			z = fixZ(z);
			_tree->nodep[i]->z[0] = z;

			count++;

			if (i > _tree->mxtips) {
				z = multiplier * storedBranchLengths[count];
				z = exp(-z);
				z = fixZ(z);
				_tree->nodep[i]->next->z[0] = z;

				count++;

				z = multiplier * storedBranchLengths[count];
				z = exp(-z);
				z = fixZ(z);
				_tree->nodep[i]->next->next->z[0] = z;

				count++;
			}
		}
		assert(count == (2 * (size_t) _tree->mxtips - 3) * 2);
		evaluateLikelihood(PLL_TRUE);
		return getLikelihood();
}

}
