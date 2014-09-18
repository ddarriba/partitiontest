/*
 * PLLTreeManager.cpp
 *
 *  Created on: Sep 18, 2014
 *      Author: diego
 */

#include "PllTreeManager.h"
#include "util/Utilities.h"

#include <cassert>
#include <cstring>
#include <iostream>

#include <parsePartition.h>

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

PllTreeManager::PllTreeManager(const pllAlignmentData * phylip,
		const vector<PEsection> & sections, size_t numberOfSites) :
		TreeManager(numberOfSites, numberOfSites) {

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

//vector<double> PllTreeManager::getBranchLengths(void) {
//	if (starting_topology == StartTopoFIXED && (branchLengths.size() == 0)
//			&& _tree) {
//		/* if topology is fixed for each element, keep the branch lengths */
//		branchLengths.resize((size_t) Utilities::numberOfBranches(num_taxa));
//		for (int i = 0; i < Utilities::numberOfBranches(num_taxa); i++) {
//			if (isnan(_tree->nodep[i + 1]->z[0])) {
//				cout << "ERROR: NAN branch in " << i + 1 << endl;
//				exit_partest(EX_SOFTWARE);
//			}
//			branchLengths[i] = _tree->nodep[i + 1]->z[0];
//		}
//	}
//	return branchLengths;
//}
//
//void PllTreeManager::setBranchLengths(double * values) {
//	cerr << "CKP SET"  << endl;
//	cerr << "CKP SIZE " <<  branchLengths.size() << endl;
//	cerr << "CKP RESIZE " << (size_t) Utilities::numberOfBranches(num_taxa) << endl;
//	branchLengths.resize((size_t) Utilities::numberOfBranches(num_taxa));
//	cerr << "CKP MEMCPY? " << endl;
////	for (int i = 0; i < Utilities::numberOfBranches(num_taxa); i++) {
////		branchLengths[i] = values[i];
////	}
//	memcpy(&(branchLengths[0]), values,
//			(size_t) Utilities::numberOfBranches(num_taxa) * sizeof(double));
//	cerr << "CKP DONE "  << endl;
//}

}

