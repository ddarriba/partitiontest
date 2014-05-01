/*
 * PartitionElement.cpp
 *
 *  Created on: Apr 8, 2014
 *      Author: diego
 */

#include "PartitionElement.h"

#include "util/Utilities.h"
#include <parsePartition.h>

#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <assert.h>
#include <fstream>
#include <sstream>

using namespace std;

namespace partest {

PartitionElement::PartitionElement(t_partitionElementId id) :
		id(id) {

	this->bestModel = 0;
	models.reserve(number_of_models);
	numberOfSections = id.size();
	numberOfSites = 0;
	numberOfPatterns = 0;
	sections = (PEsection *) malloc(numberOfSections * sizeof(PEsection));
	name.append("(");
	for (unsigned int i = 0; i < numberOfSections; i++) {
		unsigned int part = id.at(i);

		sections[i].start = pllPartitions->partitionData[part]->lower + 1;
		sections[i].end = pllPartitions->partitionData[part]->upper;
		sections[i].id = part;
		name.append(*(singleGeneNames[part]));
		if (i < numberOfSections - 1)
			name.append(",");
		//this->id.push_back(part);
		numberOfSites += pllPartitions->partitionData[part]->width;
	}
	name.append(")");

	stringstream ss;
	/* prefix */
	ss << "pt_";
	unsigned int maxId = id.at(id.size() - 1);
	int numchars = ceil((maxId + 1) / 6.0);
	int curPos = 0;
	for (int currentChar = numchars - 1; currentChar >= 0; currentChar--) {
		unsigned int charRangeStart = 6 * (numchars - currentChar - 1);
		unsigned int charRangeEnd = 6 * (numchars - currentChar);
		int value = 0;
		for (unsigned int k = charRangeStart; k < charRangeEnd; k++) {
			if (id[curPos] == k) {
				value += Utilities::binaryPow(k - charRangeStart);
				curPos++;
				if (k == maxId)
					break;
			}
		}
		ss << Utilities::toBase64(value);
	}
	string sstr = ss.str();
	if (sstr.length() < MAX_FILE_LENGTH) {
		ckpname = ss.str();
		ckphash = "";
	} else {
		ckpname = sstr.substr(0, MAX_FILE_LENGTH);
		ckphash = sstr.substr(MAX_FILE_LENGTH, MAX_FILE_LENGTH - sstr.length());
	}

	sampleSize = numberOfSites;
	loadData();
}

int PartitionElement::setupStructures(void) {
	if (!isOptimized()) {
		pllInstanceAttr attr;
		attr.fastScaling = PLL_FALSE;
		attr.randomNumberSeed = 12345;
		attr.rateHetModel = PLL_GAMMA;
		attr.saveMemory = PLL_FALSE;
		attr.useRecom = PLL_FALSE;
		attr.numberOfThreads = number_of_threads;
		_tree = pllCreateInstance(&attr);

		_alignData = pllInitAlignmentData(phylip->sequenceCount, numberOfSites);
		_alignData->siteWeights = (int *) malloc(numberOfSites * sizeof(int));
		for (int seq = 1; seq <= _alignData->sequenceCount; seq++) {
			_alignData->sequenceLabels[seq] = (char *) malloc(
					strlen(phylip->sequenceLabels[seq]) + 1);
			strcpy(_alignData->sequenceLabels[seq],
					phylip->sequenceLabels[seq]);
			int nextSite = 0;
			for (unsigned int i = 0; i < numberOfSections; i++) {
				unsigned int part = id.at(i);
				int lower = pllPartitions->partitionData[part]->lower;
				int width = pllPartitions->partitionData[part]->width;
				memcpy(&(_alignData->sequenceData[seq][nextSite]),
						&(phylip->sequenceData[seq][lower]),
						width * sizeof(unsigned char));
				nextSite += width;
			}
		}
		for (unsigned int site = 0; site < numberOfSites; site++) {
			_alignData->siteWeights[site] = 1;
		}
		pllQueue * partsQueue;
		pllQueueInit(&partsQueue);

		pllPartitionInfo * pinfo = (pllPartitionInfo *) malloc(
				sizeof(pllPartitionInfo));
		pllQueueInit(&(pinfo->regionList));
		pinfo->partitionModel = (char *) malloc(1);
		pinfo->protModels = -1;
		pinfo->protFreqs = -1;
		pinfo->dataType = PLL_DNA_DATA;
		pinfo->optimizeBaseFrequencies = PLL_TRUE;
		pinfo->partitionName = (char *) malloc(8 * sizeof(char));
		strcpy(pinfo->partitionName, "NewGene");

		int nextStart = 1;
		for (int part = 0; part <= pllPartitions->numberOfPartitions; part++) {
			if (find(id.begin(), id.end(), part) != id.end()) {
				pllPartitionRegion * pregion = (pllPartitionRegion *) malloc(
						sizeof(pllPartitionRegion));
				pregion->start = nextStart;
				pregion->end = nextStart
						+ pllPartitions->partitionData[part]->width - 1;
				pregion->stride = 1;
				nextStart = pregion->end + 1;
				pllQueueAppend(pinfo->regionList, (void *) pregion);
			}
		}
		pllQueueAppend(partsQueue, (void *) pinfo);

		_partitions = pllPartitionsCommit(partsQueue, _alignData);

		pllQueuePartitionsDestroy(&partsQueue);

		assert(_alignData->sequenceLength == (int) numberOfSites);
		pllAlignmentRemoveDups(_alignData, _partitions);
		numberOfPatterns = _alignData->sequenceLength;

		switch (data_type) {
		case DT_PROTEIC:
			for (int cur_part = 0; cur_part < _partitions->numberOfPartitions;
					cur_part++) {
				pInfo * current_part = _partitions->partitionData[cur_part];
				current_part->dataType = PLL_AA_DATA;
				current_part->states = 20;
				current_part->protFreqs = PLL_FALSE;
				current_part->optimizeBaseFrequencies = PLL_FALSE;
				current_part->protModels = PLL_DAYHOFF;
				current_part->alpha = 0.0;
			}
			break;
		case DT_NUCLEIC:
			for (int cur_part = 0; cur_part < _partitions->numberOfPartitions;
								cur_part++) {
				pInfo * current_part = _partitions->partitionData[cur_part];
				current_part->dataType = PLL_DNA_DATA;
				current_part->states = 4;
			}
			break;
		default:
			exit_partest(EX_UNAVAILABLE);
		}

		pllTreeInitTopologyForAlignment(_tree, _alignData);
		pllLoadAlignment(_tree, _alignData, _partitions, PLL_DEEP_COPY);

		pllNewickTree * nt;
		switch (starting_topology) {
		case StartTopoMP:
			pllComputeRandomizedStepwiseAdditionParsimonyTree(_tree,
					_partitions);
			_tree->start = _tree->nodep[1];
			break;
		case StartTopoFIXED:
			nt = pllNewickParseString(starting_tree);
			pllTreeInitTopologyNewick(_tree, nt, PLL_FALSE);
			pllNewickParseDestroy(&nt);
			break;
		case StartTopoUSER:
			cerr << "User Topo Not Available" << endl;
			exit_partest(EX_UNAVAILABLE);
			break;
		}

		pllInitModel(_tree, _partitions, _alignData);

		if (models.size() == 0) {
			/* build model set */
			switch (data_type) {
			case DT_NUCLEIC:
				switch (optimize_mode) {
				case OPT_GTR:
					models.push_back(
							new NucleicModel(NUC_MATRIX_GTR, RateVarG | RateVarF,
									num_taxa));
					break;
				case OPT_SEARCH:
					for (int current_model = 0; current_model < NUC_MATRIX_SIZE;
							current_model += 2) {
						models.push_back(
								new NucleicModel(
										static_cast<NucMatrix>(current_model),
										RateVarG, num_taxa));
						if (do_rate & RateVarF) {
							models.push_back(
									new NucleicModel(
											static_cast<NucMatrix>(current_model+1),
											RateVarG | RateVarF, num_taxa));
						}
					}
					break;
				case OPT_CUSTOM:
				case OPT_DEFAULT:
					exit_partest(EX_SOFTWARE);
				}
				break;
			case DT_PROTEIC:
				switch (optimize_mode) {
				case OPT_SEARCH:
				case OPT_CUSTOM:
					for (int current_model = 0;
							current_model < PROT_MATRIX_SIZE; current_model++) {
						if (Utilities::binaryPow(current_model) & protModels) {
							models.push_back(
									new ProteicModel(
											static_cast<ProtMatrix>(current_model),
											RateVarG, num_taxa));
							if (do_rate & RateVarF) {
								models.push_back(
										new ProteicModel(
												static_cast<ProtMatrix>(current_model),
												RateVarG | RateVarF, num_taxa));
							}
						}
					}
					break;
				case OPT_GTR:
					models.push_back(
							new ProteicModel(PROT_MATRIX_AUTO, RateVarG,
									num_taxa));
					break;
				case OPT_DEFAULT:
					exit_partest(EX_SOFTWARE);
				}
				break;
			default:
				exit_partest(EX_SOFTWARE);
			}
		}
	}

	ready = true;
	return EX_OK;
}

int PartitionElement::destroyStructures(void) {

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

	ready = false;

	return EX_OK;
}

PartitionElement::~PartitionElement() {

	free (sections);

	for (Model * model : models) {
		delete model;
	}
	delete bestModel;

	if (_tree) {
		if (_partitions) {
			pllPartitionsDestroy(_tree, &_partitions);
		}
		pllDestroyInstance(_tree);
	}
	if (_alignData) {
		pllAlignmentDataDestroy(_alignData);
	}
}

pllAlignmentData * PartitionElement::getAlignData(void) {
	if (!ready) {
		cerr << "[ERROR] Alignment data is not ready" << endl;
		exit_partest(EX_SOFTWARE);
	}
	return _alignData;
}

pllInstance * PartitionElement::getTree(void) {
	if (!ready) {
		cerr << "[ERROR] Tree is not ready" << endl;
		exit_partest(EX_SOFTWARE);
	}
	return _tree;
}

partitionList * PartitionElement::getPartitions(void) {
	if (!ready) {
		cerr << "[ERROR] Partitions structure is not ready" << endl;
		exit_partest(EX_SOFTWARE);
	}
	return _partitions;
}

vector<Model *> PartitionElement::getModels(void) const {
	if (models.size() == 0) {
		cerr << "[ERROR] Models are not ready" << endl;
		exit_partest(EX_SOFTWARE);
	}
	return models;
}

int PartitionElement::getNumberOfModels(void) const {
	return models.size();
}

int PartitionElement::getNumberOfSites(void) const {
	return numberOfSites;
}

int PartitionElement::getNumberOfSections(void) const {
	return numberOfSections;
}

PEsection PartitionElement::getSection(unsigned int i) {
	if (i >= numberOfSections) {
		exit_partest(EX_SOFTWARE);
	}
	return sections[i];
}

int PartitionElement::getNumberOfPatterns(void) const {
	return numberOfPatterns;
}

Model * PartitionElement::getModel(unsigned int index) {
	if (index >= models.size()) {
		cerr << "[ERROR] Model " << index << " is above the number of models ("
				<< models.size() << ")" << endl;
		exit_partest(EX_SOFTWARE);
	}
	return models.at(index);
}

void PartitionElement::setBestModel(SelectionModel * model) {
	if (bestModel) {
		cerr << "[ERROR] BestModel was already set." << endl;
		exit_partest(EX_SOFTWARE);
	}
	this->bestModel = model->clone();

	if (!ckpLoaded) {
		storeData();
	}
}

SelectionModel * PartitionElement::getBestModel(void) {
	if (!bestModel) {
		cerr << "[ERROR] BestModel has not been set." << endl;
		exit_partest(EX_SOFTWARE);
	}
	return bestModel;
}

double PartitionElement::getLnL(void) const {
	return bestModel->getModel()->getLnL();
}

double PartitionElement::getSampleSize(void) {
	return sampleSize;
}

bool PartitionElement::isReady(void) {
	return ready;
}

bool PartitionElement::isOptimized(void) {
	return (models.size() > 0) && (models.at(0)->isOptimized());
}

/* checkpointing stuff */

int PartitionElement::loadData(void) {
	if (!ckpAvailable)
		return CHECKPOINT_UNAVAILABLE;

	fstream ofs((ckpPath + os_separator + ckpname).c_str(), ios::in);
	ofs.seekg(0, ios_base::beg);

	if (!ofs)
		return CHECKPOINT_UNEXISTENT;

	streampos curpos = 0;
	while (!ofs.eof() && !ckpLoaded) {
		ofs.seekg(curpos, ios_base::cur);
		int ckpSize;
		ofs.read((char *) &(ckpSize), sizeof(int));

		int hashlen;
		ofs.read((char *) &hashlen, sizeof(int));
		if (hashlen > 0) {
			char * charhash = (char *) malloc(hashlen + 1);
			ofs.read((char *) charhash, hashlen);
			charhash[hashlen] = '\0';
			if (!strcmp(charhash, ckphash.c_str())) {
				ckpLoaded = true;
			}
		} else if (ckphash == "") {
			ckpLoaded = true;
		}

		if (!ckpLoaded) {
			if (ckpSize <= 0) {
				cerr << "ERROR LOADING CHECKPOINTS" << endl;
				exit_partest(EX_SOFTWARE);
			}
			curpos = ckpSize - hashlen - sizeof(int);
		}
	}
	if (ckpLoaded) {
		ofs.read((char *) &(numberOfSections), sizeof(int));
		ofs.read((char *) &(numberOfSites), sizeof(int));
		ofs.read((char *) &(numberOfPatterns), sizeof(int));

		size_t modelSize =
				data_type == DT_NUCLEIC ?
						sizeof(NucleicModel) : sizeof(ProteicModel);
		for (unsigned int i = 0; i < number_of_models; i++) {
			Model * model = 0;
			switch (data_type) {
			case DT_NUCLEIC:
				model = (NucleicModel *) alloca(modelSize);
				break;
			case DT_PROTEIC:
				model = (ProteicModel *) alloca(modelSize);
				break;
			case DT_DEFAULT:
			default:
				cerr << "[INTERNAL_ERROR] Unrecognized datatype" << endl;
				exit_partest(EX_SOFTWARE);
			}
			ofs.read((char *) model, modelSize);
			double * freqs = (double *) alloca(
					model->getNumberOfFrequencies() * sizeof(double));
			ofs.read((char *) freqs,
					model->getNumberOfFrequencies() * sizeof(double));

			/* checked */

			double * rates = 0;
			if (data_type == DT_NUCLEIC) {
				rates = (double *) alloca(NUM_RATES * sizeof(double));
				ofs.read((char *) rates, NUM_RATES * sizeof(double));
				//model->setRates(rates);
			}
			int len_tree;
			ofs.read((char *) &len_tree, sizeof(size_t));
			char ctree[len_tree];
			ofs.read((char *) &ctree, len_tree);
			if (data_type == DT_NUCLEIC) {
				NucleicModel * finalModel = new NucleicModel(
						((NucleicModel *) model)->getMatrix(),
						model->getRateVariation(), num_taxa);
				finalModel->setFrequencies(freqs);
				finalModel->setRates(rates);
				if (model->isGamma())
					finalModel->setAlpha(model->getAlpha());
				if (model->isPInv())
					finalModel->setpInv(model->getpInv());
				finalModel->setLnL(model->getLnL());
				finalModel->setTree(ctree);
				models.push_back(finalModel);
			} else {
				ProteicModel * finalModel = new ProteicModel(
						((ProteicModel *) model)->getMatrix(),
						model->getRateVariation(), num_taxa);
				finalModel->setFrequencies(freqs);
				if (model->isGamma())
					finalModel->setAlpha(model->getAlpha());
				if (model->isPInv())
					finalModel->setpInv(model->getpInv());
				finalModel->setLnL(model->getLnL());
				finalModel->setTree(ctree);
				models.push_back(finalModel);
			}
		}
		int bestModelIndex;
		ofs.read((char *) &bestModelIndex, sizeof(int));

		SelectionModel * sm = (SelectionModel *) alloca(sizeof(SelectionModel));
		ofs.read((char *) sm, sizeof(SelectionModel));
		SelectionModel * selectionmodel = new SelectionModel(
				models.at(bestModelIndex), sm->getValue());
		selectionmodel->setWeight(sm->getWeight());
		selectionmodel->setCumWeight(sm->getCumWeight());
		selectionmodel->setDelta(sm->getDelta());
		selectionmodel->setIndex(bestModelIndex);
		setBestModel(selectionmodel);
	}

	ofs.close();
	return ckpLoaded ? CHECKPOINT_LOADED : CHECKPOINT_UNEXISTENT;
}

int PartitionElement::storeData(void) {
	if (!ckpAvailable)
		return CHECKPOINT_UNAVAILABLE;

	if (!isOptimized()) {
		cerr
				<< "[INTERNAL_ERROR] Attempting to save unoptimized Partition Element"
				<< endl;
		exit_partest(EX_SOFTWARE);
	}

	fstream ofs((ckpPath + os_separator + ckpname).c_str(),
			ios::in | ios::out | ios::app);

	//ofs.open((ckpPath + os_separator + ckpname).c_str());
	size_t modelSize =
			data_type == DT_NUCLEIC ?
					sizeof(NucleicModel) : sizeof(ProteicModel);
	int numberOfFrequencies = data_type == DT_NUCLEIC ?
	NUM_NUC_FREQS :
														NUM_PROT_FREQS;
	int numberOfRates = data_type == DT_NUCLEIC ?
	NUM_DNA_RATES :
													0;
	streampos startPos = ofs.tellp();
	int hashlen = ckphash.length();
	int ckpSize = 4 * sizeof(int) + hashlen
			+ models.size()
					* (modelSize + sizeof(size_t) + tree->treeStringLength
							+ (numberOfFrequencies + numberOfRates)
									* sizeof(double)) + sizeof(int)
			+ sizeof(SelectionModel);
	ofs.write((char *) &ckpSize, sizeof(int));
	ofs.write((char *) &hashlen, sizeof(int));
	if (hashlen > 0) {
		ofs.write((char *) ckphash.c_str(), hashlen);
	}
	ofs.write((char *) &numberOfSections, sizeof(int));
	ofs.write((char *) &numberOfSites, sizeof(int));
	ofs.write((char *) &numberOfPatterns, sizeof(int));
	int bestModelIndex = 0;
	for (unsigned int i = 0; i < models.size(); i++) {
		Model * model = models.at(i);
		if (model == getBestModel()->getModel()) {
			bestModelIndex = i;
		}
		model->setFrequencies(model->getFrequencies());
		ofs.write((char *) model, modelSize);
		ofs.write((char *) model->getFrequencies(),
				model->getNumberOfFrequencies() * sizeof(double));
		if (data_type == DT_NUCLEIC) {
			ofs.write((char *) model->getRates(), NUM_RATES * sizeof(double));
		}
		int name_len, matrixname_len, tree_len;
		name_len = model->getName().length() + 1;
		matrixname_len = model->getMatrixName().length() + 1;
		tree_len = strlen(model->getTree().c_str()) + 1;
		ofs.write((char *) &tree_len, sizeof(size_t));
		ofs.write((char *) model->getTree().c_str(), tree_len);
	}
	/* best model */
	SelectionModel * selectionmodel = getBestModel();
	ofs.write((char *) &bestModelIndex, sizeof(int));
	ofs.write((char *) selectionmodel, sizeof(SelectionModel));
	ofs.close();

	return CHECKPOINT_SAVED;
}

void PartitionElement::print(ostream & out) {
	cout << name << endl;
	cout << "Best model: " << bestModel->getModel()->getName() << endl;
	cout << "Num.Params: " << bestModel->getModel()->getNumberOfFreeParameters()
			<< endl;
	cout << "BIC score:  " << bestModel->getValue() << endl;
	cout << "BIC weight: " << bestModel->getWeight() << endl;
	bestModel->getModel()->print(out, "  ");
}

}

