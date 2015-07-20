/*  PartitionTest, fast selection of the best fit partitioning scheme for
 *  multi-gene data sets.
 *  Copyright May 2013 by Diego Darriba
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 3 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 *
 *  For any other inquiries send an Email to Diego Darriba
 *  ddarriba@udc.es
 */

#include "PartitionElement.h"

#include "util/Utilities.h"

#include <pll/parsePartition.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <assert.h>
#include <fstream>
#include <sstream>
#include <assert.h>

using namespace std;

namespace partest {

PartitionElement::PartitionElement(t_partitionElementId _id) :
		ready(false), id(_id), sampleSize(0.0), treeManager(0), sections(
				id.size()), ckpLoaded(false), tag(false), branchLengths(0) {

	this->bestModel = 0;
	models.reserve(number_of_models);
	numberOfSections = id.size();
	numberOfSites = 0;
	numberOfPatterns = 0;
	name.append("(");
	for (size_t i = 0; i < numberOfSections; i++) {
		size_t part = id.at(i);

		sections[i].start = (size_t) pllPartitions->partitionData[part]->lower + 1;
		sections[i].end = (size_t) pllPartitions->partitionData[part]->upper;
		sections[i].id = part;
		name.append(*(singleGeneNames[part]));
		if (i < numberOfSections - 1)
			name.append(",");
		//this->id.push_back(part);
		numberOfSites += (size_t) pllPartitions->partitionData[part]->width;
	}
	name.append(")");

	stringstream ss;
	/* prefix */
	ss << "pt_";
	size_t maxId = id.at(id.size() - 1);
	size_t numchars = (size_t) ceil((maxId + 1) / 6.0);
	size_t curPos = 0;
	for (long currentChar = (long) numchars - 1; currentChar >= 0; currentChar--) {
		size_t charRangeStart = 6 * (numchars - (size_t)currentChar - 1);
		size_t charRangeEnd = 6 * (numchars - (size_t)currentChar);
		int value = 0;
		for (size_t k = charRangeStart; k < charRangeEnd; k++) {
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
		treeManager = new PllTreeManager(id, phylip, sections, numberOfSites);
		numberOfPatterns = treeManager->getNumberOfPatterns();

		if (models.size() == 0) {
			/* build model set */
			switch (data_type) {
			case DT_NUCLEIC:
				switch (optimize_mode) {
				case OPT_GTR:
					models.push_back(
							new NucleicModel(NUC_MATRIX_GTR,
									RateVarG | RateVarF, (int) num_taxa));
					break;
				case OPT_SEARCH:
					for (int current_model = 0; current_model < NUC_MATRIX_SIZE;
							current_model += 2) {
						models.push_back(
								new NucleicModel(
										static_cast<NucMatrix>(current_model),
										RateVarG, (int) num_taxa));
						if (do_rate & RateVarF) {
							models.push_back(
									new NucleicModel(
											static_cast<NucMatrix>(current_model
													+ 1), RateVarG | RateVarF,
											(int) num_taxa));
						}
					}
					break;
				case OPT_CUSTOM:
					assert(0);
				}
				break;
			case DT_PROTEIC:
				switch (optimize_mode) {
				case OPT_SEARCH:
				case OPT_CUSTOM:
					for (size_t current_model = 0;
							current_model < PROT_MATRIX_SIZE; current_model++) {
						if (Utilities::binaryPow(current_model) & protModels) {
							models.push_back(
									new ProteicModel(
											static_cast<ProtMatrix>(current_model),
											RateVarG, (int) num_taxa));
							if (do_rate & RateVarF) {
								models.push_back(
										new ProteicModel(
												static_cast<ProtMatrix>(current_model),
												RateVarG | RateVarF, (int) num_taxa));
							}
						}
					}
					break;
				case OPT_GTR:
					models.push_back(
							new ProteicModel(PROT_MATRIX_AUTO, RateVarG,
									(int) num_taxa));
					break;
				default:
					assert(0);
				}
				break;
			default:
				assert(0);
			}
		}
	}

	ready = true;
	return EX_OK;
}

int PartitionElement::destroyStructures(void) {

	if (starting_topology == StartTopoFIXED && !branchLengths) {
		/* keep branch lengths */
		branchLengths = treeManager->getBranchLengths();
	}
	delete treeManager;
	treeManager = 0;

	ready = false;

	return EX_OK;

}

PartitionElement::~PartitionElement() {

	if (treeManager)
		delete treeManager;

	if (branchLengths)
		free(branchLengths);

	for (size_t i = 0; i < models.size(); i++) {
		Model * model = models.at(i);
		delete model;
	}
	delete bestModel;

}

TreeManager * PartitionElement::getTreeManager(void) {
	if (!ready) {
			cerr << "[ERROR] Alignment data is not ready" << endl;
			exit_partest(EX_SOFTWARE);
		}
	return treeManager;
}

double * PartitionElement::getBranchLengths(void) {
	if (starting_topology == StartTopoFIXED && !branchLengths && ready) {
		/* if topology is fixed for each element, keep the branch lengths */
		branchLengths = treeManager->getBranchLengths();
	}
	return branchLengths;
}

vector<Model *> PartitionElement::getModels(void) const {
	if (models.size() == 0) {
		cerr << "[ERROR] Models are not ready" << endl;
		exit_partest(EX_SOFTWARE);
	}
	return models;
}

size_t PartitionElement::getNumberOfModels(void) const {
	return models.size();
}

size_t PartitionElement::getNumberOfSites(void) const {
	return numberOfSites;
}

size_t PartitionElement::getNumberOfSections(void) const {
	return numberOfSections;
}

PEsection PartitionElement::getSection(size_t i) {
	if (i >= numberOfSections) {
		exit_partest(EX_SOFTWARE);
	}
	return sections[i];
}

size_t PartitionElement::getNumberOfPatterns(void) const {
	return numberOfPatterns;
}

Model * PartitionElement::getModel(size_t index) {
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

	const char * ckpFilename = (ckpPath + os_separator + ckpname).c_str();
	fstream ofs(ckpFilename, ios::in);
	ofs.seekg(0, ios_base::beg);

	if (!ofs)
		return CHECKPOINT_UNEXISTENT;

	streampos curpos = 0;
	while (!ofs.eof() && !ckpLoaded) {
		ofs.seekg(curpos, ios_base::cur);
		size_t ckpSize;
		ofs.read((char *) &(ckpSize), (streamsize) sizeof(size_t));

		size_t hashlen;
		ofs.read((char *) &hashlen, (streamsize) sizeof(size_t));
		if (hashlen > 0) {
			char * charhash = (char *) malloc(hashlen + 1);
			ofs.read((char *) charhash, (streamsize) hashlen);
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
			curpos = (streampos)ckpSize - (streampos)hashlen - (streampos)sizeof(size_t);
		}
	}

	if (ckpLoaded) {
		ofs.read((char *) &(numberOfSections), (streamsize) sizeof(size_t));
		ofs.read((char *) &(numberOfSites), (streamsize) sizeof(size_t));
		ofs.read((char *) &(numberOfPatterns), (streamsize) sizeof(size_t));
		size_t ckpNumberOfModels;
		ofs.read((char *) &(ckpNumberOfModels), (streamsize) sizeof(size_t));

		if (ckpNumberOfModels != number_of_models) {
			if (!force_overriding) {
				cerr << "[ERROR]   ";
			} else {
				cerr << "[WARNING] ";
			}
			cerr
					<< "Number of models in checkpoint files does not match the expected number"
					<< endl;
			cerr
					<< "          of models. Configuration might changed. Checkpoint loading for this "
					<< endl;
			cerr    << number_of_models << " instead of " << ckpNumberOfModels << endl;
			cerr << "          element was aborted." << endl;
			if (force_overriding) {
				ckpLoaded = false;
				ofs.close();
				/* deleting checkpoint file */
				if (remove(ckpFilename)) {
					cerr
							<< "[ERROR] There was an error removing checkointing file "
							<< ckpFilename
							<< ". Please remove it manually and try again."
							<< endl;
					exit_partest(EX_IOERR);
				}
			} else {
				cerr << endl
						<< "If you really want to proceed, remove checkpointing files or re-run with --force-override, --disable-output or --disable-ckp arguments."
						<< endl << endl;
				exit_partest(EX_IOERR);
			}
			return CHECKPOINT_UNEXISTENT;
		}
		size_t modelSize =
				data_type == DT_NUCLEIC ?
						sizeof(NucleicModel) : sizeof(ProteicModel);
		for (size_t i = 0; i < number_of_models; i++) {
			Model * model = 0;
			switch (data_type) {
			case DT_NUCLEIC:
				model = (NucleicModel *) alloca(modelSize);
				break;
			case DT_PROTEIC:
				model = (ProteicModel *) alloca(modelSize);
				break;
			default:
				assert(0);
			}

			ofs.read((char *) model, (streamsize) modelSize);
			double * freqs = (double *) alloca(
					(size_t)model->getNumberOfFrequencies() * sizeof(double));
			ofs.read((char *) freqs,
					(streamsize) ((size_t)model->getNumberOfFrequencies() * sizeof(double)));

			double * rates = 0;
			if (data_type == DT_NUCLEIC) {
				rates = (double *) alloca(NUM_DNA_RATES * sizeof(double));
				ofs.read((char *) rates, (streamsize) (NUM_DNA_RATES * sizeof(double)));
				//model->setRates(rates);
			}
			size_t len_tree;
			ofs.read((char *) &len_tree, sizeof(size_t));

			char * ctree = (char *) malloc(len_tree);

			ofs.read((char *) ctree, (streamsize) len_tree);
			if (data_type == DT_NUCLEIC) {
				NucleicModel * finalModel = new NucleicModel(
						((NucleicModel *) model)->getMatrix(),
						model->getRateVariation(), (int) num_taxa);
				finalModel->setFrequencies(freqs);
				finalModel->setRates(rates);
				if (model->isGamma())
					finalModel->setAlpha(model->getAlpha());
#ifdef _IG_MODELS
				if (model->isPInv())
				finalModel->setpInv(model->getpInv());
#endif
				finalModel->setLnL(model->getLnL());
				finalModel->setTree(ctree);
				models.push_back(finalModel);
			} else {
				ProteicModel * finalModel = new ProteicModel(
						((ProteicModel *) model)->getMatrix(),
						model->getRateVariation(), (int) num_taxa);
				finalModel->setFrequencies(freqs);
				if (model->isGamma())
					finalModel->setAlpha(model->getAlpha());
#ifdef _IG_MODELS
				if (model->isPInv())
				finalModel->setpInv(model->getpInv());
#endif
				finalModel->setLnL(model->getLnL());
				finalModel->setTree(ctree);
				models.push_back(finalModel);
			}
			free(ctree);
		}
		int bestModelIndex;
		ofs.read((char *) &bestModelIndex, (streamsize) sizeof(int));

		SelectionModel * sm = (SelectionModel *) alloca(sizeof(SelectionModel));
		ofs.read((char *) sm, sizeof(SelectionModel));
		SelectionModel * selectionmodel = new SelectionModel(
				models.at((size_t) bestModelIndex), sm->getValue());
		selectionmodel->setWeight(sm->getWeight());
		selectionmodel->setCumWeight(sm->getCumWeight());
		selectionmodel->setDelta(sm->getDelta());
		selectionmodel->setBicScore(sm->getBicScore());
		selectionmodel->setAicScore(sm->getAicScore());
		selectionmodel->setAiccScore(sm->getAiccScore());
		selectionmodel->setIndex(bestModelIndex);
		setBestModel(selectionmodel);
		delete selectionmodel;
	}

	/* branch lengths */
	size_t numBranches;
	ofs.read((char *) &numBranches, (streamsize) sizeof(size_t));
	if (numBranches) {
		if (branchLengths) {
			free(branchLengths);
		}
		branchLengths = (double *) malloc(numBranches * sizeof(double));
		ofs.read((char *) branchLengths, (streamsize) (numBranches * sizeof(double)));
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

	if (starting_topology == StartTopoFIXED && !branchLengths) {
		/* if topology is fixed for each element, keep the branch lengths */
		branchLengths = treeManager->getBranchLengths();
//		branchLengths = (double *) malloc(
//				(size_t) Utilities::numberOfBranches(num_taxa)
//						* sizeof(double));
//		for (int i = 0; i < Utilities::numberOfBranches(num_taxa); i++) {
//			if (isnan(treeManager->_tree->nodep[i + 1]->z[0])) {
//				cerr << "ERROR: NAN branch in " << i + 1 << endl;
//				exit_partest(EX_SOFTWARE);
//			}
//			branchLengths[i] = treeManager->_tree->nodep[i + 1]->z[0];
//		}
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

	size_t treeStrLen = 0;
	for (size_t i = 0; i < models.size(); i++) {
		treeStrLen += models[i]->getTree().size() + 1;
	}
	size_t hashlen = ckphash.length();
	size_t ckpSize = 6 * sizeof(size_t) + hashlen
			+ treeStrLen + models.size()
					* (modelSize + sizeof(size_t)
							+ ((size_t) (numberOfFrequencies + numberOfRates))
									* sizeof(double)) + sizeof(int)
			+ sizeof(SelectionModel)
			+ (branchLengths ?
					(size_t)Utilities::numberOfBranches((int)num_taxa) * sizeof(double) : sizeof(size_t));
	ofs.write((char *) &ckpSize, sizeof(size_t));
	ofs.write((char *) &hashlen, sizeof(size_t));
	if (hashlen > 0) {
		ofs.write((char *) ckphash.c_str(), hashlen);
	}
	ofs.write((char *) &numberOfSections, sizeof(size_t));
	ofs.write((char *) &numberOfSites, (streamsize) sizeof(size_t));
	ofs.write((char *) &numberOfPatterns, (streamsize) sizeof(size_t));
	ofs.write((char *) &number_of_models, (streamsize) sizeof(size_t));
	int bestModelIndex = 0;
	for (size_t i = 0; i < models.size(); i++) {
		Model * model = models.at(i);
		if (model == getBestModel()->getModel()) {
			bestModelIndex = (int)i;
		}
		model->setFrequencies(model->getFrequencies());
		ofs.write((char *) model, (streamsize) modelSize);

		assert(numberOfFrequencies == model->getNumberOfFrequencies());
		ofs.write((char *) model->getFrequencies(),
				(streamsize) ((size_t)model->getNumberOfFrequencies() * sizeof(double)));

		if (data_type == DT_NUCLEIC) {
			ofs.write((char *) model->getRates(),
					(streamsize) (NUM_DNA_RATES * sizeof(double)));
		}
		size_t tree_len = strlen(model->getTree().c_str()) + 1;
		ofs.write((char *) &tree_len, (streamsize) sizeof(size_t));
		ofs.write((char *) model->getTree().c_str(), (streamsize) tree_len);
	}
	/* best model */
	SelectionModel * selectionmodel = getBestModel();
	ofs.write((char *) &bestModelIndex, (streamsize) sizeof(int));
	ofs.write((char *) selectionmodel, (streamsize) sizeof(SelectionModel));

	/* branch lengths */
	if (branchLengths) {
		size_t numBranches = (size_t)Utilities::numberOfBranches((int)num_taxa);
		ofs.write((char *) &numBranches, (streamsize) sizeof(size_t));
		ofs.write((char *) branchLengths, (streamsize) (numBranches * sizeof(double)));
	} else {
		size_t zero = 0;
		ofs.write((char *) &zero, (streamsize) sizeof(size_t));
	}

	ofs.close();

	return CHECKPOINT_SAVED;
}

void PartitionElement::print(ostream & out) {
	cout << name << endl;
	if (isOptimized()) {
		cout << "Best model: " << bestModel->getModel()->getName() << endl;
		cout << "Num.Params: "
				<< bestModel->getModel()->getNumberOfFreeParameters() << endl;
		cout << "BIC score:  " << bestModel->getValue() << endl;
		cout << "BIC weight: " << bestModel->getWeight() << endl;
		bestModel->getModel()->print(out, "  ");
	}
}

}

