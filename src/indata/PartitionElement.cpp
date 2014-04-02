/*
 * PartitionElemet.cpp
 *
 *  Created on: Jan 14, 2013
 *      Author: diego
 */

#include "PartitionElement.h"
#include <algorithm>
#include <fstream>
#ifdef _PLL
#include "PLLAlignment.h"
#endif

#define CHECKPOINT_LOADED 0
#define CHECKPOINT_UNAVAILABLE 1
#define CHECKPOINT_UNEXISTENT 2

namespace partest {

int PartitionElement::loadData() {
	if (!ckpAvailable) return CHECKPOINT_UNAVAILABLE;

	fstream ofs ((ckpPath + os_separator + ckpname).c_str(), ios::in);

	if (!ofs) return CHECKPOINT_UNEXISTENT;

	alignment = 0;

	ofs.seekg (0);
	int numTaxa;
	ofs.read((char *) &(numTaxa), sizeof(int));
	int treeLen;
	ofs.read((char *) &(treeLen), sizeof(int));
	char treeStr[treeLen];
	ofs.read((char *) treeStr, treeLen);
	ModelSet * ms = (ModelSet *) alloca(sizeof(ModelSet));
	if (!ms) {
		cerr << "[INTERNAL_ERROR] Error allocating memory for modelset" << endl;
		Utilities::exit_partest(EX_SOFTWARE);
	}
	ofs.read((char *) ms, sizeof(ModelSet));

	modelset = new ModelSet(rateVariation, ms->getDataType(), numTaxa, optimizeMode, false, true);
	ms->allocateModels(ms->getNumberOfModels());
	size_t modelSize = ms->getDataType()==DT_NUCLEIC?sizeof(NucleicModel):sizeof(ProteicModel);
	for (unsigned int i = 0; i < ms->getNumberOfModels(); i++) {
			Model * model = 0;
			switch (ms->getDataType()) {
			case DT_NUCLEIC:
			  model = (NucleicModel *) alloca(modelSize);
			  break;
			case DT_PROTEIC:
			  model = (ProteicModel *) alloca(modelSize);
			  break;
			default:
				cerr << "[INTERNAL_ERROR] Unrecognized datatype" << endl;
				Utilities::exit_partest(EX_SOFTWARE);
			}
			ofs.read((char *) model, modelSize);
			double * freqs = (double *) alloca(model->getNumberOfFrequencies() * sizeof(double));
			ofs.read((char *) freqs, model->getNumberOfFrequencies() * sizeof(double));

			//model->allocateFrequencies(model->getNumberOfFrequencies());
			//model->setFrequencies(freqs);

			double * rates = 0;
			if (ms->getDataType()==DT_NUCLEIC) {
				rates = (double *) alloca(NUM_RATES * sizeof(double));
				ofs.read((char *) rates, NUM_RATES * sizeof(double));
				//model->setRates(rates);
			}
			int len_tree;
			ofs.read((char *) &len_tree, sizeof(size_t));
			char ctree[len_tree];
			ofs.read((char *) &ctree, len_tree);

			if (ms->getDataType() == DT_NUCLEIC) {
				NucleicModel * finalModel = new NucleicModel(((NucleicModel *)model)->getMatrix(), model->getRateVariation(), numTaxa);
				finalModel->setFrequencies(freqs);
				finalModel->setRates(rates);
				if (model->isGamma())
					finalModel->setAlpha(model->getAlpha());
				if (model->isPInv())
					finalModel->setpInv(model->getpInv());
				finalModel->setLnL(model->getLnL());
				finalModel->setTree(ctree);
				modelset->setModel(finalModel, i);
			} else {
				ProteicModel * finalModel = new ProteicModel(((ProteicModel *)model)->getMatrix(), model->getRateVariation(), numTaxa);
				finalModel->setFrequencies(freqs);
				if (model->isGamma())
					finalModel->setAlpha(model->getAlpha());
				if (model->isPInv())
					finalModel->setpInv(model->getpInv());
				finalModel->setLnL(model->getLnL());
				finalModel->setTree(ctree);
				modelset->setModel(finalModel, i);
			}
	}
	int bestModelIndex;
	ofs.read((char *) &bestModelIndex, sizeof(int));
	SelectionModel * sm = (SelectionModel *) alloca(sizeof(SelectionModel));
	ofs.read((char *) sm, sizeof(SelectionModel));
	SelectionModel * selectionmodel = new SelectionModel(modelset->getModel(bestModelIndex), sm->getValue());
	selectionmodel->setWeight(sm->getWeight());
	selectionmodel->setCumWeight(sm->getCumWeight());
	selectionmodel->setDelta(sm->getDelta());
	selectionmodel->setIndex(bestModelIndex);
	setBestModel(selectionmodel, false);
	ofs.close();
	return CHECKPOINT_LOADED;
}

int PartitionElement::storeData() {

	if (!ckpAvailable) return -1;
	if (!isOptimized()) {
		cerr << "[INTERNAL_ERROR] Attempting to save unoptimized Partition Element" << endl;
		Utilities::exit_partest(EX_SOFTWARE);
	}

	fstream ofs ((ckpPath + os_separator + ckpname).c_str(), ios::out);

	//ofs.open((ckpPath + os_separator + ckpname).c_str());
	ofs.seekg (0);
	PLLAlignment * align = ((PLLAlignment *) alignment);
	pllInstance * tree = align->getTree();
	ofs.write((char *) &(tree->mxtips), sizeof(int));
	ofs.write((char *) &(tree->treeStringLength), sizeof(int));
	ofs.write((char *) tree->tree_string, tree->treeStringLength);
	ofs.write((char *) modelset, sizeof(ModelSet));
	size_t modelSize = modelset->getDataType()==DT_NUCLEIC?sizeof(NucleicModel):sizeof(ProteicModel);
	int bestModelIndex = 0;
	for (unsigned int i = 0; i < modelset->getNumberOfModels(); i++) {
		Model * model = modelset->getModel(i);
		if (model == getBestModel()->getModel()) {
			bestModelIndex = i;
		}
		model->setFrequencies(model->getFrequencies());
		ofs.write((char *) model, modelSize);
		ofs.write((char *) model->getFrequencies(), model->getNumberOfFrequencies() * sizeof(double));
		if (modelset->getDataType()==DT_NUCLEIC)
			ofs.write((char *) model->getRates(), NUM_RATES * sizeof(double));
		int name_len, matrixname_len, tree_len;
		name_len = model->getName().length() + 1;
		matrixname_len = model->getMatrixName().length() + 1;
		tree_len = model->getTree().length() + 1;
		ofs.write((char *) &tree_len, sizeof(size_t));
		ofs.write((char *) model->getTree().c_str(), tree_len);
	}
	/* best model */
	SelectionModel * selectionmodel = getBestModel();
	ofs.write((char *) &bestModelIndex, sizeof(int));
	ofs.write((char *) selectionmodel, sizeof(SelectionModel));
	ofs.close();

	delete alignment;
	alignment = 0;

	return 0;
}

PartitionElement::PartitionElement(t_partitionElementId id, string name,
		Alignment * alignment, int start, int end, int stride,
		bitMask rateVariation, DataType dataType, OptimizeMode optimizeMode) :
		id(id), name(name), ckpname(name), optimizeMode(optimizeMode), rateVariation(rateVariation) {

#ifdef DEBUG
	cout << "[TRACE] PartitionElement: Creating " << name << endl;
#endif
	numberOfSections = 1;
	this->bestModel = 0;
	this->start = (int *) malloc(sizeof(int));
	this->end = (int *) malloc(sizeof(int));
	this->stride = (int *) malloc(sizeof(int));
	this->start[0] = start;
	this->end[0] = end;
	this->stride[0] = stride;

#ifdef _PLL
	partitionInfo = 0;
#endif
	ckpname.erase(remove(ckpname.begin(), ckpname.end(), '('), ckpname.end());
	ckpname.erase(remove(ckpname.begin(), ckpname.end(), ')'), ckpname.end());
	replace( ckpname.begin(), ckpname.end(), ' ', '_');

	if (loadData() != CHECKPOINT_LOADED) {
		modelset = new ModelSet(rateVariation, dataType, alignment->getNumSeqs(), optimizeMode);
		this->alignment = alignment->splitAlignment(start, end);
	}
}

PartitionElement::PartitionElement(t_partitionElementId id, string name,
		Alignment * alignment, int * start, int * end, int * stride,
		int numberOfSections, bitMask rateVariation, DataType dataType, OptimizeMode optimizeMode) :
		id(id), name(name), ckpname(name), numberOfSections(numberOfSections), optimizeMode(optimizeMode),
		rateVariation(rateVariation) {

#ifdef DEBUG
	cout << "[TRACE] PartitionElement: Creating " << name <<  "  Sections: " << numberOfSections << endl;
#endif
	this->start = (int *) malloc(numberOfSections * sizeof(int));
	this->end = (int *) malloc(numberOfSections * sizeof(int));
	this->stride = (int *) malloc(numberOfSections * sizeof(int));

	for (int i=0; i< numberOfSections; i++) {
#ifdef DEBUG
	cout << "[TRACE] PartitionElement: Section " << i+1 << "/" << numberOfSections << ":    (" << start[i] << "-" << end[i] << "\\" << stride[i] << ")" << endl;
#endif
		this->start[i] = start[i];
		this->end[i] = end[i];
		this->stride[i] = stride[i];
	}

	this->bestModel = 0;

#ifdef DEBUG
	cout << "[TRACE] PartitionElement: Done" << endl;
#endif

#ifdef _PLL
	partitionInfo = 0;
#endif
	ckpname.erase(remove(ckpname.begin(), ckpname.end(), '('), ckpname.end());
	ckpname.erase(remove(ckpname.begin(), ckpname.end(), ')'), ckpname.end());
	replace( ckpname.begin(), ckpname.end(), ' ', '_');

#ifdef DEBUG
	cout << "[TRACE] PartitionElement: Creating modelset" << endl;
#endif
	if (loadData() != CHECKPOINT_LOADED) {
		modelset = new ModelSet(rateVariation, dataType, alignment->getNumSeqs(), optimizeMode);
#ifdef DEBUG
	cout << "[TRACE] PartitionElement: Creating alignment" << endl;
#endif
		this->alignment = alignment->splitAlignment(start, end, numberOfSections);
	}
}

PartitionElement::~PartitionElement() {
	free(start);
	free(end);
	free(stride);

	if (alignment)
		delete alignment;
	if (modelset)
		delete modelset;

	if (bestModel != 0)
		delete bestModel;
}

SelectionModel * PartitionElement::getBestModel(void) {
	return bestModel;
}

void PartitionElement::setBestModel(SelectionModel * bestModel, bool save) {
	if (this->bestModel) {
		delete this->bestModel;
	}
	this->bestModel = bestModel->clone();
	if (save)
		storeData();
}

bool PartitionElement::isOptimized() {
	if (bestModel) return true;
	bool optimized = true;
	for (unsigned int i = 0; i < modelset->getNumberOfModels(); i++) {
		optimized &= modelset->getModel(i)->isOptimized();
	}
	return optimized;
}

void PartitionElement::buildCompleteModelSet(bool clearAll) {
	modelset->buildCompleteModelSet(clearAll);
	bestModel = 0;
}

int PartitionElement::getStart(int section) {
	return start[section];
}

int PartitionElement::getEnd(int section) {
	return end[section];
}

int PartitionElement::getStride(int section) {
	return stride[section];
}

int PartitionElement::getNumberOfSections() {
	return numberOfSections;
}

string PartitionElement::getName() {
	return name;
}

} /* namespace partest */
