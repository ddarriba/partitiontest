/*
 * PartitionElemet.cpp
 *
 *  Created on: Jan 14, 2013
 *      Author: diego
 */

#include "PartitionElement.h"

namespace partest {

PartitionElement::PartitionElement(t_partitionElementId id, string name,
		Alignment * alignment, int start, int end, int stride,
		bitMask rateVariation, DataType dataType) :
		id(id), name(name) {

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
	this->alignment = alignment->splitAlignment(start, end);
	modelset = new ModelSet(rateVariation, dataType, alignment->getNumSeqs());
#ifdef _PLL
	partitionInfo = 0;
#endif
}

PartitionElement::PartitionElement(t_partitionElementId id, string name,
		Alignment * alignment, int * start, int * end, int * stride,
		int numberOfSections, bitMask rateVariation, DataType dataType) :
		id(id), name(name), numberOfSections(numberOfSections) {

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
	cout << "[TRACE] PartitionElement: Creating alignment" << endl;
#endif
	this->alignment = alignment->splitAlignment(start, end, numberOfSections);

#ifdef DEBUG
	cout << "[TRACE] PartitionElement: Creating modelset" << endl;
#endif
	modelset = new ModelSet(rateVariation, dataType, alignment->getNumSeqs());
#ifdef DEBUG
	cout << "[TRACE] PartitionElement: Done" << endl;
#endif

#ifdef _PLL
	partitionInfo = 0;
#endif
}

PartitionElement::~PartitionElement() {
	free(start);
	free(end);
	free(stride);

	delete alignment;
	delete modelset;

	if (bestModel != 0)
		delete bestModel;
}

SelectionModel * PartitionElement::getBestModel(void) {
	return bestModel;
}

void PartitionElement::setBestModel(SelectionModel * bestModel) {
	if (this->bestModel) {
		delete this->bestModel;
	}
	this->bestModel = bestModel->clone();
}

bool PartitionElement::isOptimized() {
	if (bestModel) return true;
	bool optimized = true;
	for (int i = 0; i < modelset->getNumberOfModels(); i++) {
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
