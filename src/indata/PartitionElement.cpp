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
	cout << "[TRACE] PartitionElement: Creating " << id.at(0) << endl;
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
}

PartitionElement::PartitionElement(t_partitionElementId id, string name,
		Alignment * alignment, int * start, int * end, int * stride,
		int numberOfSections, bitMask rateVariation, DataType dataType) :
		id(id), name(name), numberOfSections(numberOfSections) {
	this->bestModel = 0;
	this->start = start;
	this->end = end;
	this->stride = stride;

	this->alignment = alignment->splitAlignment(start, end, numberOfSections);
	modelset = new ModelSet(rateVariation, dataType, alignment->getNumSeqs());
}

PartitionElement::~PartitionElement() {
	free(start);
	free(end);
	free(stride);

	delete (alignment);
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
