/*
 * SelectionModel.cpp
 *
 *  Created on: Jan 14, 2013
 *      Author: diego
 */

#include "SelectionModel.h"
#include <iostream>

namespace partest {

SelectionModel::SelectionModel(Model * model, double value) :
		model(model), value(value) {
	/* initialize to invalid values */
	index = -1;
	delta = -1.0;
	weight = -1.0;
	cumWeight = -1.0;
}

SelectionModel::~SelectionModel() {
	// NOTHING
}

bool SelectionModel::operator<(const SelectionModel& other) const {
	return weight < other.weight;
}

bool SelectionModel::operator>(const SelectionModel& other) const {
	return weight > other.weight;
}

void SelectionModel::setIndex(int index) {
	this->index = index;
}

void SelectionModel::setWeight(double weight) {
	this->weight = weight;
}

void SelectionModel::setCumWeight(double cumWeight) {
	this->cumWeight = cumWeight;
}

void SelectionModel::setValue(double delta) {
	this->value = value;
}

void SelectionModel::setDelta(double delta) {
	this->delta = delta;
}

int SelectionModel::getIndex() {
	return index;
}

Model * SelectionModel::getModel() {
	return model;
}

double SelectionModel::getValue() {
	return value;
}

double SelectionModel::getWeight() {
	return weight;
}

double SelectionModel::getCumWeight() {
	return cumWeight;
}

double SelectionModel::getDelta() {
	return delta;
}

SelectionModel * SelectionModel::clone(void) {
	SelectionModel * cloneModel = new SelectionModel(*this);
	return cloneModel;
}

} /* namespace partest */
