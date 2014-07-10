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

#include "SelectionModel.h"

namespace partest {

SelectionModel::SelectionModel(Model * model, double value) :
		model(model), value(value) {
	/* initialize to invalid values */
	index = -1;
	delta = -1.0;
	weight = -1.0;
	cumWeight = -1.0;

	bicScore = -1.0;
	aicScore = -1.0;
	aiccScore = -1.0;
	dtScore = -1.0;
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

double SelectionModel::getBicScore(void) {
	return bicScore;
}

double SelectionModel::getAicScore(void) {
	return aicScore;
}

double SelectionModel::getAiccScore(void) {
	return aiccScore;
}

double SelectionModel::getDTScore(void) {
	return dtScore;
}

void SelectionModel::setBicScore(double value) {
	bicScore = value;
}

void SelectionModel::setAicScore(double value) {
	aicScore = value;
}

void SelectionModel::setAiccScore(double value) {
	aiccScore = value;
}

void SelectionModel::setDTScore(double value) {
	dtScore = value;
}

SelectionModel * SelectionModel::clone(void) {
	SelectionModel * cloneModel = new SelectionModel(*this);
	return cloneModel;
}

} /* namespace partest */
