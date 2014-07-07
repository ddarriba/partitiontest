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

/**
 * @file SelectionModel.h
 *
 * @brief Wrapper of a substitution model including model selection information
 */

#ifndef SELECTIONMODEL_H_
#define SELECTIONMODEL_H_

#include "Model.h"

namespace partest {

class SelectionModel {
public:
	SelectionModel(Model * model, double value);
	int getIndex(void);
	Model * getModel(void);
	double getValue(void);
	double getWeight(void);
	double getCumWeight(void);
	double getDelta(void);
	void setIndex(int index);
	void setDelta(double delta);
	void setValue(double value);
	void setWeight(double weight);
	void setCumWeight(double cumWeight);
	virtual ~SelectionModel();
	bool operator<( const SelectionModel& other ) const;
	bool operator>( const SelectionModel& other ) const;
	SelectionModel * clone(void);
private:
	/** Selection index 0..(N-1) */
	int index;
	/** Inner substitution model */
	Model * model;
	/** Absolute criterion value */
	double value;
	/** Absolute criterion weight */
	double weight;
	/** Difference with the weight of the best-fit model */
	double delta;
	/** Cumulative weight of the models better than this */
	double cumWeight;
};

} /* namespace partest */
#endif /* SELECTIONMODEL_H_ */
