/*
 * SelectionModel.h
 *
 *  Created on: Jan 14, 2013
 *      Author: diego
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
