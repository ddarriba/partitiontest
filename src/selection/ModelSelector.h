/*
 * ModelSelector.h
 *
 *  Created on: Jan 14, 2013
 *      Author: diego
 */

#ifndef MODELSELECTOR_H_
#define MODELSELECTOR_H_

#include "SelectionModel.h"
#include "model/ModelSet.h"
#include "util/GlobalDefs.h"
#include "indata/PartitionElement.h"
#include <vector>

namespace partest {

class ModelSelector {
public:
	ModelSelector(ModelSet * modelset, InformationCriterion ic,
			double sampleSize);
	ModelSelector(PartitionElement * partitionElement, InformationCriterion ic,
				double sampleSize);
	virtual ~ModelSelector();
	double getAlphaImportance(void) const;
	double getInvImportance(void) const;
	double getAlphaInvImportance(void) const;
	double getFImportance(void) const;
	double getOverallAlpha(void) const;
	double getOverallInv(void) const;
	double getOverallInvAlpha(void) const;
	double getOverallAlphaInv(void) const;
	SelectionModel * getBestModel(void) {
		return bestModel;
	}
	static double computeIc(InformationCriterion ic, double lnL,
			int freeParameters, double sampleSize);
	static double computeBic(double lnL, int freeParameters, double sampleSize);
	static double computeAic(double lnL, int freeParameters);
	static double computeAicc(double lnL, int freeParameters, double sampleSize);
	void print(ostream& out);
private:
	void doSelection(ModelSet * modelset, InformationCriterion ic,
			double sampleSize);
	vector<SelectionModel *> * selectionModels;
	SelectionModel * bestModel;
	PartitionElement * partitionElement;
	InformationCriterion ic;
	double sampleSize;
	double minValue;

	double alphaImportance;
	double invImportance;
	double alphaInvImportance;
	double fImportance;
	double overallAlpha;
	double overallInv;
	double overallAlphaInv;
	double overallInvAlpha;
	DataType dataType;
};

} /* namespace partest */
#endif /* MODELSELECTOR_H_ */
