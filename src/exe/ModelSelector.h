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
 * @file ModelSelector.h
 *
 * @brief Perform selection among a set of models
 */

#ifndef ModelSelector_H_
#define ModelSelector_H_

#include "model/SelectionModel.h"
#include "util/GlobalDefs.h"
#include "indata/PartitionElement.h"
#include <vector>

namespace partest {

using namespace std;

class ModelSelector {
public:
	ModelSelector(PartitionElement * partitionElement, InformationCriterion ic,
			double sampleSize);
	virtual ~ModelSelector();
	double getAlphaImportance(void) const;
	double getFImportance(void) const;
	double getOverallAlpha(void) const;
#ifdef _IG_MODELS
	double getInvImportance(void) const;
	double getAlphaInvImportance(void) const;
	double getOverallInv(void) const;
	double getOverallInvAlpha(void) const;
	double getOverallAlphaInv(void) const;
#endif
	SelectionModel * getBestModel(void) {
		return bestModel;
	}
	static double computeIc(InformationCriterion ic, double lnL,
			int freeParameters, double sampleSize);
	static double computeBic(double lnL, int freeParameters, double sampleSize);
	static double computeAic(double lnL, int freeParameters);
	static double computeAicc(double lnL, int freeParameters,
			double sampleSize);
	void print(ostream& out);
private:
	void doSelection(vector<Model *> modelset, InformationCriterion ic,
			double sampleSize);
	vector<SelectionModel *> * selectionModels;
	SelectionModel * bestModel;
	PartitionElement * partitionElement;
	InformationCriterion ic;
	double sampleSize;
	double minValue;

	double alphaImportance;
	double fImportance;
	double overallAlpha;
#ifdef _IG_MODELS
	double invImportance;
	double alphaInvImportance;
	double overallInv;
	double overallAlphaInv;
	double overallInvAlpha;
#endif

};

} /* namespace partest */
#endif /* ModelSelector_H_ */
