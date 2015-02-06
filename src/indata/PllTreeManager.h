/*
 * PLLTreeManager.h
 *
 *  Created on: Sep 18, 2014
 *      Author: diego
 */

#include "util/GlobalDefs.h"
#include "indata/TreeManager.h"

#include <vector>

#ifndef PLLTREEMANAGER_H_
#define PLLTREEMANAGER_H_

using namespace std;

namespace partest {

class PllTreeManager: public TreeManager {
public:
	PllTreeManager(const t_partitionElementId id,
			const pllAlignmentData * phylip, const vector<PEsection> & sections,
			size_t numberOfSites);
	virtual ~PllTreeManager();

	virtual double * getBranchLengths(void);
	virtual void setModelParameters(const Model * _model, int index,
			bool setAlphaFreqs);
	virtual double searchMlTopology(bool estimateModel);
	virtual double getLikelihood();
	virtual void optimizeBranchLengths(int smoothIterations);
	virtual void optimizeModelParameters(double epsilon);
	virtual void optimizeBaseFreqs(double epsilon);
	virtual void optimizeRates(double epsilon);
	virtual void optimizeAlphas(double epsilon);
	virtual double evaluateLikelihood(bool fullTraversal);
	virtual const char * getNewickTree();

	virtual double * getFrequencies(size_t partition = 0);
	virtual double * getRates(size_t partition = 0);
	virtual double getAlpha(size_t partition = 0);

	virtual int getAutoProtModel(size_t partition = 0);

private:
	void scaleBranchLengthsSymmetric( int smoothIterations );
	double scaleBranchLengths( double multiplier );
	vector<double> storedBranchLengths;

	pllInstance * _tree;
	pllAlignmentData * _alignData;
	partitionList * _partitions;
};

}

#endif /* PLLTREEMANAGER_H_ */
