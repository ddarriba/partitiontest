/*
 * PartitionSelector.h
 *
 *  Created on: Mar 25, 2013
 *      Author: diego
 */

#ifndef PARTITIONSELECTOR_H_
#define PARTITIONSELECTOR_H_

#include "indata/PartitioningScheme.h"
#include "util/GlobalDefs.h"

namespace partest {

typedef struct {
	PartitioningScheme * scheme;
	double lnL;
	double numParameters;
	double value;
	double delta;
	double weight;
} SelectionPartitioningScheme;

class PartitionSelector {
public:
	PartitionSelector(PartitioningScheme ** schemesArray,
			int numberOfPartitions, InformationCriterion ic,
			SampleSize sampleSize, double sampleSizeValue = 0.0);
	virtual ~PartitionSelector();
	PartitioningScheme * getBestScheme(void) { return bestSelectionScheme->scheme; }
	SelectionPartitioningScheme * getBestSelectionScheme(void) { return bestSelectionScheme; }
	void print(void);
private:
	PartitioningScheme ** schemesArray;
	SelectionPartitioningScheme * bestSelectionScheme;
	int numberOfSchemes;
	InformationCriterion ic;
	SampleSize sampleSize;
	double sampleSizeValue;
	vector<SelectionPartitioningScheme *> * schemesVector;
};

} /* namespace partest */
#endif /* PARTITIONSELECTOR_H_ */
