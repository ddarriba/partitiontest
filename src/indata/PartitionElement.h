/*
 * PartitionElemet.h
 *
 *  Created on: Jan 14, 2013
 *      Author: diego
 */

#ifndef PARTITION_ELEMENT_H_
#define PARTITION_ELEMENT_H_

#include "model/ModelSet.h"
#include "util/GlobalDefs.h"
#include "indata/Alignment.h"
#include "selection/SelectionModel.h"

namespace partest {

class PartitionElement {
public:
	PartitionElement(t_partitionElementId id, string name, Alignment * alignment, int start, int end, int stride, bitMask rateVariation, DataType dataType);
	PartitionElement(t_partitionElementId id, string name, Alignment * alignment, int * start, int * end, int * stride, int numberOfSections, bitMask rateVariation, DataType dataType);
	ModelSet * getModelset() { return modelset; }
	Alignment * getAlignment() { return alignment; }
	PartitionElement * splitPartition(int first, int last);
	int isOptimized();
	virtual ~PartitionElement();
	t_partitionElementId getId() { return id; }
	int getStart(int section = 0);
	int getEnd(int section = 0);
	int getStride(int section = 0);
	int getNumberOfSections(void);
	string getName(void);
	SelectionModel * getBestModel(void);
	void setBestModel(SelectionModel *);
private:
	t_partitionElementId id;
	Alignment * alignment;
	int * start;
	int * end;
	int * stride;
	int numberOfSections;
	string name;
	SelectionModel * bestModel;
	ModelSet * modelset;
};

} /* namespace partest */
#endif /* PARTITION_ELEMENT_H_ */
