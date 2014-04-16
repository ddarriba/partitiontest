/*
 * PartitionElement.h
 *
 *  Created on: Apr 8, 2014
 *      Author: diego
 */

#ifndef PARTITIONELEMENT_H_
#define PARTITIONELEMENT_H_

#include "util/GlobalDefs.h"
#include "model/SelectionModel.h"
#include "model/NucleicModel.h"
#include "model/ProteicModel.h"
#include <string>
#include <vector>
#include <pll.h>

using namespace std;

namespace partest {

class PartitionElement {

public:
	PartitionElement(t_partitionElementId id);
	t_partitionElementId getId(void) { return id; }
	string & getName(void) { return name; }

	int setupStructures(void);
	int destroyStructures(void);
	virtual ~PartitionElement();

	int getNumberOfSites(void) const;
	int getNumberOfPatterns(void) const;
	int getNumberOfSections(void) const;

	int getNumberOfModels(void) const;
	Model * getModel(unsigned int index);
	vector<Model *> getModels(void) const;

	double getLnL(void) const;
	SelectionModel * getBestModel(void);
	void setBestModel(SelectionModel * model);

	bool isReady(void);
	bool isOptimized(void);

	pllInstance * getTree(void);
	partitionList * getPartitions(void);
	pllAlignmentData * getAlignData(void);

	double getSampleSize(void);

	void print(ostream & out);
private:

	int loadData(void);
	int storeData(void);

	bool ready = false;

	t_partitionElementId id;
	int numberOfSections;

	int numberOfSites;
	int numberOfPatterns;

	vector<Model *> models;
	SelectionModel * bestModel;

	string name, ckpname, ckphash;
	double sampleSize = 0.0;

	pllAlignmentData * _alignData = 0;
	pllInstance * _tree = 0;
	partitionList * _partitions = 0;

	bool ckpLoaded = false;
};

}
#endif /* PARTITIONELEMENT_H_ */
