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
 * @file PartitionSelector.h
 *
 * @brief Perform selection among partitioning schemes
 */

#ifndef PARTITIONSELECTOR_H_
#define PARTITIONSELECTOR_H_

#include "indata/PartitioningScheme.h"
#include "util/GlobalDefs.h"
#include <iostream>

namespace partest {

typedef struct {
	PartitioningScheme * scheme;
	double lnL;
	double numParameters;
	double value;
	double delta;
	double weight;
	void print(ostream & out);
} SelectionPartitioningScheme;

class PartitionSelector {
public:
	PartitionSelector(vector<PartitioningScheme *> schemesArray);
	virtual ~PartitionSelector();
	PartitioningScheme * getBestScheme(void) {
		return bestSelectionScheme->scheme;
	}
	SelectionPartitioningScheme * getBestSelectionScheme(void) {
		return bestSelectionScheme;
	}
	void print(ostream& out, int limit=-1);
private:
	vector<PartitioningScheme *> schemesArray;
	SelectionPartitioningScheme * bestSelectionScheme;
	vector<SelectionPartitioningScheme *> * schemesVector;
};

} /* namespace partest */
#endif /* PARTITIONSELECTOR_H_ */
