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
 * @file SearchAlgorithm.h
 * @author Diego Darriba
 * @brief Interface for search algorithms
 */

#ifndef SEARCHALGORITHM_H_
#define SEARCHALGORITHM_H_

#include "indata/PartitioningScheme.h"
#include "exe/ModelOptimize.h"

namespace partest {

class SearchAlgorithm {
public:
	SearchAlgorithm();
	virtual ~SearchAlgorithm();
	virtual PartitioningScheme * start(PartitioningScheme * startingPoint = 0) = 0;
protected:
	ModelOptimize mo;
	class SchemeManager {
	public:
		SchemeManager();
		~SchemeManager();

		int addSchemes(vector<PartitioningScheme *> schemesToAdd);
		int addScheme(PartitioningScheme * schemeToAdd);
		int optimize(ModelOptimize &mo);
	private:
		vector<PartitioningScheme *> * nextSchemes;
	};
};

} /* namespace partest */

#endif /* SEARCHALGORITHM_H_ */
